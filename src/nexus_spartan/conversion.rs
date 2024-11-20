#![allow(warnings)]

use ark_ff::One;
use ark_ff::Zero;
use ark_std::{end_timer, start_timer};
use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::R1CSShape;
use crate::gadgets::sparse::SparseMatrix;
use crate::nexus_spartan::analyze_vector_sparseness;
use crate::nexus_spartan::crr1cs::{CRR1CSInstance, CRR1CSShape, CRR1CSWitness};
use crate::nexus_spartan::errors::R1CSError;
use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
use crate::nexus_spartan::{Assignment, Instance};
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_relations::r1cs::ConstraintSystemRef;
use std::error::Error;
use std::fmt::Display;

#[derive(Debug)]
pub enum ConversionError {
    ConversionError(R1CSError),
}

impl From<R1CSError> for ConversionError {
    fn from(error: R1CSError) -> Self {
        Self::ConversionError(error)
    }
}

impl Error for ConversionError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ConversionError::ConversionError(e) => Some(e),
        }
    }
}

impl Display for ConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ConversionError(e) => write!(f, "Conversion error: {e}"),
        }
    }
}

impl<F: PrimeField + Absorb> CRR1CSShape<F> {
    pub(crate) fn convert<G: SWCurveConfig<ScalarField=F>>(cs: ConstraintSystemRef<F>) -> Self {
        // convert constraint system into R1CS shape
        let shape: R1CSShape<G> = R1CSShape::from(cs);
        // extract R1CS field
        let R1CSShape {
            num_constraints,
            num_vars,
            // This includes the leading `u` entry
            num_io,
            A,
            B,
            C,
        } = shape;
        // Spartan arranges the R1CS matrices using Z = [w, u, x], rather than [u, x, w]
        let rearrange =
            |matrix: SparseMatrix<F>| -> Vec<(usize, usize, F)> {
                matrix.iter().map(|(row, col, val)|
                    // this is a witness entry
                    if col >= num_io {
                        (row, col - num_io, val)
                    } else {
                        // this is an IO entry
                        (row, col + num_vars, val)
                    }).collect()
            };
        CRR1CSShape {
            inst: Instance::new(
                num_constraints,
                num_vars,
                // Spartan does not include the leading `u` entry in `num_inputs`.
                num_io - 1,
                rearrange(A).as_slice(),
                rearrange(B).as_slice(),
                rearrange(C).as_slice(),
            ).unwrap(),
        }
    }
}

impl<PC, E> CRR1CSInstance<E, PC>
where
    E: Pairing,
    PC: PolyCommitmentScheme<E>,
    <E as Pairing>::ScalarField: Absorb,
{
    pub(crate) fn convert<G: SWCurveConfig>(
        cs: ConstraintSystemRef<G::ScalarField>,
        key: &PC::PolyCommitmentKey,
    ) -> Self where
        E: Pairing<G1Affine=Affine<G>, ScalarField=G::ScalarField>,
    {
        let cs_borrow = cs.borrow().unwrap();
        let witness = cs_borrow.witness_assignment.clone();
        let pub_io = cs_borrow.instance_assignment.clone();

        assert!(!pub_io.is_empty(), "instance is empty");

        analyze_vector_sparseness("witness", &witness);

        let poly_W = MultilinearPolynomial::new(witness);
        let commit_timer = start_timer!(|| "Instance conversion (commit to witness)");
        let comm_W = PC::commit(&poly_W, &key);
        end_timer!(commit_timer);

        CRR1CSInstance {
            input: Assignment::new(&pub_io[1..]).unwrap(),
            comm_W,
        }
    }
}

impl<F: PrimeField + Absorb> CRR1CSWitness<F> {
    pub(crate) fn convert(cs: ConstraintSystemRef<F>) -> Self {
        let cs_borrow = cs.borrow().unwrap();
        let mut witness = cs_borrow.witness_assignment.clone();

        // Calculate the next power of two for the length of the witness
        let next_power_of_two = witness.len().next_power_of_two();

        // Pad the witness vector with F::zero() to reach the next power of two
        witness.resize(next_power_of_two, F::zero());

        CRR1CSWitness { W: Assignment::new(&witness).unwrap() }
    }
}


#[cfg(test)]
pub mod tests {
    use crate::constant_for_curves::{ScalarField, E, G1};
    use crate::hash::pederson::PedersenCommitment;
    use crate::nexus_spartan::crr1cs::{is_sat, CRR1CSInstance, CRR1CSKey, CRR1CSShape, CRR1CSWitness};
    use crate::nexus_spartan::crr1csproof::CRR1CSProof;
    use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
    use crate::pcs::kzh2::KZH2SRS;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use crate::transcript::transcript::Transcript;
    use ark_ec::short_weierstrass::Projective;
    use ark_ec::CurveConfig;
    use ark_ff::PrimeField;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::boolean::Boolean;
    use ark_r1cs_std::eq::EqGadget;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError, SynthesisMode};
    use rand::thread_rng;

    type Pedersen = PedersenCommitment<Projective<G1>>;

    // Define a simple constraint system struct
    // for the circuit: a + b == a^2
    pub struct TrivialCircuit<F: PrimeField> {
        pub a: F,  // Public input
        pub b: F,  // Public input
    }

    impl<F: PrimeField> ConstraintSynthesizer<F> for TrivialCircuit<F> {
        fn generate_constraints(
            self,
            cs: ConstraintSystemRef<F>,
        ) -> Result<(), SynthesisError> {
            // Allocate the public inputs as FpVar
            let a_var = FpVar::new_input(cs.clone(), || Ok(self.a))?;
            let b_var = FpVar::new_input(cs.clone(), || Ok(self.b))?;

            // Perform arithmetic operations
            let sum = &a_var + &b_var; // Sum a and b

            // Add a constraint: a + b == a^2
            let a_squared = &a_var * &a_var;
            let a_squared = &a_var * &a_var;
            let a_squared = &a_var * &a_var;
            let a_squared = &a_var * &a_var;
            for i in 0..12 {
                let _ = Boolean::new_witness(cs.clone(), || Ok(false));
            }

            sum.enforce_equal(&a_squared).expect("OMG");
            sum.enforce_equal(&a_squared).expect("OMG");
            sum.enforce_equal(&a_squared).expect("OMG");
            sum.enforce_equal(&a_squared).expect("OMG");
            sum.enforce_equal(&a_squared).expect("OMG");

            Ok(())
        }
    }


    #[test]
    fn test_r1cs_conversion() {
        // Create a new constraint system for: a + b == a^2
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // Example public inputs
        // 4 + 12 == 16 == 4^2
        let a = ScalarField::from(4u32);
        let b = ScalarField::from(12u32);

        // Instantiate the trivial circuit with inputs
        let circuit = TrivialCircuit { a, b };

        // Generate the constraints
        circuit.generate_constraints(cs.clone()).unwrap();
        assert!(cs.is_satisfied().unwrap());
        println!("{}", cs.num_constraints());

        cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
        cs.finalize();

        // convert to the corresponding Spartan types
        let shape = CRR1CSShape::<ScalarField>::convert::<G1>(cs.clone());
        let SRS: KZH2SRS<E> = MultilinearPolynomial::setup(4, &mut thread_rng()).unwrap();
        let key: CRR1CSKey<E, MultilinearPolynomial<ScalarField>> = CRR1CSKey::new(&SRS, shape.get_num_cons(), shape.get_num_vars());
        let instance: CRR1CSInstance<E, MultilinearPolynomial<ScalarField>> = CRR1CSInstance::convert(cs.clone(), &key.keys.ck);
        let witness = CRR1CSWitness::<ScalarField>::convert(cs.clone());
        // check that the Spartan instance-witness pair is still satisfying
        assert!(is_sat(&shape, &instance, &witness, &key).unwrap());

        let (num_cons, num_vars, _num_inputs) = (
            shape.get_num_cons(),
            shape.get_num_vars(),
            shape.get_num_inputs(),
        );

        let mut prover_transcript = Transcript::new(b"example");

        let (proof, rx, ry) = CRR1CSProof::prove(
            &shape,
            &instance,
            witness,
            &key,
            &mut prover_transcript,
        );

        let inst_evals = shape.inst.inst.evaluate(&rx, &ry);

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(
                num_vars,
                num_cons,
                &instance,
                &inst_evals,
                &mut verifier_transcript,
            )
            .is_ok());
    }
}

#![allow(warnings)]

use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::{R1CSInstance, R1CSShape, R1CSWitness};
use crate::gadgets::sparse::SparseMatrix;
use crate::nexus_spartan::crr1cs::{CRR1CSInstance, CRR1CSKey, CRR1CSShape, CRR1CSWitness};
use crate::nexus_spartan::errors::R1CSError;
use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
use crate::nexus_spartan::{Assignment, Instance};
use crate::polynomial::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
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

impl<F: PrimeField> CRR1CSShape<F> {
    pub(crate) fn convert<G: SWCurveConfig<ScalarField=F>>(cs: ConstraintSystemRef<F>) -> Self{
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

        let poly_W = MultilinearPolynomial::new(witness);
        let comm_W = PC::commit(&poly_W, &key);

        CRR1CSInstance {
            input: Assignment::new(&pub_io[1..]).unwrap(),
            comm_W,
        }
    }
}

impl<F: PrimeField> CRR1CSWitness<F> {
    pub(crate) fn convert(cs: ConstraintSystemRef<F>, ) -> Self {
        let cs_borrow = cs.borrow().unwrap();
        let witness = cs_borrow.witness_assignment.clone();
        CRR1CSWitness { W: Assignment::new(&witness).unwrap() }
    }
}


#[cfg(test)]
mod tests {
    use crate::commitment::CommitmentScheme;
    use crate::constant_for_curves::{ScalarField, E, G1};
    use crate::gadgets::r1cs::r1cs::{R1CSInstance, R1CSWitness};
    use crate::gadgets::r1cs::R1CSShape;
    use crate::hash::pederson::PedersenCommitment;
    use crate::nexus_spartan::crr1cs::{is_sat, CRR1CSInstance, CRR1CSKey, CRR1CSShape, CRR1CSWitness};
    use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
    use crate::pcs::multilinear_pcs::SRS;
    use crate::polynomial::multilinear_poly::MultilinearPolynomial;
    use ark_ec::short_weierstrass::Projective;
    use ark_ec::CurveConfig;
    use ark_ff::PrimeField;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError};
    use rand::thread_rng;

    type Pedersen = PedersenCommitment<Projective<G1>>;

    // Define a simple constraint system struct
    struct TrivialCircuit<F: PrimeField> {
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
            for i in 0..10 {
                let product = &a_var * &b_var; // Multiply a and b
            }

            // Add a constraint: sum + product == a^2 + b^2
            let a_squared = &a_var * &a_var;
            let b_squared = &b_var * &b_var;
            let _sum_of_squares = a_squared + b_squared;

            Ok(())
        }
    }


    #[test]
    fn test_conversion() {
        // Create a new constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // Example public inputs
        let a = ScalarField::from(3u32);
        let b = ScalarField::from(4u32);

        // Instantiate the trivial circuit with inputs
        let circuit = TrivialCircuit { a, b };

        // Generate the constraints
        circuit.generate_constraints(cs.clone()).unwrap();
        assert!(cs.is_satisfied().unwrap());
        println!("{}", cs.num_constraints());



        // convert to the corresponding Spartan types
        let shape = CRR1CSShape::<ScalarField>::convert::<G1>(cs.clone());
        let SRS: SRS<E> = MultilinearPolynomial::setup(5, &mut thread_rng()).unwrap();
        let key: CRR1CSKey<E, MultilinearPolynomial<ScalarField>> = CRR1CSKey::new(&SRS, shape.get_num_cons(), shape.get_num_vars());
        let instance: CRR1CSInstance<E, MultilinearPolynomial<ScalarField>> = CRR1CSInstance::convert(cs.clone(), &key.keys.ck);
        let witness = CRR1CSWitness::<ScalarField>::convert(cs.clone());

        // check that the Spartan instance-witness pair is still satisfying
        assert!(is_sat(&shape, &instance, &witness, &key).unwrap());
    }
}

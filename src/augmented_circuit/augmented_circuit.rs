use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_ec::{AffineRepr, CurveConfig};
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

use crate::math::Math;
use crate::commitment::CommitmentScheme;
use crate::nexus_spartan::partial_verifier::partial_verifier_var::PartialVerifierVar;
use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
use crate::nexus_spartan::sparse_polynomial::sparse_polynomial::SparsePoly;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::Transcript;

use crate::nexus_spartan::partial_verifier::partial_verifier::PartialVerifier;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AugmentedCircuit<F: PrimeField + Absorb>
{
    spartan_partial_verifier: PartialVerifier<F>,
}

pub struct AugmentedCircuitVar<F: PrimeField + Absorb>
{
    spartan_partial_verifier: PartialVerifierVar<F>,
}

impl<F: PrimeField + Absorb> AllocVar<AugmentedCircuit<F>, F> for AugmentedCircuitVar<F> {
    fn new_variable<T: Borrow<AugmentedCircuit<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        // Convert to Namespace<F>
        let ns = cs.into();
        // Get the constraint system reference
        let cs = ns.cs();

        // Fetch the instance of `AugmentedCircuit<F>`
        let binding = f()?;
        let data = binding.borrow();

        // Allocate the sumcheck proof phase 1
        let partial_verifier = PartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.spartan_partial_verifier),
            mode,
        )?;

        Ok(AugmentedCircuitVar {
            spartan_partial_verifier: partial_verifier
        })
    }
}

#[cfg(test)]
mod tests {
    use ark_relations::r1cs::ConstraintSystem;

    use super::*;
    use crate::nexus_spartan::crr1cs::is_sat;
    use crate::nexus_spartan::crr1cs::produce_synthetic_crr1cs;
    use crate::nexus_spartan::crr1csproof::CRR1CSProof;
    use crate::constant_for_curves::{ScalarField, E};

    pub fn get_test_proof<E, PC, F>() -> (CRR1CSProof<E, PC, F>, PartialVerifier<F>)
    where
        F: PrimeField + Absorb,
        PC: PolyCommitmentScheme<E>,
        E: Pairing<ScalarField=F>,
    {
        let num_vars = 1024;
        let num_cons = num_vars;
        let num_inputs = 10;
        let (shape, instance, witness, gens) = produce_synthetic_crr1cs::<E, PC>(num_cons, num_vars, num_inputs);
        assert!(is_sat(&shape, &instance, &witness, &gens.gens_r1cs_sat).unwrap());

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
            &gens.gens_r1cs_sat,
            &mut prover_transcript,
        );

        let inst_evals = shape.inst.inst.evaluate(&rx, &ry);

        let mut verifier_transcript = Transcript::new(b"example");
        let mut verifier_transcript_clone1 = verifier_transcript.clone();
        let partial_verifier = PartialVerifier::initialise(
            &proof,
            num_vars,
            num_cons,
            instance.input.assignment,
            &inst_evals,
            &mut verifier_transcript,
        );

        partial_verifier.verify::<E>(&mut verifier_transcript_clone1);

        (proof, partial_verifier)
    }

    #[test]
    pub fn test_augmented_circuit() {
        let (proof, partial_verifier) = get_test_proof::<E, MultilinearPolynomial<ScalarField>, ScalarField>();


        let cs = ConstraintSystem::<ScalarField>::new_ref();
        let augmented_circuit = AugmentedCircuit {
            spartan_partial_verifier: partial_verifier.clone()
        };

        let augmented_circuit_var = AugmentedCircuitVar::new_variable(
            cs.clone(),
            || Ok(augmented_circuit.clone()),
            AllocationMode::Input,
        ).unwrap();

        let kzh_opening_proof = proof.eval_vars_at_ry;
        // Turn it into an accumulator
        


        // assert_eq!(augmented_circuit, augmented_circuit_var.value().unwrap());

        // let (_r_x, _r_y) = partial_verifier_var.verify(&mut transcript);
        // println!("constraint count: {} {}", cs.num_instance_variables(), cs.num_witness_variables());
        // assert!(cs.is_satisfied().unwrap());
    }
}




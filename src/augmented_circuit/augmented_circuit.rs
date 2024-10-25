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
    use crate::constant_for_curves::{ScalarField, E};
    use crate::nexus_spartan::partial_verifier::partial_verifier::tests::partial_verifier_test_helper;

    #[test]
    pub fn test_augmented_circuit() {
        let (partial_verifier, _transcript) = partial_verifier_test_helper::<E, MultilinearPolynomial<ScalarField>, ScalarField>();
        let cs = ConstraintSystem::<ScalarField>::new_ref();
        let augmented_circuit = AugmentedCircuit {
            spartan_partial_verifier: partial_verifier.clone()
        };

        let augmented_circuit_var = AugmentedCircuitVar::new_variable(
            cs.clone(),
            || Ok(augmented_circuit.clone()),
            AllocationMode::Input,
        ).unwrap();

        // assert_eq!(augmented_circuit, augmented_circuit_var.value().unwrap());

        // let (_r_x, _r_y) = partial_verifier_var.verify(&mut transcript);
        // println!("constraint count: {} {}", cs.num_instance_variables(), cs.num_witness_variables());
        // assert!(cs.is_satisfied().unwrap());
    }
}




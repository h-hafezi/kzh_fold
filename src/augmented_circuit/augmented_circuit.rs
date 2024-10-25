use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_ec::{AffineRepr, CurveConfig};
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

use crate::accumulation_circuit::verifier_circuit::AccumulatorVerifier;
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
pub struct AugmentedCircuit<G1, G2, C2, E>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField>,
{
    spartan_partial_verifier: PartialVerifier<G1::ScalarField>,
    kzh_acc_verifier: AccumulatorVerifier<G1, G2, C2, E>,
}

pub struct AugmentedCircuitVar<F: PrimeField + Absorb>
{
    spartan_partial_verifier: PartialVerifierVar<F>,
    // Hossein: add AccumulatorVerifierVar here
}

impl<G1, G2, C2, E> AllocVar<AugmentedCircuit<G1, G2, C2, E>, G1::ScalarField> for AugmentedCircuitVar<G1::ScalarField>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField>,
{
    fn new_variable<T: Borrow<AugmentedCircuit<G1, G2, C2, E>>>(
        cs: impl Into<Namespace<G1::ScalarField>>,
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

        // Allocate the Spartan partial verifier
        let partial_verifier = PartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.spartan_partial_verifier),
            mode,
        )?;

        // Allocate the KZH AccVerifier
        // Hossein: AccumulatorVerifierVar

        Ok(AugmentedCircuitVar {
            spartan_partial_verifier: partial_verifier
        })
    }
}

#[cfg(test)]
mod tests {
    use ark_relations::r1cs::ConstraintSystem;
    use rand::thread_rng;

    use super::*;
    use crate::accumulation::accumulator::Accumulator;
    use crate::hash::pederson::PedersenCommitment;
    use crate::nexus_spartan::crr1cs::is_sat;
    use crate::nexus_spartan::crr1cs::produce_synthetic_crr1cs;
    use crate::nexus_spartan::crr1csproof::CRR1CSProof;
    use crate::pcs::multilinear_pcs::PolyCommit;
    use crate::pcs::multilinear_pcs::SRS;
    use crate::constant_for_curves::{ScalarField, E, G1, G2};
    type C2 = PedersenCommitment<Projective<G2>>;

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

    /// Take as input `proof_i` and `running_accumulator_{i}` and produce `proof_{i+1}` and `running_accumulator_{i+1}`.
    #[test]
    pub fn test_augmented_circuit() {
        let degree_x = 512;
        let degree_y = 512;
        let srs_pcs: SRS<E> = PolyCommit::<E>::setup(degree_x, degree_y, &mut thread_rng());
        let srs = Accumulator::setup(srs_pcs.clone(), &mut thread_rng());

        // Get `running_accumulator_{i}`
        let running_accumulator = Accumulator::random_satisfying_accumulator(&srs, &mut thread_rng());

        // Get a dummy `proof_i`
        let (proof, partial_verifier) = get_test_proof::<E, MultilinearPolynomial<ScalarField>, ScalarField>();

        // Extract the KZH accumulator out of `proof_i`:
        let kzh_opening_proof = proof.eval_vars_at_ry;
        // Hossein: Turn kzh_opening_proof into accumulator: kzh_opening_proof_accumulator
        // For now let's use a random accumulator:
        let kzh_opening_proof_accumulator = Accumulator::random_satisfying_accumulator(&srs, &mut thread_rng());

        // Aggregate running_accumulator with kzh_opening_proof_accumulator into acc_{i+1}
        let (acc_instance, acc_witness, Q) = Accumulator::prove(&srs,
                                                                &running_accumulator,
                                                                &kzh_opening_proof_accumulator);

        // Hossein: Get an AccumulatorVerifier representing the accumulation above
        // We need something like `PartialVerifier::initialise` but for `AccumulatorVerifier`, so that we can do:
        // `AccumulatorVerifier::initialise(running_accumulator, kzh_opening_proof_accumulator, acc_instance)`
        // and it will create a valid `AccumulatorVerifier` object.
        let kzh_acc_verifier = todo!();

        // Let's build the circuit
        let cs = ConstraintSystem::<ScalarField>::new_ref();
        let augmented_circuit = AugmentedCircuit::<G1, G2, C2, E> {
            spartan_partial_verifier: partial_verifier.clone(),
            kzh_acc_verifier: kzh_acc_verifier,
        };

        let augmented_circuit_var = AugmentedCircuitVar::new_variable(
            cs.clone(),
            || Ok(augmented_circuit.clone()),
            AllocationMode::Input,
        ).unwrap();

        // Prove it...

    }
}




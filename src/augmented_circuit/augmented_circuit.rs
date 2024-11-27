#![allow(dead_code)]

use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use crate::accumulation_circuit::verifier_circuit::{AccumulatorVerifier, AccumulatorVerifierVar};
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::hash::poseidon::PoseidonHashVar;
use crate::kzh::kzh2::split_between_x_and_y;
use crate::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationAccVerifier, MatrixEvaluationAccVerifierVar};
use crate::nexus_spartan::partial_verifier::partial_verifier::SpartanPartialVerifier;
use crate::nexus_spartan::partial_verifier::partial_verifier_var::SpartanPartialVerifierVar;
use crate::nova::cycle_fold::coprocessor_constraints::RelaxedOvaInstanceVar;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use itertools::izip;
use rand::thread_rng;
use std::borrow::Borrow;

type Output<'a, G2, C2, G1, F> = (
    (RelaxedOvaInstanceVar<G2, C2>, &'a AccumulatorInstanceVar<G1>),  // accumulator final instance, Ova final instance
    (Vec<FpVar<F>>, Vec<FpVar<F>>), // r_x, r_y
    (Vec<FpVar<F>>, Vec<FpVar<F>>, (FpVar<F>, FpVar<F>, FpVar<F>)), // (vector_x, vector_y, evaluations)
);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AugmentedCircuit<G1, G2, C2, E, F>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    F: PrimeField,
{
    pub spartan_partial_verifier: SpartanPartialVerifier<F, E>,
    pub kzh_acc_verifier: AccumulatorVerifier<G1, G2, C2, E>,
    pub matrix_evaluation_verifier: MatrixEvaluationAccVerifier<F>,
}

pub struct AugmentedCircuitVar<G1, G2, C2, F>
where
    F: PrimeField + Absorb,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
{
    pub spartan_partial_verifier: SpartanPartialVerifierVar<F, G1>,
    pub kzh_acc_verifier: AccumulatorVerifierVar<G1, G2, C2>,
    pub matrix_evaluation_verifier: MatrixEvaluationAccVerifierVar<F>,
}

impl<G1, G2, C2, E, F> AllocVar<AugmentedCircuit<G1, G2, C2, E, F>, F> for AugmentedCircuitVar<G1, G2, C2, F>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    F: PrimeField,
{
    fn new_variable<T: Borrow<AugmentedCircuit<G1, G2, C2, E, F>>>(
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

        // Allocate the Spartan partial verifier
        let spartan_partial_verifier = SpartanPartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.spartan_partial_verifier),
            mode,
        )?;

        // Allocate the accumulator verifier
        let kzh_acc_verifier = AccumulatorVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.kzh_acc_verifier),
            mode,
        )?;

        // Allocate the accumulator verifier
        let matrix_evaluation_verifier = MatrixEvaluationAccVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.matrix_evaluation_verifier),
            mode,
        )?;

        Ok(AugmentedCircuitVar {
            spartan_partial_verifier,
            kzh_acc_verifier,
            matrix_evaluation_verifier,
        })
    }
}

impl<G1, G2, C2, F> AugmentedCircuitVar<G1, G2, C2, F>
where
    F: PrimeField + Absorb,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
{
    pub fn verify<E: Pairing>(&self, cs: ConstraintSystemRef<F>, transcript: &mut TranscriptVar<F>, poseidon_num: usize) -> Output<G2, C2, G1, F> {
        let (rx, ry) = self.spartan_partial_verifier.verify(transcript);
        let (final_cycle_fold_instance, final_accumulator_instance) = self.kzh_acc_verifier.accumulate(transcript);

        // also return these later
        let ((vector_x, vector_y), evaluations) = self.matrix_evaluation_verifier.accumulate(transcript);

        // ************* do the consistency checks *************
        let length_x = self.kzh_acc_verifier.current_accumulator_instance_var.x_var.len();
        let length_y = self.kzh_acc_verifier.current_accumulator_instance_var.y_var.len();

        let (expected_x_var, expected_y_var) = split_between_x_and_y(length_x, length_y, &ry[1..], FpVar::zero());
        for (e1, e2) in izip!(&self.kzh_acc_verifier.current_accumulator_instance_var.x_var, expected_x_var) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }

        for (e1, e2) in izip!(&self.kzh_acc_verifier.current_accumulator_instance_var.y_var, expected_y_var) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }

        // enforce equal eval_Z_at_ry and accumulator.z_var
        self.spartan_partial_verifier.eval_vars_at_ry.enforce_equal(
            &self.kzh_acc_verifier
                .current_accumulator_instance_var
                .z_var
        ).expect("error while enforcing equality");

        // enforce the commitment in spartan verifier and the accumulator new instance
        NonNativeAffineVar::enforce_equal(
            &self.spartan_partial_verifier.instance.1,
            &self.kzh_acc_verifier.current_accumulator_instance_var.C_var,
        ).expect("error while enforcing equality");

        // pad it with some random poseidon hash
        let mut hash = PoseidonHashVar::new(cs.clone());
        for _ in 0..poseidon_num {
            // get a random element
            let r = FpVar::new_variable(cs.clone(), || Ok(F::rand(&mut thread_rng())), AllocationMode::Witness).unwrap();
            // update sponge with this random element
            hash.update_sponge(vec![r]);
            // output the hash
            let _ = hash.output();
        }

        ((final_cycle_fold_instance, final_accumulator_instance), (rx, ry), (vector_x, vector_y, evaluations))
    }
}
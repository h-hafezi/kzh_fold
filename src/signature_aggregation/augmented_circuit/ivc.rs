/*use crate::accumulation_circuit::verifier_circuit::{AccumulatorVerifier, AccumulatorVerifierVar};
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationAccVerifier, MatrixEvaluationAccVerifierVar};
use crate::nexus_spartan::partial_verifier::partial_verifier::SpartanPartialVerifier;
use crate::nexus_spartan::partial_verifier::partial_verifier_var::SpartanPartialVerifierVar;
use crate::pcs::multilinear_pcs::split_between_x_and_y;
use crate::signature_aggregation::verifier_circuit::verifier_circuit::SignatureVerifierCircuit;
use crate::signature_aggregation::verifier_circuit::verifier_circuit_var::SignatureVerifierCircuitVar;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_relations::r1cs::{Namespace, SynthesisError};
use itertools::izip;
use std::borrow::Borrow;
use rand::Rng;
use crate::nexus_spartan::crr1cs::CRR1CSShape;
use crate::nexus_spartan::matrix_evaluation_accumulation::prover::fold_matrices_evaluations;
use crate::transcript::transcript::Transcript;

pub struct SignatureAugmentedCircuit<G1, G2, C2, E, F>
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
    pub kzh_acc_verifier_1: AccumulatorVerifier<G1, G2, C2, E>,
    pub matrix_evaluation_verifier: MatrixEvaluationAccVerifier<F>,
    pub kzh_acc_verifier_2: AccumulatorVerifier<G1, G2, C2, E>,
    pub signature_verifier_circuit: SignatureVerifierCircuit<E, F, G1, G2, C2>,
}

pub struct SignatureAugmentedCircuitVar<G1, G2, C2, F>
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
    pub kzh_acc_verifier_1: AccumulatorVerifierVar<G1, G2, C2>,
    pub matrix_evaluation_verifier: MatrixEvaluationAccVerifierVar<F>,
    pub kzh_acc_verifier_2: AccumulatorVerifierVar<G1, G2, C2>,
    pub signature_verifier_circuit: SignatureVerifierCircuitVar<F, G1, G2, C2>,
}

impl<G1, G2, C2, E, F> AllocVar<SignatureAugmentedCircuit<G1, G2, C2, E, F>, F> for SignatureAugmentedCircuitVar<G1, G2, C2, F>
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
    fn new_variable<T: Borrow<SignatureAugmentedCircuit<G1, G2, C2, E, F>>>(
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
        let kzh_acc_verifier_1 = AccumulatorVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.kzh_acc_verifier_1),
            mode,
        )?;

        // Allocate the accumulator verifier
        let kzh_acc_verifier_2 = AccumulatorVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.kzh_acc_verifier_2),
            mode,
        )?;


        // Allocate the accumulator verifier
        let matrix_evaluation_verifier = MatrixEvaluationAccVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.matrix_evaluation_verifier),
            mode,
        )?;

        // Allocate the accumulator verifier
        let signature_verifier_circuit = SignatureVerifierCircuitVar::new_variable(
            cs.clone(),
            || Ok(&data.signature_verifier_circuit),
            mode,
        )?;


        Ok(SignatureAugmentedCircuitVar {
            spartan_partial_verifier,
            matrix_evaluation_verifier,
            kzh_acc_verifier_2,
            kzh_acc_verifier_1,
            signature_verifier_circuit,
        })
    }
}

impl<G1, G2, C2, F> SignatureAugmentedCircuitVar<G1, G2, C2, F>
where
    F: PrimeField + Absorb,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
{
    fn verify<E: Pairing>(&self, transcript: &mut TranscriptVar<F>) {
        let (rx, ry) = self.spartan_partial_verifier.verify(transcript);
        let _ = self.kzh_acc_verifier_1.accumulate(transcript);
        let _ = self.kzh_acc_verifier_2.accumulate(transcript);

        // verify the signature circuit
        let _ = self.signature_verifier_circuit.verify(transcript);

        // also return these later
        let _ = self.matrix_evaluation_verifier.accumulate(transcript);

        // ************* new consistency checks *************
        //self.kzh_acc_verifier_2.running_accumulator_instance_var.enforce_equal(
        //    &self.kzh_acc_verifier_1.final_accumulator_instance_var
        //).expect("error while enforcing equality");

        //self.kzh_acc_verifier_2.current_accumulator_instance_var.enforce_equal(
        //    &self.signature_verifier_circuit.commitment
        //).expect("error while enforcing equality");

        //self.kzh_acc_verifier_1.final_cycle_fold_instance_var.enforce_equal(
        //    &self.kzh_acc_verifier_2.running_cycle_fold_instance_var
        //).expect("error while enforcing equality");

        // ************* do the consistency checks *************
        let length_x = self.kzh_acc_verifier_1.current_accumulator_instance_var.x_var.len();
        let length_y = self.kzh_acc_verifier_1.current_accumulator_instance_var.y_var.len();

        let (expected_x_var, expected_y_var) = split_between_x_and_y(length_x, length_y, &ry[1..], FpVar::zero());
        for (e1, e2) in izip!(&self.kzh_acc_verifier_1.current_accumulator_instance_var.x_var, expected_x_var) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }

        for (e1, e2) in izip!(&self.kzh_acc_verifier_1.current_accumulator_instance_var.y_var, expected_y_var) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }


        // enforce equal eval_Z_at_ry and accumulator.z_var
        self.spartan_partial_verifier.eval_vars_at_ry.enforce_equal(
            &self.kzh_acc_verifier_1
                .current_accumulator_instance_var
                .z_var
        ).expect("error while enforcing equality");

        // enforce the commitment in spartan verifier and the accumulator new instance
        NonNativeAffineVar::enforce_equal(
            &self.spartan_partial_verifier.instance.1,
            &self.kzh_acc_verifier_1.current_accumulator_instance_var.C_var,
        ).expect("error while enforcing equality");
    }
}
 */

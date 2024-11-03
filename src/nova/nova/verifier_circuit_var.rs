use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::gadgets::r1cs::r1cs_var::{R1CSInstanceVar, RelaxedR1CSInstanceVar};
use crate::gadgets::r1cs::{OvaInstance, R1CSInstance, RelaxedOvaInstance, RelaxedR1CSInstance};
use crate::nova::cycle_fold::coprocessor_constraints::{OvaInstanceVar, RelaxedOvaInstanceVar};
use crate::nova::nova::verifier_circuit::NovaAugmentedCircuit;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

pub struct NovaAugmentedCircuitVar<F, G1, G2, C1, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>, Commitment=Projective<G1>, SetupAux=()>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>, Commitment=Projective<G2>, SetupAux=()>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    F: PrimeField,
{
    pub running_instance: RelaxedR1CSInstanceVar<G1, C1>,
    pub final_instance: RelaxedR1CSInstanceVar<G1, C1>,
    pub current_instance: R1CSInstanceVar<G1, C1>,

    /// nova accumulation prof (cross term error)
    pub nova_cross_term_error: NonNativeAffineVar<G1>,

    /// this is hash of two instance and nova_cross_term_error
    pub beta: FpVar<F>,
    pub beta_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// beta_1 is the randomness used to fold cycle fold running instance with auxiliary_input_W
    pub beta_1: FpVar<F>,
    pub beta_1_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// beta_2 is the randomness used to fold cycle fold running instance with auxiliary_input_E
    pub beta_2: FpVar<F>,
    pub beta_2_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// running cycle fold instance
    pub running_cycle_fold_instance: RelaxedOvaInstanceVar<G2, C2>,
    pub final_cycle_fold_instance: RelaxedOvaInstanceVar<G2, C2>,

    /// auxiliary input which helps to have W = W_1 + beta * W_2 without scalar multiplication
    pub auxiliary_input_W_var: OvaInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E = E_1 + beta * com_T without scalar multiplication
    pub auxiliary_input_E_var: OvaInstanceVar<G2, C2>,

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub cross_term_error_commitment_w: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub cross_term_error_commitment_e: ProjectiveVar<G2, FpVar<G2::BaseField>>,
}

impl<F, G1, G2, C1, C2> AllocVar<NovaAugmentedCircuit<F, G1, G2, C1, C2>, G1::ScalarField> for NovaAugmentedCircuitVar<F, G1, G2, C1, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>, Commitment=Projective<G1>, SetupAux=()>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>, Commitment=Projective<G2>, SetupAux=()>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    F: PrimeField,
{
    fn new_variable<T: Borrow<NovaAugmentedCircuit<F, G1, G2, C1, C2>>>(
        cs: impl Into<Namespace<G1::ScalarField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        // Allocation for running_instance
        let running_instance = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "running_instance"),
            || circuit.map(|e| e.running_instance.clone()),
            mode,
        )?;

        // Allocation for final_instance
        let final_instance = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "final_instance"),
            || circuit.map(|e| e.final_instance.clone()),
            mode,
        )?;

        // Allocation for current_instance
        let current_instance = R1CSInstanceVar::new_variable(
            ns!(cs, "current_instance"),
            || circuit.map(|e| e.current_instance.clone()),
            mode,
        )?;

        // Allocation for nova_cross_term_error
        let nova_cross_term_error = NonNativeAffineVar::new_variable(
            ns!(cs, "nova_cross_term_error"),
            || circuit.map(|e| e.nova_cross_term_error.clone()),
            mode,
        )?;

        // Allocation for beta
        let beta = FpVar::new_variable(
            ns!(cs, "beta"),
            || circuit.map(|e| e.beta.clone()),
            mode,
        )?;

        // Allocation for beta_non_native
        let beta_non_native = NonNativeFieldVar::new_variable(
            ns!(cs, "beta_non_native"),
            || circuit.map(|e| e.beta_non_native.clone()),
            mode,
        )?;

        // Allocation for beta_1
        let beta_1 = FpVar::new_variable(
            ns!(cs, "beta_1"),
            || circuit.map(|e| e.beta_1.clone()),
            mode,
        )?;

        // Allocation for beta_1_non_native
        let beta_1_non_native = NonNativeFieldVar::new_variable(
            ns!(cs, "beta_1_non_native"),
            || circuit.map(|e| e.beta_1_non_native.clone()),
            mode,
        )?;

        // Allocation for beta_2
        let beta_2 = FpVar::new_variable(
            ns!(cs, "beta_2"),
            || circuit.map(|e| e.beta_2.clone()),
            mode,
        )?;

        // Allocation for beta_2_non_native
        let beta_2_non_native = NonNativeFieldVar::new_variable(
            ns!(cs, "beta_2_non_native"),
            || circuit.map(|e| e.beta_2_non_native.clone()),
            mode,
        )?;

        // Allocation for running_cycle_fold_instance
        let running_cycle_fold_instance = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "running_cycle_fold_instance"),
            || circuit.map(|e| e.running_cycle_fold_instance.clone()),
            mode,
        )?;

        // Allocation for final_cycle_fold_instance
        let final_cycle_fold_instance = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "final_cycle_fold_instance"),
            || circuit.map(|e| e.final_cycle_fold_instance.clone()),
            mode,
        )?;

        // Allocation for auxiliary_input_W_var
        let auxiliary_input_W_var = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_W_var"),
            || circuit.map(|e| e.auxiliary_input_W_var.clone()),
            mode,
        )?;

        // Allocation for auxiliary_input_E_var
        let auxiliary_input_E_var = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_var"),
            || circuit.map(|e| e.auxiliary_input_E_var.clone()),
            mode,
        )?;

        // Allocation for cross_term_error_commitment_w
        let cross_term_error_commitment_w = ProjectiveVar::new_variable(
            ns!(cs, "cross_term_error_commitment_w"),
            || circuit.map(|e| e.cross_term_error_commitment_w.clone()),
            mode,
        )?;

        // Allocation for cross_term_error_commitment_e
        let cross_term_error_commitment_e = ProjectiveVar::new_variable(
            ns!(cs, "cross_term_error_commitment_e"),
            || circuit.map(|e| e.cross_term_error_commitment_e.clone()),
            mode,
        )?;

        Ok(NovaAugmentedCircuitVar {
            running_instance,
            final_instance,
            current_instance,
            nova_cross_term_error,
            beta,
            beta_non_native,
            beta_1,
            beta_1_non_native,
            beta_2,
            beta_2_non_native,
            running_cycle_fold_instance,
            final_cycle_fold_instance,
            auxiliary_input_W_var,
            auxiliary_input_E_var,
            cross_term_error_commitment_w,
            cross_term_error_commitment_e,
        })
    }
}

impl<F, G1, G2, C1, C2> ConstraintSynthesizer<F> for NovaAugmentedCircuitVar<F, G1, G2, C1, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>, Commitment=Projective<G1>, SetupAux=()>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>, Commitment=Projective<G2>, SetupAux=()>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    F: PrimeField,
{
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        // check that the linear combination of the running instance, new instance and finally instance for X

        // check that secondary circuit W is consistent with the current/new/final instances

        // check that secondary circuit E consistency too

        // fold the cycle fold instances and check consistency with the provided one

        // compute beta, beta_1, beta_2 and check their correctness as well consistency with their non-native representation

        Ok(())
    }
}
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::gadgets::r1cs::r1cs_var::{R1CSInstanceVar, RelaxedR1CSInstanceVar};
use crate::nova::cycle_fold::coprocessor_constraints::{OvaInstanceVar, RelaxedOvaInstanceVar};
use crate::nova::nova::get_affine_var_coords;
use crate::nova::nova::verifier_circuit::NovaAugmentedCircuit;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::ToBitsGadget;
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
    pub ova_running_instance: RelaxedOvaInstanceVar<G2, C2>,

    /// auxiliary input which helps to have W = W_1 + beta * W_2 without scalar multiplication
    pub ova_auxiliary_input_W: OvaInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E = E_1 + beta * com_T without scalar multiplication
    pub ova_auxiliary_input_E: OvaInstanceVar<G2, C2>,

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub ova_cross_term_error_commitment_w: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub ova_cross_term_error_commitment_e: ProjectiveVar<G2, FpVar<G2::BaseField>>,
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
            || circuit.map(|e| e.running_ova_instance.clone()),
            mode,
        )?;

        // Allocation for auxiliary_input_W_var
        let ova_auxiliary_input_W = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_W_var"),
            || circuit.map(|e| e.auxiliary_input_W_var.clone()),
            mode,
        )?;

        // Allocation for auxiliary_input_E_var
        let ova_auxiliary_input_E = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_var"),
            || circuit.map(|e| e.auxiliary_input_E_var.clone()),
            mode,
        )?;

        // Allocation for cross_term_error_commitment_w
        let ova_cross_term_error_commitment_w = ProjectiveVar::new_variable(
            ns!(cs, "cross_term_error_commitment_w"),
            || circuit.map(|e| e.cross_term_error_commitment_w.clone()),
            mode,
        )?;

        // Allocation for cross_term_error_commitment_e
        let ova_cross_term_error_commitment_e = ProjectiveVar::new_variable(
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
            ova_running_instance: running_cycle_fold_instance,
            ova_auxiliary_input_W,
            ova_auxiliary_input_E,
            ova_cross_term_error_commitment_w,
            ova_cross_term_error_commitment_e,
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
        for i in 0..self.current_instance.X.len() {
            let res = &self.running_instance.X[i] + &self.beta * &self.current_instance.X[i];
            self.final_instance.X[i].enforce_equal(&res).expect("error while enforcing equality");
        }

        assert_eq!(self.current_instance.X.len(), self.final_instance.X.len(), "instances must have equal size");
        assert_eq!(self.current_instance.X.len(), self.running_instance.X.len(), "instances must have equal size");


        // check that secondary circuit W is consistent with the current/new/final instances
        let (flag,
            r,
            g1,
            g2,
            g_out
        ) = self.ova_auxiliary_input_W.parse_secondary_io::<G1>().unwrap();

        r.enforce_equal(&self.beta_non_native).expect("error while enforcing equality");
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        g1.enforce_equal(&self.current_instance.commitment_W).expect("error while enforcing equality");
        g2.enforce_equal(&self.running_instance.commitment_W).expect("error while enforcing equality");
        g_out.enforce_equal(&self.final_instance.commitment_W).expect("error while enforcing equality");


        // check that secondary circuit E consistency too
        let (flag,
            r,
            g1,
            g2,
            g_out
        ) = self.ova_auxiliary_input_E.parse_secondary_io::<G1>().unwrap();

        r.enforce_equal(&self.beta_non_native).expect("error while enforcing equality");
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        g1.enforce_equal(&self.nova_cross_term_error).expect("error while enforcing equality");
        g2.enforce_equal(&self.running_instance.commitment_E).expect("error while enforcing equality");
        g_out.enforce_equal(&self.final_instance.commitment_E).expect("error while enforcing equality");


        // fold the cycle fold instances and check consistency with the provided one
        let beta_1_bits: Vec<Boolean<F>> = self.beta_1_non_native.to_bits_le().unwrap();
        let beta_2_bits: Vec<Boolean<F>> = self.beta_2_non_native.to_bits_le().unwrap();

        let _ova_final_instance = self.ova_running_instance.fold(
            &[((&self.ova_auxiliary_input_W, None), &self.ova_cross_term_error_commitment_w, &self.beta_1_non_native, &beta_1_bits),
                ((&self.ova_auxiliary_input_E, None), &self.ova_cross_term_error_commitment_e, &self.beta_2_non_native, &beta_2_bits)
            ]
        ).unwrap();

        // compute beta, beta_1, beta_2 and check their correctness as well consistency with their non-native representation

        // compute beta
        let mut transcript = TranscriptVar::new(cs.clone(), b"new transcript");
        transcript.append_scalars(b"label", self.running_instance.to_sponge_field_elements().as_slice());
        transcript.append_scalars(b"label", self.current_instance.to_sponge_field_elements().as_slice());
        transcript.append_scalars(b"label", self.nova_cross_term_error.to_sponge_field_elements().unwrap().as_slice());
        let beta = transcript.challenge_scalar(b"challenge");

        // consistency checks of beta
        beta.enforce_equal(&self.beta).expect("error while enforcing equality");
        beta.enforce_equal(&{
            let beta_bits = self.beta_non_native.to_bits_le().unwrap();
            Boolean::le_bits_to_fp_var(beta_bits.as_slice()).unwrap()
        }).expect("error while enforcing equality");


        // derive beta_1 and its consistency checks
        transcript.append_scalars(
            b"add scalars",
            get_affine_var_coords(&self.ova_cross_term_error_commitment_w.to_affine()?).as_slice(),
        );
        let beta_1 = transcript.challenge_scalar(b"challenge");
        beta_1.enforce_equal(&self.beta_1).expect("error while enforcing equality");
        beta_1.enforce_equal(&Boolean::le_bits_to_fp_var(beta_1_bits.as_slice())?).expect("error while enforcing equality");


        // derive beta_2 and its consistency checks
        transcript.append_scalars(
            b"add scalars",
            get_affine_var_coords(&self.ova_cross_term_error_commitment_e.to_affine()?).as_slice(),
        );
        let beta_2 = transcript.challenge_scalar(b"challenge");
        beta_2.enforce_equal(&self.beta_2).expect("error while enforcing equality");
        beta_2.enforce_equal(&Boolean::le_bits_to_fp_var(beta_2_bits.as_slice()).unwrap()).expect("error while enforcing equality");

        // this randomness should be used to get the randomness beta_1 and beta_2, we currently simply add them but don't use them
        transcript.append_scalars_non_native(
            b"label",
            self.ova_running_instance.X.as_slice(),
        );
        transcript.append_scalars_non_native(
            b"label",
            self.ova_auxiliary_input_E.X.as_slice(),
        );
        transcript.append_scalars_non_native(
            b"label",
            self.ova_auxiliary_input_W.X.as_slice(),
        );
        transcript.append_scalars(
            b"label",
            &[
                self.ova_auxiliary_input_W.commitment.x.clone(),
                self.ova_auxiliary_input_W.commitment.y.clone(),
                self.ova_auxiliary_input_W.commitment.z.clone(),
                self.ova_auxiliary_input_E.commitment.x.clone(),
                self.ova_auxiliary_input_E.commitment.y.clone(),
                self.ova_auxiliary_input_E.commitment.z.clone(),
                self.ova_running_instance.commitment.x.clone(),
                self.ova_running_instance.commitment.y.clone(),
                self.ova_running_instance.commitment.z.clone(),
            ],
        );

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::constant_for_curves::{ScalarField, C1, C2, G1, G2};
    use crate::nova::nova::prover::NovaProver;
    use crate::nova::nova::verifier_circuit::NovaAugmentedCircuit;
    use crate::nova::nova::verifier_circuit_var::NovaAugmentedCircuitVar;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem};

    type F = ScalarField;

    #[test]
    fn test_synthesize() {
        let prover: NovaProver<F, G1, G2, C1, C2> = NovaProver::rand((10, 3, 17));
        let cs = ConstraintSystem::<F>::new_ref();

        let augmented_circuit: NovaAugmentedCircuit<F, G1, G2, C1, C2> = NovaAugmentedCircuit::initialise(prover);

        augmented_circuit.verify();

        let augmented_circuit_var: NovaAugmentedCircuitVar<F, G1, G2, C1, C2> = NovaAugmentedCircuitVar::new_variable(
            cs.clone(),
            || Ok(augmented_circuit.clone()),
            AllocationMode::Input,
        ).unwrap();

        augmented_circuit_var.generate_constraints(cs.clone()).expect("error");
        println!("{}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap());
    }
}

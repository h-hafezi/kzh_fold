use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::util::convert_field_one_to_field_two;
use crate::gadgets::r1cs::{OvaInstance, R1CSInstance, RelaxedOvaInstance, RelaxedR1CSInstance};
use crate::nova::nova::prover::NovaProver;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::CurveConfig;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;

pub struct NovaAugmentedCircuit<F, G1, G2, C1, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>>,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    F: PrimeField,
{
    pub running_instance: RelaxedR1CSInstance<G1, C1>,
    pub final_instance: RelaxedR1CSInstance<G1, C1>,
    pub current_instance: R1CSInstance<G1, C1>,

    /// nova accumulation prof (cross term error)
    pub nova_cross_term_error: C1::Commitment,

    // this is hash of two instance and nova_cross_term_error
    pub beta: F,
    pub beta_non_native: G1::BaseField,

    /// running cycle fold instance
    pub running_cycle_fold_instance: RelaxedOvaInstance<G2, C2>,
    pub final_cycle_fold_instance: RelaxedOvaInstance<G2, C2>,

    /// auxiliary input which helps to have W = W_1 + beta * W_2 without scalar multiplication
    pub auxiliary_input_W_var: OvaInstance<G2, C2>,
    /// auxiliary input which helps to have E = E_1 + beta * com_T without scalar multiplication
    pub auxiliary_input_E_var: OvaInstance<G2, C2>,

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub cross_term_error_commitment_w: C2::Commitment,
    pub cross_term_error_commitment_e: C2::Commitment,
}

impl<F, G1, G2, C1, C2> NovaAugmentedCircuit<F, G1, G2, C1, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>, Commitment=Projective<G1>, SetupAux = ()>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>, SetupAux = ()>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    F: PrimeField,
{
    pub fn verify(&self) {
        let expected_final_instance = self.running_instance.fold(
            &self.current_instance,
            &self.nova_cross_term_error,
            &self.beta,
        ).unwrap();

        assert_eq!(expected_final_instance, self.final_instance);

        let secondary_circuit_W = self.auxiliary_input_W_var.parse_secondary_io::<G1>().unwrap();

        assert_eq!(secondary_circuit_W.r, self.beta_non_native);
        assert_eq!(secondary_circuit_W.g1, self.current_instance.commitment_W);
        assert_eq!(secondary_circuit_W.g2, self.running_instance.commitment_W);
        assert_eq!(secondary_circuit_W.g_out, self.final_instance.commitment_W);
        assert_eq!(secondary_circuit_W.flag, true);

        let secondary_circuit_E = self.auxiliary_input_E_var.parse_secondary_io::<G1>().unwrap();

        assert_eq!(secondary_circuit_E.r, self.beta_non_native);
        assert_eq!(secondary_circuit_E.g1, self.nova_cross_term_error);
        assert_eq!(secondary_circuit_E.g2, self.running_instance.commitment_E);
        assert_eq!(secondary_circuit_E.g_out, self.final_instance.commitment_E);
        assert_eq!(secondary_circuit_E.flag, true);

        let expected_final_cycle_fold_instance = {
            let temp = self.running_cycle_fold_instance.fold(&self.auxiliary_input_W_var, &self.cross_term_error_commitment_w, &self.beta_non_native).expect("TODO: panic message");
            temp.fold(&self.auxiliary_input_E_var, &self.cross_term_error_commitment_e, &self.beta_non_native).expect("TODO: panic message")
        };

        assert_eq!(expected_final_cycle_fold_instance, self.final_cycle_fold_instance);

        // add beta tests too
    }

    pub fn initialise(prover: NovaProver<F, G1, G2, C1, C2>) -> NovaAugmentedCircuit<F, G1, G2, C1, C2> where <G2 as CurveConfig>::ScalarField: Absorb {
        let beta = prover.compute_beta();
        let beta_non_native = convert_field_one_to_field_two::<G1::ScalarField, G1::BaseField>(beta);
        let (final_instance, _, nova_cross_term_error) = prover.compute_final_accumulator(&beta);
        let (final_cycle_fold_instance, _, cross_term_error_commitment_w, cross_term_error_commitment_e) = prover.compute_final_cycle_fold_instance(&beta);

        NovaAugmentedCircuit {
            running_instance: prover.running_accumulator.0.clone(),
            final_instance,
            current_instance: prover.current_accumulator.0.clone(),
            nova_cross_term_error,
            beta,
            beta_non_native,
            running_cycle_fold_instance: prover.cycle_fold_running_instance.clone(),
            final_cycle_fold_instance,
            auxiliary_input_W_var: prover.compute_auxiliary_input_W(&beta).0,
            auxiliary_input_E_var: prover.compute_auxiliary_input_E(&beta).0,
            cross_term_error_commitment_w,
            cross_term_error_commitment_e,
        }
    }
}

#[cfg(test)]
mod test {
    use crate::constant_for_curves::{ScalarField, C1, C2, G1, G2};
    use crate::nova::nova::prover::NovaProver;
    use crate::nova::nova::verifier_circuit::NovaAugmentedCircuit;

    type F = ScalarField;

    #[test]
    fn test() {
        let prover: NovaProver<F, G1, G2, C1, C2> = NovaProver::rand((10, 3, 7));

        let augmented_circuit: NovaAugmentedCircuit<F, G1, G2, C1, C2> = NovaAugmentedCircuit::initialise(prover);

        augmented_circuit.verify();
    }
}



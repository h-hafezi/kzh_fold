use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::{OvaInstance, R1CSInstance, RelaxedOvaInstance, RelaxedR1CSInstance};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
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
    pub com_T: C1::Commitment,

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
    pub com_W: C2::Commitment,
    pub com_E: C2::Commitment,
}

impl<F, G1, G2, C1, C2> NovaAugmentedCircuit<F, G1, G2, C1, C2>
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
    pub fn verify(&self) {
        let expected_final_instance = self.running_instance.fold(
            &self.current_instance,
            &self.com_T,
            &self.beta,
        ).unwrap();

        assert_eq!(expected_final_instance, self.final_instance);

        let secondary_circuit_W = self.auxiliary_input_W_var.parse_secondary_io::<G1>().unwrap();

        assert_eq!(secondary_circuit_W.r, self.beta_non_native);
        assert_eq!(secondary_circuit_W.g1, self.current_instance.commitment_W.into());
        assert_eq!(secondary_circuit_W.g2, self.running_instance.commitment_W.into());
        assert_eq!(secondary_circuit_W.g_out, self.final_instance.commitment_W.into());
        assert_eq!(secondary_circuit_W.flag, true);

        let secondary_circuit_E = self.auxiliary_input_E_var.parse_secondary_io::<G1>().unwrap();

        assert_eq!(secondary_circuit_E.r, self.beta_non_native);
        assert_eq!(secondary_circuit_E.g1, self.com_T.into());
        assert_eq!(secondary_circuit_E.g2, self.running_instance.commitment_E.into());
        assert_eq!(secondary_circuit_E.g_out, self.final_instance.commitment_E.into());
        assert_eq!(secondary_circuit_E.flag, true);

        let expected_final_cycle_fold_instance = {
            let temp = self.running_cycle_fold_instance.fold(&self.auxiliary_input_W_var, &self.com_W, &self.beta_non_native).expect("TODO: panic message");
            temp.fold(&self.auxiliary_input_E_var, &self.com_E, &self.beta_non_native).expect("TODO: panic message")
        };

        assert_eq!(expected_final_cycle_fold_instance, self.final_cycle_fold_instance);
    }
}



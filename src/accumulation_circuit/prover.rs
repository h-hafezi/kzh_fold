use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::PrimeField;

use crate::accumulation::accumulator::Accumulator;
use crate::gadgets::r1cs::{R1CSInstance, R1CSWitness, RelaxedR1CSInstance, RelaxedR1CSWitness};
use crate::nova::commitment::CommitmentScheme;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccumulatorVerifierCircuitProver<G1, G2, C2, E>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing,
{
    /// the randomness used for taking linear combination
    pub beta: G1::ScalarField,

    /// the instance to be folded
    pub instance: Accumulator<E>,
    /// the running accumulator
    pub acc: Accumulator<E>,

    /// running cycle fold instance
    pub cycle_fold_running_instance: RelaxedR1CSInstance<G2, C2>,
    pub cycle_fold_running_witness: RelaxedR1CSWitness<G2>,

    // these are constant values
    pub n: u32,
    pub m: u32,
}

pub trait AccumulatorVerifierCircuitProverTrait<E, G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing,
{
    fn compute_auxiliary_input_C(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>);

    fn compute_auxiliary_input_T(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>);

    fn compute_auxiliary_input_E_1(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>);

    fn compute_auxiliary_input_E_2(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>);

    fn compute_proof_Q(&self) -> Projective<G1>;

    fn compute_cycle_fold_proofs(&self) -> (Projective<G2>, Projective<G2>, Projective<G2>, Projective<G2>);
}

impl<G1, G2, C2, E> AccumulatorVerifierCircuitProverTrait<E, G1, G2, C2> for AccumulatorVerifierCircuitProver<G1, G2, C2, E> {
    fn compute_auxiliary_input_C(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>) {
        unimplemented!()
    }

    fn compute_auxiliary_input_T(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>) {
        unimplemented!()
    }

    fn compute_auxiliary_input_E_1(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>) {
        unimplemented!()
    }

    fn compute_auxiliary_input_E_2(&self) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>) {
        unimplemented!()
    }

    fn compute_proof_Q(&self) -> Projective<G1> {
        unimplemented!()
    }

    fn compute_cycle_fold_proofs(&self) -> (Projective<G2>, Projective<G2>, Projective<G2>, Projective<G2>) {
        unimplemented!()
    }
}
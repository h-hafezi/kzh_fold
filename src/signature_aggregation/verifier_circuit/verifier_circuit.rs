use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::{OvaInstance, RelaxedOvaInstance};
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::CurveConfig;
use ark_ff::PrimeField;

pub struct SignatureVerifierCircuit<F, G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<
        BaseField=G2::ScalarField,
        ScalarField=G2::BaseField
    >,
    F: PrimeField + Absorb,
{
    /// public keys pk_t = pk_1 + pk_2
    pub running_pk: Projective<G1>,
    pub current_pk: Projective<G1>,
    pub final_pk: Projective<G1>,

    /// randomness for homomorphically combining bitfield commitments
    pub c_0: F,
    pub c_0_non_native: G1::BaseField,
    pub c_1: F,
    pub c_1_non_native: G1::BaseField,

    /// commitments to bitfields
    pub com_bitfield_C: Projective<G1>,
    pub com_bitfield_B_1: Projective<G1>,
    pub com_bitfield_B_2: Projective<G1>,

    /// com_homomorphic_bitfield = com_bitfield_B_1 + c_0 * com_bitfield_B_2 + c1 * com_bitfield_C
    pub com_homomorphic_bitfield: Projective<G1>,

    /// beta used to take linear combination for ova instances
    pub beta: F,
    pub beta_non_native: G1::BaseField,


    /// adding pk
    pub ova_cross_term_error_pk: Projective<G2>,
    pub ova_auxiliary_input_pk: OvaInstance<G2, C2>,

    /// adding temp = B1 + c_0 * B2
    pub ova_cross_term_error_bitfield_1: Projective<G2>,
    pub ova_auxiliary_input_bitfield_1: OvaInstance<G2, C2>,

    /// adding res = temp + c1 * C
    pub ova_cross_term_error_bitfield_2: Projective<G2>,
    pub ova_auxiliary_input_bitfield_2: OvaInstance<G2, C2>,

    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub ova_running_instance: RelaxedOvaInstance<G2, C2>,
    pub ova_final_instance: RelaxedOvaInstance<G2, C2>,

    /// the sumcheck proof
    pub sumcheck_proof: SumcheckCircuit<F>,

    /// Evaluations of the inner polynomials at rho:
    pub b_1_at_rho: F,
    pub b_2_at_rho: F,
    pub c_at_rho: F,

    /// size of the bitfield
    pub bitfield_num_variables: usize,
}

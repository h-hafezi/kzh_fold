use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::{OvaInstance, RelaxedOvaInstance};
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::CurveConfig;
use ark_ff::PrimeField;

pub struct SignatureVerifierCircuit<E, F, G1, G2, C2>
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
    E: Pairing<
        G1Affine=Affine<G1>,
        ScalarField=<G1 as CurveConfig>::ScalarField,
    >,
    F: PrimeField + Absorb,
{
    /// public keys pk_t = pk_1 + pk_2
    pub pk_1: E::G1Affine,
    pub pk_2: E::G1Affine,
    pub pk_t: E::G1Affine,

    /// bitfield commitment, solely used for Fiat-Shamir
    pub com_bitfield: Projective<G1>,

    /// XXX: All the cyclefold things should be ideally combined into a single thing
    /// beta used to take linear combination for cycle fold
    pub beta: G1::BaseField,
    /// the cross term error
    pub com_pk: Projective<G2>,
    pub cycle_fold_fresh_instance: OvaInstance<G2, C2>,
    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub cycle_fold_running_instance: RelaxedOvaInstance<G2, C2>,
    pub cycle_fold_final_instance: RelaxedOvaInstance<G2, C2>,

    /// the sumcheck proof
    pub sumcheck_proof: SumcheckCircuit<F>,

    /// Evaluations of the inner polynomials at rho:
    pub b_1_at_rho: F,
    pub b_2_at_rho: F,
    pub c_at_rho: F,

    /// size of the bitfield
    pub bitfield_num_variables: usize,
}

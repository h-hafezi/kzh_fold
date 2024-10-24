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

    /// the cross term error
    pub com_pk: Projective<G2>,
    pub auxiliary_input_pk: OvaInstance<G2, C2>,
    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub running_auxiliary_input_pk: RelaxedOvaInstance<G2, C2>,
    pub final_auxiliary_input_pk: RelaxedOvaInstance<G2, C2>,

    /// the bitfield polynomial
    pub bitfield_poly: MultilinearPolynomial<F>,
    /// commitment to bitfield (only for Fiat-Shamir)
    pub bitfield_poly_commitment: E::G1Affine,

    /// the sumcheck proof
    pub sumcheck_proof: SumcheckCircuit<F>,

    /// Evaluations of the inner polynomials at rho:
    pub b_1_at_rho: F,
    pub b_2_at_rho: F,
    pub c_at_rho: F,
}


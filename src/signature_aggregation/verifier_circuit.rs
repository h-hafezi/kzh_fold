use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::pcs::multilinear_pcs::{Commitment};

pub struct SignatureVerifierCircuit<E: Pairing<ScalarField=F>, F: PrimeField> {
    // public keys pk_t = pk_1 + pk_2
    pk_1: E::G1Affine,
    pk_2: E::G1Affine,
    pk_t: E::G1Affine,

    // signatures sig_t = sig_1 + sig_2
    sig_1: E::G2Affine,
    sig_2: E::G2Affine,
    sig_t: E::G2Affine,

    // the bitfield polynomial
    bitfield_poly: MultilinearPolynomial<F>,

    // Commitment to c(x)
    bitfield_commitment: Commitment<E>,
    sumcheck_proof: SumcheckInstanceProof<F>,

    // Evaluations of the inner polynomials at rho:
    b_1_at_rho: F, // b_1(rho)
    b_2_at_rho: F, // b_2(rho)
    c_at_rho: F, // c(rho)
}
use ark_crypto_primitives::sponge::Absorb;
use crate::nexus_spartan::unipoly::unipoly::CompressedUniPoly;
use ark_ff::PrimeField;
use crate::nexus_spartan::errors::ProofVerifyError;
use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::transcript::transcript::Transcript;

pub struct SumcheckCircuit<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPoly<F>>,
    claim: F,
    num_rounds: usize,
    degree_bound: usize,
    transcript: Transcript<F>,
}


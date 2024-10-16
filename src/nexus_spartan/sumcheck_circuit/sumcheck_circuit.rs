use crate::nexus_spartan::unipoly::CompressedUniPoly;
use ark_ff::PrimeField;
use crate::nexus_spartan::errors::ProofVerifyError;
use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::transcript::transcript::Transcript;

pub struct SumcheckCircuit<F: PrimeField> {
    pub compressed_polys: Vec<CompressedUniPoly<F>>,
    claim: F,
    num_rounds: usize,
    degree_bound: usize,
    transcript: Transcript<F>,
}


use crate::nexus_spartan::unipoly::CompressedUniPoly;
use ark_ff::PrimeField;
use merlin::Transcript;
use crate::nexus_spartan::errors::ProofVerifyError;
use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;

pub struct SumcheckCircuit<F: PrimeField> {
    pub compressed_polys: Vec<CompressedUniPoly<F>>,
    claim: F,
    num_rounds: usize,
    degree_bound: usize,
    transcript: Transcript,
}


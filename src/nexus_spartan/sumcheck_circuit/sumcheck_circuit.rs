use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::nexus_spartan::unipoly::unipoly::CompressedUniPoly;
use crate::nexus_spartan::unipoly::unipoly_var::{CompressedUniPolyVar, UniPolyVar};
use crate::transcript::transcript::Transcript;
use crate::transcript::transcript_var::{AppendToTranscriptVar, TranscriptVar};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use std::borrow::Borrow;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit_var::SumcheckCircuitVar;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SumcheckCircuit<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPoly<F>>,
    pub claim: F,
    pub num_rounds: usize,
    pub degree_bound: usize,
}

impl<F: PrimeField + Absorb> SumcheckCircuit<F> {
    pub fn verify<E: Pairing<ScalarField=F>>(&self, transcript: &mut Transcript<F>) -> (F, Vec<F>) {
        let proof = SumcheckInstanceProof::new(self.compressed_polys.clone());
        proof.verify::<E>(self.claim, self.num_rounds, self.degree_bound, transcript).unwrap()
    }
}


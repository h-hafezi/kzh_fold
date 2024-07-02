use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::nexus_spartan::unipoly::unipoly::CompressedUniPoly;
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;

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


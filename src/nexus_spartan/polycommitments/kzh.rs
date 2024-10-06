use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_poly_commit::Error;
use merlin::Transcript;
use rand::RngCore;

use crate::nexus_spartan::polycommitments::{PCSKeys, PolyCommitmentScheme};
use crate::nexus_spartan::polycommitments::error::PCSError;
use crate::nexus_spartan::transcript::{AppendToTranscript, ProofTranscript};
use crate::pcs::multilinear_pcs::{Commitment, OpeningProof, PolyCommit, SRS};
use crate::polynomial::multilinear_poly::MultilinearPolynomial;

impl<E: Pairing> AppendToTranscript<E> for Commitment<E> {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        <Transcript as ProofTranscript<E>>::append_point(transcript, label, &self.C);
    }
}


impl<F: PrimeField, E: Pairing<ScalarField=F>> PolyCommitmentScheme<E> for MultilinearPolynomial<F> {
    type SRS = SRS<E>;
    type PolyCommitmentKey = PolyCommit<E>;
    type EvalVerifierKey = PolyCommit<E>;
    type Commitment = Commitment<E>;
    type PolyCommitmentProof = OpeningProof<E>;

    fn commit(poly: &MultilinearPolynomial<E::ScalarField>, ck: &Self::PolyCommitmentKey) -> Self::Commitment {
        ck.commit(&poly)
    }

    fn prove(C: Option<&Self::Commitment>, poly: &MultilinearPolynomial<E::ScalarField>, r: &[E::ScalarField], ck: &Self::PolyCommitmentKey) -> Self::PolyCommitmentProof {
        let (x, _) = ck.srs.split_between_x_and_y(r);
        ck.open(&poly, C.unwrap().clone(), x.as_slice())
    }

    fn verify(commitment: &Self::Commitment, proof: &Self::PolyCommitmentProof, ck: &Self::EvalVerifierKey, r: &[E::ScalarField], eval: &E::ScalarField) -> Result<(), PCSError> {
        let (x, y) = ck.srs.split_between_x_and_y(r);

        // verify the proof
        assert!(ck.verify(commitment, proof, x.as_slice(), y.as_slice(), eval));

        Ok(())
    }

    fn setup(mut max_poly_vars: usize, rng: &mut impl RngCore) -> Result<Self::SRS, Error> {
        let x = max_poly_vars / 2;
        let y = max_poly_vars - x;
        let degree_x = 2usize.pow(x as u32);
        let degree_y = 2usize.pow(y as u32);
        Ok(PolyCommit::<E>::setup(degree_x, degree_y, rng))
    }

    fn trim(srs: &Self::SRS) -> PCSKeys<E, Self> {
        PCSKeys {
            ck: PolyCommit { srs: srs.clone() },
            vk: PolyCommit { srs: srs.clone() },
        }
    }
}

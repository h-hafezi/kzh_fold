use std::io::{Read, Write};
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_poly_commit::Error;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate};
use merlin::Transcript;
use rand::{RngCore, thread_rng};
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
    type PolyCommitmentKey = SRS<E>;
    type EvalVerifierKey = SRS<E>;
    type Commitment = Commitment<E>;
    type PolyCommitmentProof = OpeningProof<E>;

    fn commit(poly: &MultilinearPolynomial<E::ScalarField>, ck: &Self::PolyCommitmentKey) -> Self::Commitment {
        // assert if the poly length is 2^num of variables otherwise append zeros
        //let poly = poly.extend_number_of_variables(ck.get_y_length()+ck.get_x_length());

        let poly_commit: PolyCommit<E> = PolyCommit { srs: ck.clone() };
        poly_commit.commit(&poly)
    }

    fn prove(C: Option<&Self::Commitment>, poly: &MultilinearPolynomial<E::ScalarField>, r: &[E::ScalarField], eval: &E::ScalarField, ck: &Self::PolyCommitmentKey, transcript: &mut Transcript) -> Self::PolyCommitmentProof {
        // assert if the poly length is 2^num of variables otherwise append zeros
        //let poly = poly.extend_number_of_variables(ck.get_y_length()+ck.get_x_length());

        let poly_commit: PolyCommit<E> = PolyCommit { srs: ck.clone() };
        let (x, _) = ck.split_between_x_and_y(r);

        // assert if the poly length is 2^num of variables otherwise append zeros
        poly_commit.open(&poly, C.unwrap().clone(), x.as_slice())
    }

    fn verify(commitment: &Self::Commitment, proof: &Self::PolyCommitmentProof, ck: &Self::EvalVerifierKey, transcript: &mut Transcript, r: &[E::ScalarField], eval: &E::ScalarField) -> Result<(), PCSError> {
        let poly_commit: PolyCommit<E> = PolyCommit { srs: ck.clone() };
        let (x, y) = ck.split_between_x_and_y(r);

        // verify the proof
        assert!(poly_commit.verify(commitment,proof, x.as_slice(), y.as_slice(), eval));

        Ok(())
    }

    fn setup(mut max_poly_vars: usize, label: &'static [u8], rng: &mut impl RngCore) -> Result<Self::SRS, Error> {
        let x = max_poly_vars / 2;
        let y = max_poly_vars - x;
        let degree_x = 2usize.pow(x as u32);
        let degree_y = 2usize.pow(y as u32);
        Ok(PolyCommit::<E>::setup(degree_x, degree_y, rng))
    }

    fn trim(srs: &Self::SRS, supported_num_vars: usize) -> PCSKeys<E, Self> {
        PCSKeys {
            ck: srs.clone(),
            vk: srs.clone(),
        }
    }
}

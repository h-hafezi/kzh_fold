use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_poly_commit::Error;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::RngCore;
use core::fmt::Debug;

use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::AppendToTranscript;

pub mod error;
pub mod kzh;

pub trait ToAffine<E: Pairing> {
    fn to_affine(self) -> E::G1Affine;
}

pub trait VectorCommitmentScheme<E: Pairing>
where
    <E as Pairing>::ScalarField: Absorb,
{
    type VectorCommitment: AppendToTranscript<E::ScalarField>
    + Sized
    + Sync
    + CanonicalSerialize
    + CanonicalDeserialize;

    type CommitmentKey;

    fn commit(vec: &[E::ScalarField], ck: &Self::CommitmentKey) -> Self::VectorCommitment;
}

pub trait PolyCommitmentScheme<E: Pairing>: Send + Sync
where
    <E as Pairing>::ScalarField: Absorb,
{
    type SRS: CanonicalSerialize + CanonicalDeserialize + Clone;

    type Commitment: AppendToTranscript<E::ScalarField>
    + Debug
    + CanonicalSerialize
    + CanonicalDeserialize
    + PartialEq
    + Eq
    + Clone
    + AppendToTranscript<E::ScalarField>
    + ToAffine<E>;

    // The commitments should be compatible with a homomorphic vector commitment valued in G
    type PolyCommitmentProof: Sync + CanonicalSerialize + CanonicalDeserialize + Debug;

    // Optionally takes `vector_comm` as a "hint" to speed up the commitment process if a
    // commitment to the vector of evaluations has already been computed
    fn commit(
        poly: &MultilinearPolynomial<E::ScalarField>,
        srs: &Self::SRS,
    ) -> Self::Commitment;

    fn prove(
        C: Option<&Self::Commitment>,
        poly: &MultilinearPolynomial<E::ScalarField>,
        r: &[E::ScalarField],
        srs: &Self::SRS,
    ) -> Self::PolyCommitmentProof;

    fn verify(
        commitment: &Self::Commitment,
        proof: &Self::PolyCommitmentProof,
        srs: &Self::SRS,
        r: &[E::ScalarField],
        eval: &E::ScalarField,
    ) -> Result<(), error::PCSError>;

    // Generate a SRS using the provided RNG; this is just for testing purposes, since in reality
    // we need to perform a trusted setup ceremony and then read the SRS from a file.
    fn setup(
        max_poly_vars: usize,
        rng: &mut impl RngCore,
    ) -> Result<Self::SRS, Error>;
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> VectorCommitmentScheme<E> for PC
where
    <E as Pairing>::ScalarField: Absorb,
{
    type VectorCommitment = PC::Commitment;
    type CommitmentKey = PC::SRS;
    fn commit(vec: &[<E>::ScalarField], srs: &Self::CommitmentKey) -> Self::VectorCommitment {
        let poly = MultilinearPolynomial::new(vec.to_vec());
        PC::commit(&poly, srs)
    }
}
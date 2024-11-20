use std::fmt::Debug;
use ark_crypto_primitives::sponge::Absorb;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use rand::Rng;
use crate::nexus_spartan::commitment_traits::ToAffine;
use crate::transcript::transcript::AppendToTranscript;

pub mod kzh2;

pub mod kzh3;

pub mod kzh4;

pub trait KZH<E: Pairing> where <E as Pairing>::ScalarField: Absorb {
    type Degree;
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
    type Opening: Sync + CanonicalSerialize + CanonicalDeserialize + Debug;

    fn split_input(srs: &Self::SRS, input: &[E::ScalarField]) -> Vec<Vec<E::ScalarField>>;

    fn get_degree_from_maximum_supported_degree(n: usize) -> Self::Degree;

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS;

    fn commit(
        srs: &Self::SRS,
        poly: &MultilinearPolynomial<E::ScalarField>,
    ) -> Self::Commitment;

    fn open(
        srs: &Self::SRS,
        input: &[E::ScalarField],
        com: &Self::Commitment,
        poly: &MultilinearPolynomial<E::ScalarField>,
    ) -> Self::Opening;

    fn verify(
        srs: &Self::SRS,
        input: &[E::ScalarField],
        output: &E::ScalarField,
        com: &Self::Commitment,
        open: &Self::Opening,
    );
}

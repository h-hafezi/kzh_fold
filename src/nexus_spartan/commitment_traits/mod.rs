use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_poly_commit::Error;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::RngCore;
use core::fmt::Debug;
use crate::kzh::KZH;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::AppendToTranscript;
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

impl<E: Pairing, PC: KZH<E>> VectorCommitmentScheme<E> for PC
where
    <E as Pairing>::ScalarField: Absorb,
{
    type VectorCommitment = PC::Commitment;
    type CommitmentKey = PC::SRS;
    fn commit(vec: &[<E>::ScalarField], srs: &Self::CommitmentKey) -> Self::VectorCommitment {
        let poly = MultilinearPolynomial::new(vec.to_vec());
        PC::commit(srs, &poly)
    }
}
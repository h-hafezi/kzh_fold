use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use rand::thread_rng;
use crate::kzh::KZH;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::AppendToTranscript;

/// Convert to affine representation. Use `&self` (non-consuming).
pub trait ToAffine<E: Pairing> {
    fn to_affine(&self) -> E::G1Affine;
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

    type VectorAux: AppendToTranscript<E::ScalarField>
    + Sized
    + Sync
    + CanonicalSerialize
    + CanonicalDeserialize;

    type CommitmentKey;

    /// Commit returns both the commitment and its auxiliary/opening information.
    fn commit(
        vec: &[E::ScalarField],
        ck: &Self::CommitmentKey,
    ) -> (Self::VectorCommitment, Self::VectorAux);
}

impl<E: Pairing, PC: KZH<E>> VectorCommitmentScheme<E> for PC
where
    <E as Pairing>::ScalarField: Absorb,
{
    type VectorCommitment = PC::Commitment;
    type VectorAux = PC::Aux;
    type CommitmentKey = PC::SRS;

    fn commit(
        vec: &[<E>::ScalarField],
        srs: &Self::CommitmentKey,
    ) -> (Self::VectorCommitment, Self::VectorAux) {
        let poly = MultilinearPolynomial::new(vec.to_vec());
        // KZH::commit returns (Commitment, Aux)
        PC::commit(srs, &poly, &mut thread_rng())
    }
}

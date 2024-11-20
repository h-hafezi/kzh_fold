use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use rand::Rng;

pub mod kzh2;

pub mod kzh3;

pub mod kzh4;

pub trait KZH<E: Pairing> {
    type Degree;
    type SRS;
    type Commitment;
    type Opening;

    fn split_input(srs: &Self::SRS, input: &[E::ScalarField]) -> Vec<Vec<E::ScalarField>>;

    fn setup<R: Rng>(degree: &Self::Degree, rng: &mut R) -> Self::SRS;

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

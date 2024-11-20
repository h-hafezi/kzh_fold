use ark_ec::pairing::Pairing;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

pub mod kzh2;

pub mod kzh3;

pub mod kzh4;

/*pub trait KZH<E: Pairing>{
    type Degree;
    type Input;
    type ProverKey;
    type VerifierKey;
    type Commitment;
    type Opening;

    fn split_input(input: &[E::ScalarField]) -> Self::Input;

    fn setup(degree: &Self::Degree) -> Self::SRS;

    fn commit(degree: &Self::Degree, poly: &MultilinearPolynomial<E::ScalarField>) -> Self::Commitment;

    fn open(srs: &Self::SRS, input: &Self::Input, poly: &MultilinearPolynomial<E::ScalarField>) -> Self::Opening;

    fn verify(srs: &Self::SRS, input: &Self::Input, output: &E::ScalarField,poly: &MultilinearPolynomial<E::ScalarField>);
}
 */
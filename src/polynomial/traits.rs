use ark_ec::pairing::Pairing;

use crate::polynomial::multilinear_polynomial::bivariate_multilinear::BivariateMultiLinearPolynomial;
use crate::polynomial::multilinear_polynomial::multilinear_poly::MultilinearPolynomial;

/// trait defined for lagrange_base and EqPolynomial
pub trait Evaluable<E: Pairing> {
    type Input;

    fn evaluate(&self, point: &Self::Input) -> Vec<E::ScalarField>;
}

/// trait define for Univariate and MultiLinear polynomials
pub trait OneDimensionalPolynomial<E: Pairing> {
    type Input;

    fn evaluate(&self, input: &Self::Input) -> E::ScalarField;
    fn evaluations_over_boolean_domain(&self) -> Vec<E::ScalarField>;
    fn from_multilinear_polynomial(multi_poly: MultilinearPolynomial<E::ScalarField, E>) -> Self;
}

/// trait define for Bivariate and BivariateMultiLinear polynomials
pub trait TwoDimensionalPolynomial<E: Pairing> {
    type Input;
    type PartialEvalType: OneDimensionalPolynomial<E>;

    fn partial_evaluation(&self, input: &Self::Input) -> Self::PartialEvalType;
    fn partial_evaluations_over_boolean_domain(&self, i: usize) -> Vec<E::ScalarField>;
    fn from_bivariate_multilinear_polynomial(multi_poly: BivariateMultiLinearPolynomial<E::ScalarField, E>) -> Self;
}


pub trait FromPartialEvaluation<E: Pairing, P> {
    fn from_partial_evaluation(evaluation: P) -> Self;
}

impl<E: Pairing> FromPartialEvaluation<E, MultilinearPolynomial<E::ScalarField, E>> for MultilinearPolynomial<E::ScalarField, E> {
    fn from_partial_evaluation(evaluation: MultilinearPolynomial<E::ScalarField, E>) -> Self {
        evaluation
    }
}

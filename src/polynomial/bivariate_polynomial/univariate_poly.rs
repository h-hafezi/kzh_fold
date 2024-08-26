use std::fmt;
use std::fmt::Display;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign};

use ark_ec::pairing::Pairing;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::CanonicalSerialize;

use crate::polynomial::bivariate_polynomial::lagrange_basis::{Evaluatable, LagrangeBasis};

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct UnivariatePolynomial<F: FftField, E: Pairing> {
    pub evaluations: Vec<F>,
    pub lagrange_basis: LagrangeBasis<F>,
    pub phantom: PhantomData<E>,
}

impl<F: FftField, E: Pairing> Display for UnivariatePolynomial<F, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "UnivariatePolynomial {{ evaluations: {:?} }}", self.evaluations)
    }
}

impl<F: FftField, E: Pairing<ScalarField=F>> UnivariatePolynomial<F, E> {
    /// Evaluate the polynomial at p(z) = L_1(z) * p(w_1) + ... + L_n(z) * p(w_n)
    /// Where L_i(z) = Z_w(z) / z - w_i
    pub fn evaluate(&self, z: &F) -> F {
        // the evaluation points p(w_i)
        let w_i = &self.evaluations;
        // the lagrange basis L_i(z)
        let l_i = <LagrangeBasis<F> as Evaluatable<E>>::evaluate(&self.lagrange_basis, *z);
        w_i.iter()
            .zip(l_i.iter())
            .map(|(&a, &b)| a * b)
            .sum()
    }

    pub fn sum_evaluations_in_domain(&self) -> F {
        self.evaluations.iter().cloned().sum()
    }

    #[inline]
    pub fn new(evaluations: Vec<F>, domain: GeneralEvaluationDomain<F>) -> Self {
        Self {
            evaluations,
            lagrange_basis: LagrangeBasis { domain },
            phantom: Default::default(),
        }
    }
}

impl<F: FftField, E: Pairing> Add for UnivariatePolynomial<F, E> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        // Ensure the domains are equal
        assert_eq!(self.lagrange_basis, other.lagrange_basis, "Domains must be equal");

        // Ensure the evaluations vectors are of the same size
        assert_eq!(self.evaluations.len(), other.evaluations.len(), "Evaluation vectors must be of the same length");

        // Perform entry-wise addition of the evaluation vectors
        let evaluations = self
            .evaluations
            .iter()
            .zip(other.evaluations.iter())
            .map(|(a, b)| *a + *b)
            .collect();

        // Return the resulting polynomial
        UnivariatePolynomial {
            evaluations,
            lagrange_basis: self.lagrange_basis,
            phantom: Default::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::pairing::Pairing;
    use ark_ff::{Field, MontBackend};
    use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{E, ScalarField};

    use super::*;

    type F = ScalarField;

    #[test]
    pub fn test_add() {
        let poly_degree = 9usize;
        let lagrange_basis = LagrangeBasis::new(poly_degree);
        let poly1: UnivariatePolynomial<ark_ff::Fp<MontBackend<ark_bn254::FrConfig, 4>, 4>, E> = UnivariatePolynomial {
            evaluations: vec![F::ONE; poly_degree],
            lagrange_basis: lagrange_basis.clone(),
            phantom: Default::default(),
        };
        let poly2: UnivariatePolynomial<ark_ff::Fp<MontBackend<ark_bn254::FrConfig, 4>, 4>, E> = UnivariatePolynomial {
            evaluations: vec![F::ONE; poly_degree],
            lagrange_basis: lagrange_basis.clone(),
            phantom: Default::default(),
        };
        println!("{}", poly1 + poly2);
    }

    #[test]
    pub fn test_evaluation() {
        let poly_degree = 9usize;
        // generate a random polynomial
        let poly = DensePolynomial::rand(poly_degree, &mut thread_rng());
        let z = F::rand(&mut thread_rng());
        // convert the random polynomial into a univariate polynomial
        let mut evaluations = vec![];
        let domain = GeneralEvaluationDomain::<F>::new(poly_degree).unwrap();
        for w_i in domain.elements() {
            let eval = poly.evaluate(&w_i);
            evaluations.push(eval);
        }
        let univariate: UnivariatePolynomial<ark_ff::Fp<MontBackend<ark_bn254::FrConfig, 4>, 4>, E> = UnivariatePolynomial::new(evaluations, domain);
        // assert equality of evaluation of both polynomials at a random point
        assert_eq!(univariate.evaluate(&z), poly.evaluate(&z));
    }
}

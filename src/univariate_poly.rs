/// UnivariatePolynomial is a shallow wrapper on top of DensePolynomial.
/// Might be better to delete this wrapper and just use DensePolynomial directly

use std::ops::Add;

use ark_ff::{Field, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_serialize::CanonicalSerialize;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct UnivariatePolynomial<F: Field> {
    pub poly: DensePolynomial<F>,
}

/// Implement addition of two polynomials
impl<F: Field> Add for UnivariatePolynomial<F> {
    type Output = UnivariatePolynomial<F>;

    fn add(self, other: UnivariatePolynomial<F>) -> Self {
        UnivariatePolynomial { poly: &self.poly + &other.poly }
    }
}

impl<F: Field> UnivariatePolynomial<F> {
    #[inline]
    pub fn new(coeffs: Vec<F>) -> Self {
        Self {
            poly: DensePolynomial {
                coeffs
            }
        }
    }

    /// Evaluate the polynomial at z
    pub fn evaluate(&self, z: &F) -> F {
        self.poly.evaluate(z)
    }
}

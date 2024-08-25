use ark_ff::PrimeField;
use crate::polynomial::multilinear_polynomial::decimal_to_boolean_vector;
use crate::polynomial::multilinear_polynomial::dense_multilinear_poly::MultilinearPolynomial;
use crate::polynomial::multilinear_polynomial::math::Math;

pub struct BivariateMultiLinearPolynomial<F: PrimeField> {
    pub poly: MultilinearPolynomial<F>,
    pub partial_multilinear: Vec<MultilinearPolynomial<F>>,
}

impl<F: PrimeField> BivariateMultiLinearPolynomial<F> {
    // this would output a bivariate multilinear polynomial consisting of n partial evaluations
    pub fn from_multilinear_to_bivariate_multilinear(poly: MultilinearPolynomial<F>, n: usize) -> BivariateMultiLinearPolynomial<F> {
        BivariateMultiLinearPolynomial {
            poly: poly.clone(),
            partial_multilinear: {
                let mut res = Vec::new();
                for i in 0..n {
                    let bits = decimal_to_boolean_vector(i, n.log_2());
                    res.push(poly.partial_evaluation(bits.as_slice()));
                }
                res
            },
        }
    }
}
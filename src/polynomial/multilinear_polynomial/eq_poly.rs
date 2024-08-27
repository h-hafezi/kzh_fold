use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, Zero};

use crate::polynomial::multilinear_polynomial::math::Math;
use crate::polynomial::traits::Evaluable;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EqPolynomial<F: Field + Copy> {
    w: Vec<F>,
}

impl<F: Field + Copy> EqPolynomial<F> {
    /// Creates a new EqPolynomial from a vector `w`
    pub fn new(w: Vec<F>) -> Self {
        EqPolynomial { w }
    }

    /// Evaluates the polynomial eq_w(r) = prod_{i} (w_i * r_i + (F::ONE - w_i) * (F::ONE - r_i))
    pub fn evaluate_at_single_point(&self, r: &[F]) -> F {
        assert_eq!(self.w.len(), r.len(), "Vectors w and r must be of the same length.");

        let mut product = F::ONE;
        for (w_i, r_i) in self.w.iter().zip(r.iter()) {
            let term = (*w_i * *r_i) + ((F::ONE - *w_i) * (F::ONE - *r_i));
            product *= term;
        }
        product
    }
}


impl<E: Pairing> Evaluable<E> for EqPolynomial<E::ScalarField> {
    type Input = Vec<E::ScalarField>;

    fn evaluate(&self, r: &Self::Input) -> Vec<E::ScalarField> {
        let n = r.len();
        let mut dp = vec![vec![E::ScalarField::zero(); 1 << n]; n + 1];
        dp[0][0] = E::ScalarField::one();

        for i in 0..n {
            for j in 0..(1 << i) {
                dp[i + 1][j] = dp[i][j] * (E::ScalarField::one() - r[i]);
                dp[i + 1][j + (1 << i)] = dp[i][j] * r[i];
            }
        }

        dp[n].clone()
    }
}


#[cfg(test)]
mod tests {
    use ark_ff::AdditiveGroup;

    use crate::constant_for_curves::{E, ScalarField};

    use super::*;

    type F = ScalarField;

    #[test]
    fn test() {
        // test for single point evaluation
        let w = vec![F::ZERO, F::ONE, F::ZERO];
        let r = vec![F::ONE, F::ZERO, F::ONE];
        let eq_poly = EqPolynomial::new(w);
        let result = eq_poly.evaluate_at_single_point(&r);
        assert_eq!(result, F::ZERO);

        // test for single point evaluation
        let w = vec![F::ZERO, F::ONE, F::ZERO];
        let r = vec![F::ZERO, F::ONE, F::ZERO];
        let eq_poly = EqPolynomial::new(w);
        let result = eq_poly.evaluate_at_single_point(&r);
        assert_eq!(result, F::ONE);

        // test for range evaluation
        let r = vec![F::ONE, F::ZERO, F::ONE];
        let eq_poly = EqPolynomial::<F>::new(vec![F::ONE]);
        let results: Vec<F> = <EqPolynomial<F> as Evaluable<E>>::evaluate(&eq_poly, &r);
        assert_eq!(results.len(), 8);
        assert_eq!(results, vec![F::ZERO, F::ZERO, F::ZERO, F::ZERO,
                                 F::ZERO, F::ONE, F::ZERO, F::ZERO]);
    }
}

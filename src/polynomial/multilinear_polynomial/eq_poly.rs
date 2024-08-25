use ark_ff::Field;
use crate::polynomial::multilinear_polynomial::math::Math;

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
    pub fn evaluate(&self, r: &[F]) -> F {
        assert_eq!(self.w.len(), r.len(), "Vectors w and r must be of the same length.");

        let mut product = F::ONE;

        for (w_i, r_i) in self.w.iter().zip(r.iter()) {
            let term = (*w_i * *r_i) + ((F::ONE - *w_i) * (F::ONE - *r_i));
            product *= term;
        }

        product
    }

    /// Evaluates the polynomial for all w_i in {0, 1}^n where |r| = n
    /// Efficiently compute eq_w(r) for all w in {0, 1}^n using dynamic programming.
    pub fn evals(r: &[F]) -> Vec<F> {
        let n = r.len();
        let mut dp = vec![vec![F::zero(); 1 << n]; n + 1];
        dp[0][0] = F::one();

        for i in 0..n {
            for j in 0..(1 << i) {
                dp[i + 1][j] = dp[i][j] * (F::one() - r[i]);
                dp[i + 1][j + (1 << i)] = dp[i][j] * r[i];
            }
        }

        dp[n].clone()
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::AdditiveGroup;

    use crate::constant_for_curves::ScalarField;

    use super::*;

    type F = ScalarField;

    #[test]
    fn test_evaluate() {
        let w = vec![F::ZERO, F::ONE, F::ZERO];
        let r = vec![F::ONE, F::ZERO, F::ONE];
        let eq_poly = EqPolynomial::new(w);
        let result = eq_poly.evaluate(&r);
        assert_eq!(result, F::ZERO);

        let w = vec![F::ZERO, F::ONE, F::ZERO];
        let r = vec![F::ZERO, F::ONE, F::ZERO];
        let eq_poly = EqPolynomial::new(w);
        let result = eq_poly.evaluate(&r);
        assert_eq!(result, F::ONE);
    }

    #[test]
    fn test_evals() {
        let r = vec![F::ONE, F::ZERO, F::ONE];
        let results = EqPolynomial::evals(&r);
        assert_eq!(results.len(), 8);
        assert_eq!(results, vec![F::ZERO, F::ZERO, F::ZERO, F::ZERO,
                                 F::ZERO, F::ONE, F::ZERO, F::ZERO]);
    }
}

// borrowed from Arkworks
use ark_ff::{Field};
use crate::math::Math;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EqPolynomial<F: Field + Copy> {
    r: Vec<F>,
}

impl<F: Field + Copy> EqPolynomial<F> {
    /// Creates a new EqPolynomial from a vector `w`
    pub fn new(r: Vec<F>) -> Self {
        EqPolynomial { r }
    }

    /// Evaluates the polynomial eq_w(r) = prod_{i} (w_i * r_i + (F::ONE - w_i) * (F::ONE - r_i))
    pub fn evaluate(&self, rx: &[F]) -> F {
        assert_eq!(self.r.len(), rx.len());
        (0..rx.len())
            .map(|i| self.r[i] * rx[i] + (F::one() - self.r[i]) * (F::one() - rx[i]))
            .product()
    }

    pub fn evals(&self) -> Vec<F> {
        let ell = self.r.len();

        let mut evals: Vec<F> = vec![F::one(); ell.pow2()];
        let mut size = 1;
        for j in 0..ell {
            // in each iteration, we double the size of chis
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2];
                evals[i] = scalar * self.r[j];
                evals[i - 1] = scalar - evals[i];
            }
        }
        evals
    }

    pub fn compute_factored_lens(ell: usize) -> (usize, usize) {
        (ell / 2, ell - ell / 2)
    }

    pub fn compute_factored_evals(&self) -> (Vec<F>, Vec<F>) {
        let ell = self.r.len();
        let (left_num_vars, _right_num_vars) = Self::compute_factored_lens(ell);

        let L = EqPolynomial::new(self.r[..left_num_vars].to_vec()).evals();
        let R = EqPolynomial::new(self.r[left_num_vars..ell].to_vec()).evals();

        (L, R)
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::AdditiveGroup;

    use crate::constant_for_curves::{ScalarField};

    use super::*;

    type F = ScalarField;

    #[test]
    fn test() {
        // test for single point evaluation
        let w = vec![F::ZERO, F::ONE, F::ZERO];
        let r = vec![F::ONE, F::ZERO, F::ONE];
        let eq_poly = EqPolynomial::new(w);
        let result = eq_poly.evaluate(&r);
        assert_eq!(result, F::ZERO);

        // test for single point evaluation
        let w = vec![F::ZERO, F::ONE, F::ZERO];
        let r = vec![F::ZERO, F::ONE, F::ZERO];
        let eq_poly = EqPolynomial::new(w);
        let result = eq_poly.evaluate(&r);
        assert_eq!(result, F::ONE);

        // test for range evaluation
        let r = vec![F::ONE, F::ZERO, F::ONE];
        let results: Vec<F> = <EqPolynomial<F>>::new(r).evals();
        assert_eq!(results.len(), 8);
        assert_eq!(results, vec![F::ZERO, F::ZERO, F::ZERO, F::ZERO,
                                 F::ZERO, F::ONE, F::ZERO, F::ZERO]);
    }
}

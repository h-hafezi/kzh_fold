use ark_ff::FftField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::CanonicalSerialize;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct LagrangeBasis<F: FftField> {
    pub domain: GeneralEvaluationDomain<F>,
}

pub trait LagrangeTraits<F: FftField> {
    fn evaluate(&self, z: &F) -> Vec<F>;

    fn evaluate_vanishing_polynomial(&self, z: &F) -> F;
}

impl<F: FftField> LagrangeTraits<F> for LagrangeBasis<F> {
    // TODO: optimize
    fn evaluate(&self, z: &F) -> Vec<F> {
        let mut evaluation_points = vec![];
        let eval = self.domain.evaluate_vanishing_polynomial(z.clone());

        for w_i in self.domain.elements() {
            if z == &w_i {
                // If z is one of the roots of unity, L_i(z) = 1 if z = w_i, otherwise 0
                evaluation_points.push(F::one());
            } else {
                // L_i(z) = w_i * eval / (z - w_i)
                evaluation_points.push((self.domain.size_inv() * w_i * eval) / (z.clone() - w_i));
            }
        }
        evaluation_points
    }

    fn evaluate_vanishing_polynomial(&self, z: &F) -> F {
        self.domain.evaluate_vanishing_polynomial(z.clone())
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::Field;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use crate::constant_for_curves::ScalarField;
    use crate::polynomial::lagrange_basis::{LagrangeBasis, LagrangeTraits};

    type F = ScalarField;

    #[test]
    fn lagrange_test() {
        let lagrange_basis = LagrangeBasis { domain: GeneralEvaluationDomain::<F>::new(10).unwrap() };
        assert_eq!(lagrange_basis.evaluate(&F::from(2u8)).len(), 16);
    }
}



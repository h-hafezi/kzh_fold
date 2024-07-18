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
            assert_ne!(z.clone(), w_i, "the value z is in the unity roots");
            // L_i(z)
            evaluation_points.push((self.domain.size_inv() * w_i * eval) / (z.clone() - w_i));
        }
        evaluation_points
    }

    fn evaluate_vanishing_polynomial(&self, z: &F) -> F {
        self.domain.evaluate_vanishing_polynomial(z.clone())
    }
}


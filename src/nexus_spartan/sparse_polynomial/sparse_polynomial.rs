use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use crate::math::Math;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparsePoly<F: Absorb> {
    pub num_vars: usize,
    pub evaluations: Vec<F>,
}

impl<F: PrimeField + Absorb> SparsePoly<F> {
    pub fn new(num_vars: usize, evaluations: Vec<F>) -> Self {
        SparsePoly { num_vars, evaluations }
    }

    pub(crate) fn compute_chi(a: &[bool], r: &[F]) -> F {
        assert_eq!(a.len(), r.len());
        let mut chi_i = F::one();
        for j in 0..r.len() {
            if a[j] {
                chi_i *= r[j];
            } else {
                chi_i *= F::one() - r[j];
            }
        }
        chi_i
    }

    pub fn evaluate(&self, r: &[F]) -> F {
        assert_eq!(self.num_vars, r.len());

        (0..self.evaluations.len())
            .map(|i| {
                let bits = i.get_bits_canonical_order(r.len());
                SparsePoly::compute_chi(&bits, r) * self.evaluations[i]
            })
            .sum()
    }
}

// Hossein
// Stolen from: https://github.com/nexus-xyz/nexus-zkvm/blob/main/spartan/src/dense_mlpoly.rs
// The bound_* functions will come useful when you want to partially evaluate something

#![allow(clippy::too_many_arguments)]
use ark_ec::CurveGroup;
use ark_ec::VariableBaseMSM;
use ark_ff::PrimeField;
use ark_serialize::*;
use ark_std::Zero;
use core::ops::Index;

use crate::spartan::math::Math;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

fn compute_dotproduct<F: PrimeField>(a: &[F], b: &[F]) -> F {
    assert_eq!(a.len(), b.len());
    (0..a.len()).map(|i| a[i] * b[i]).sum()
}

#[derive(Debug, Clone, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct DensePolynomial<F>
where
  F: Sync + CanonicalSerialize + CanonicalDeserialize,
{
  num_vars: usize, // the number of variables in the multilinear polynomial
  len: usize,
  Z: Vec<F>, // evaluations of the polynomial in all the 2^num_vars Boolean inputs
}

pub struct EqPolynomial<F> {
  r: Vec<F>,
}

impl<F: PrimeField> EqPolynomial<F> {
  pub fn new(r: Vec<F>) -> Self {
    EqPolynomial { r }
  }

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

pub struct IdentityPolynomial {
  size_point: usize,
}

impl IdentityPolynomial {
  pub fn new(size_point: usize) -> Self {
    IdentityPolynomial { size_point }
  }

  pub fn evaluate<F: PrimeField>(&self, r: &[F]) -> F {
    let len = r.len();
    assert_eq!(len, self.size_point);
    (0..len)
      .map(|i| F::from((len - i - 1).pow2() as u64) * r[i])
      .sum()
  }
}

impl<F: PrimeField> DensePolynomial<F> {
  pub fn new(Z: Vec<F>) -> Self {
    DensePolynomial {
      num_vars: Z.len().log_2(),
      len: Z.len(),
      Z,
    }
  }

  pub fn get_num_vars(&self) -> usize {
    self.num_vars
  }

  pub fn len(&self) -> usize {
    self.len
  }

  pub fn is_empty(&self) -> bool {
    self.len == 0
  }

  pub fn split(&self, idx: usize) -> (Self, Self) {
    assert!(idx < self.len());
    (
      Self::new(self.Z[..idx].to_vec()),
      Self::new(self.Z[idx..2 * idx].to_vec()),
    )
  }

  pub fn bound(&self, L: &[F]) -> Vec<F> {
    let (left_num_vars, right_num_vars) =
      EqPolynomial::<F>::compute_factored_lens(self.get_num_vars());
    let L_size = left_num_vars.pow2();
    let R_size = right_num_vars.pow2();
    (0..R_size)
      .map(|i| (0..L_size).map(|j| L[j] * self.Z[j * R_size + i]).sum())
      .collect()
  }

  pub fn bound_poly_var_top(&mut self, r: &F) {
    let n = self.len() / 2;
    for i in 0..n {
      self.Z[i] = self.Z[i] + *r * (self.Z[i + n] - self.Z[i]);
    }
    self.num_vars -= 1;
    self.len = n;
  }

  pub fn bound_poly_var_bot(&mut self, r: &F) {
    let n = self.len() / 2;
    for i in 0..n {
      self.Z[i] = self.Z[2 * i] + *r * (self.Z[2 * i + 1] - self.Z[2 * i]);
    }
    self.num_vars -= 1;
    self.len = n;
  }

  // returns Z(r) in O(n) time
  pub fn evaluate(&self, r: &[F]) -> F
  {
    // r must have a value for each variable
    assert_eq!(r.len(), self.get_num_vars());
    let chis = EqPolynomial::new(r.to_vec()).evals();
    assert_eq!(chis.len(), self.Z.len());
    compute_dotproduct(&self.Z, &chis)
  }

  pub fn vec(&self) -> &Vec<F> {
    &self.Z
  }

  pub fn extend(&mut self, other: &DensePolynomial<F>) {
    // TODO: allow extension even when some vars are bound
    assert_eq!(self.Z.len(), self.len);
    let other_vec = other.vec();
    assert_eq!(other_vec.len(), self.len);
    self.Z.extend(other_vec);
    self.num_vars += 1;
    self.len *= 2;
    assert_eq!(self.Z.len(), self.len);
  }

  pub fn merge(polys: &[DensePolynomial<F>]) -> DensePolynomial<F> {
    let mut Z: Vec<F> = Vec::new();
    for poly in polys.iter() {
      Z.extend(poly.vec().iter());
    }

    // pad the polynomial with zero polynomial at the end
    Z.resize(Z.len().next_power_of_two(), F::zero());

    DensePolynomial::new(Z)
  }

  pub fn from_usize(Z: &[usize]) -> Self {
    DensePolynomial::new(
      (0..Z.len())
        .map(|i| F::from(Z[i] as u64))
        .collect::<Vec<F>>(),
    )
  }
}

impl<F> Index<usize> for DensePolynomial<F>
where
  F: Sync + CanonicalDeserialize + CanonicalSerialize,
{
  type Output = F;

  #[inline(always)]
  fn index(&self, _index: usize) -> &F {
    &(self.Z[_index])
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use ark_bls12_381::Fr;
  use ark_bls12_381::G1Projective;
  use ark_std::test_rng;
  use ark_std::One;
  use ark_std::UniformRand;

  fn evaluate_with_LR<G: CurveGroup>(Z: &[G::ScalarField], r: &[G::ScalarField]) -> G::ScalarField {
    let eq = EqPolynomial::<G::ScalarField>::new(r.to_vec());
    let (L, R) = eq.compute_factored_evals();

    let ell = r.len();
    // ensure ell is even
    assert!(ell % 2 == 0);
    // compute n = 2^\ell
    let n = ell.pow2();
    // compute m = sqrt(n) = 2^{\ell/2}
    let m = n.square_root();

    // compute vector-matrix product between L and Z viewed as a matrix
    let LZ = (0..m)
      .map(|i| (0..m).map(|j| L[j] * Z[j * m + i]).sum())
      .collect::<Vec<G::ScalarField>>();

    // compute dot product between LZ and R
    compute_dotproduct(&LZ, &R)
  }

  #[test]
  fn check_polynomial_evaluation() {
    check_polynomial_evaluation_helper::<G1Projective>()
  }

  fn check_polynomial_evaluation_helper<G: CurveGroup>() {
    // Z = [1, 2, 1, 4]
    let Z = vec![
      G::ScalarField::one(),
      G::ScalarField::from(2u64),
      G::ScalarField::one(),
      G::ScalarField::from(4u64),
    ];

    // r = [4,3]
    let r = vec![G::ScalarField::from(4u64), G::ScalarField::from(3u64)];

    let eval_with_LR = evaluate_with_LR::<G>(&Z, &r);
    let poly = DensePolynomial::new(Z);

    let eval = poly.evaluate(&r);
    assert_eq!(eval, G::ScalarField::from(28u64));
    assert_eq!(eval_with_LR, eval);
  }

  pub fn compute_factored_chis_at_r<F: PrimeField>(r: &[F]) -> (Vec<F>, Vec<F>) {
    let mut L: Vec<F> = Vec::new();
    let mut R: Vec<F> = Vec::new();

    let ell = r.len();
    assert!(ell % 2 == 0); // ensure ell is even
    let n = ell.pow2();
    let m = n.square_root();

    // compute row vector L
    for i in 0..m {
      let mut chi_i = F::one();
      for j in 0..ell / 2 {
        let bit_j = ((m * i) & (1 << (r.len() - j - 1))) > 0;
        if bit_j {
          chi_i *= r[j];
        } else {
          chi_i *= F::one() - r[j];
        }
      }
      L.push(chi_i);
    }

    // compute column vector R
    for i in 0..m {
      let mut chi_i = F::one();
      for j in ell / 2..ell {
        let bit_j = (i & (1 << (r.len() - j - 1))) > 0;
        if bit_j {
          chi_i *= r[j];
        } else {
          chi_i *= F::one() - r[j];
        }
      }
      R.push(chi_i);
    }
    (L, R)
  }

  pub fn compute_chis_at_r<F: PrimeField>(r: &[F]) -> Vec<F> {
    let ell = r.len();
    let n = ell.pow2();
    let mut chis: Vec<F> = Vec::new();
    for i in 0..n {
      let mut chi_i = F::one();
      for j in 0..r.len() {
        let bit_j = (i & (1 << (r.len() - j - 1))) > 0;
        if bit_j {
          chi_i *= r[j];
        } else {
          chi_i *= F::one() - r[j];
        }
      }
      chis.push(chi_i);
    }
    chis
  }

  pub fn compute_outerproduct<F: PrimeField>(L: Vec<F>, R: Vec<F>) -> Vec<F> {
    assert_eq!(L.len(), R.len());
    (0..L.len())
      .map(|i| (0..R.len()).map(|j| L[i] * R[j]).collect::<Vec<F>>())
      .collect::<Vec<Vec<F>>>()
      .into_iter()
      .flatten()
      .collect::<Vec<F>>()
  }

  #[test]
  fn check_memoized_chis() {
    check_memoized_chis_helper::<G1Projective>()
  }

  fn check_memoized_chis_helper<G: CurveGroup>() {
    let mut prng = test_rng();

    let s = 10;
    let mut r: Vec<G::ScalarField> = Vec::new();
    for _i in 0..s {
      r.push(G::ScalarField::rand(&mut prng));
    }
    let chis = tests::compute_chis_at_r::<G::ScalarField>(&r);
    let chis_m = EqPolynomial::<G::ScalarField>::new(r).evals();
    assert_eq!(chis, chis_m);
  }

  #[test]
  fn check_factored_chis() {
    check_factored_chis_helper::<Fr>()
  }

  fn check_factored_chis_helper<F: PrimeField>() {
    let mut prng = test_rng();

    let s = 10;
    let mut r: Vec<F> = Vec::new();
    for _i in 0..s {
      r.push(F::rand(&mut prng));
    }
    let chis = EqPolynomial::new(r.clone()).evals();
    let (L, R) = EqPolynomial::new(r).compute_factored_evals();
    let O = compute_outerproduct(L, R);
    assert_eq!(chis, O);
  }

  #[test]
  fn check_memoized_factored_chis() {
    check_memoized_factored_chis_helper::<Fr>()
  }

  fn check_memoized_factored_chis_helper<F: PrimeField>() {
    let mut prng = test_rng();

    let s = 10;
    let mut r: Vec<F> = Vec::new();
    for _i in 0..s {
      r.push(F::rand(&mut prng));
    }
    let (L, R) = tests::compute_factored_chis_at_r(&r);
    let eq = EqPolynomial::new(r);
    let (L2, R2) = eq.compute_factored_evals();
    assert_eq!(L, L2);
    assert_eq!(R, R2);
  }
}

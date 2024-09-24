#![allow(clippy::too_many_arguments)]
use core::ops::Index;
use std::fmt::Debug;
use std::marker::PhantomData;

use ark_crypto_primitives::crh::sha256::digest::typenum::private::PrivateIntegerAdd;
use ark_ec::CurveGroup;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ff::{Field, PrimeField};
use ark_serialize::*;
use ark_std::Zero;
use itertools::Itertools;
use num::One;
use rand::{Rng, RngCore};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::polynomial::multilinear_polynomial::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_polynomial::math::Math;
use crate::utils::inner_product;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultilinearPolynomial<F: PrimeField, E: Pairing> {
    /// the number of variables in the multilinear polynomial
    pub(crate) num_variables: usize,
    /// evaluations of the polynomial in all the 2^num_vars Boolean inputs
    pub(crate) evaluation_over_boolean_hypercube: Vec<F>,
    /// length of Z = 2^num_vars
    pub(crate) len: usize,
    /// phantom data for E
    pub phantom: PhantomData<E>,
}

impl<F: PrimeField, E: Pairing> MultilinearPolynomial<F, E>
where
    E: Pairing<ScalarField=F>,
{
    pub fn new(evaluation_over_boolean_hypercube: Vec<F>) -> Self {
        MultilinearPolynomial {
            num_variables: evaluation_over_boolean_hypercube.len().log_2(),
            len: evaluation_over_boolean_hypercube.len(),
            evaluation_over_boolean_hypercube,
            phantom: Default::default(),
        }
    }

    pub fn get_num_vars(&self) -> usize {
        self.num_variables
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn clone(&self) -> Self {
        Self::new(self.evaluation_over_boolean_hypercube[0..self.len].to_vec())
    }

    pub fn split(&self, idx: usize) -> (Self, Self) {
        assert!(idx < self.len());
        (
            Self::new(self.evaluation_over_boolean_hypercube[..idx].to_vec()),
            Self::new(self.evaluation_over_boolean_hypercube[idx..2 * idx].to_vec()),
        )
    }

    pub fn bound(&self, L: &[F]) -> Vec<F> {
        let (left_num_vars, right_num_vars) =
            EqPolynomial::<F>::compute_factored_lens(self.get_num_vars());
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();
        (0..R_size)
            .map(|i| (0..L_size).map(|j| L[j] * self.evaluation_over_boolean_hypercube[j * R_size + i]).sum())
            .collect()
    }

    pub fn bound_poly_var_top(&mut self, r: &F) {
        let n = self.len() / 2;
        for i in 0..n {
            self.evaluation_over_boolean_hypercube[i] = self.evaluation_over_boolean_hypercube[i]
                + *r * (self.evaluation_over_boolean_hypercube[i + n]
                - self.evaluation_over_boolean_hypercube[i]);
        }
        self.num_variables -= 1;
        self.len = n;
        self.evaluation_over_boolean_hypercube = self.evaluation_over_boolean_hypercube[..n].to_vec();
    }

    pub fn bound_poly_var_bot(&mut self, r: &F) {
        let n = self.len() / 2;
        for i in 0..n {
            self.evaluation_over_boolean_hypercube[i] = self.evaluation_over_boolean_hypercube[2 * i]
                + *r * (self.evaluation_over_boolean_hypercube[2 * i + 1]
                - self.evaluation_over_boolean_hypercube[2 * i]);
        }
        self.num_variables -= 1;
        self.len = n;
        self.evaluation_over_boolean_hypercube = self.evaluation_over_boolean_hypercube[..n].to_vec();
    }

    // returns Z(r) in O(n) time
    pub fn evaluate(&self, r: &[F]) -> F
    {
        // r must have a value for each variable
        assert_eq!(r.len(), self.get_num_vars());
        println!("length {} {}", r.len(), self.get_num_vars());
        let chis = EqPolynomial::new(r.to_vec()).evals();
        assert_eq!(chis.len(), self.evaluation_over_boolean_hypercube.len());
        inner_product(&self.evaluation_over_boolean_hypercube, &chis)
    }

    fn vec(&self) -> &Vec<F> {
        &self.evaluation_over_boolean_hypercube
    }

    pub fn extend(&mut self, other: &MultilinearPolynomial<F, E>) {
        // TODO: allow extension even when some vars are bound
        assert_eq!(self.evaluation_over_boolean_hypercube.len(), self.len);
        let other_vec = other.vec();
        assert_eq!(other_vec.len(), self.len);
        self.evaluation_over_boolean_hypercube.extend(other_vec);
        self.num_variables += 1;
        self.len *= 2;
        assert_eq!(self.evaluation_over_boolean_hypercube.len(), self.len);
    }

    pub fn merge(polys: &[MultilinearPolynomial<F, E>]) -> MultilinearPolynomial<F, E> {
        let mut Z: Vec<F> = Vec::new();
        for poly in polys.iter() {
            Z.extend(poly.vec().iter());
        }

        // pad the polynomial with zero polynomial at the end
        Z.resize(Z.len().next_power_of_two(), F::zero());

        MultilinearPolynomial::new(Z)
    }

    pub fn from_usize(Z: &[usize]) -> Self {
        MultilinearPolynomial::new(
            (0..Z.len())
                .map(|i| F::from(Z[i] as u64))
                .collect::<Vec<F>>(),
        )
    }
}

impl<F: PrimeField, E: Pairing<ScalarField=F>> MultilinearPolynomial<F, E> {
    /// Perform partial evaluation by fixing the first `fixed_vars.len()` variables to `fixed_vars`.
    /// Returns a new MultilinearPolynomial in the remaining variables.
    pub fn partial_evaluation(&self, fixed_vars: &[F]) -> MultilinearPolynomial<F, E> {
        let mut temp = self.clone();
        for r in fixed_vars {
            temp.bound_poly_var_top(r);
        }
        temp
    }

    pub fn get_partial_evaluation_for_boolean_input(&self, index: usize, n: usize) -> Vec<E::ScalarField> {
        self.evaluation_over_boolean_hypercube[n * index..n * index +n].to_vec()
    }

    pub fn rand<T: RngCore>(num_variables: usize, rng: &mut T) -> MultilinearPolynomial<F, E> {
        MultilinearPolynomial {
            num_variables,
            evaluation_over_boolean_hypercube: {
                let size = 1 << num_variables;
                let mut vector = Vec::with_capacity(size);
                for _ in 0..size {
                    vector.push(F::rand(rng));
                }
                vector
            },
            len: 1 << num_variables,
            phantom: Default::default(),
        }
    }
}


#[cfg(test)]
mod tests {
    use ark_ff::AdditiveGroup;
    use ark_std::One;
    use ark_std::test_rng;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{E, ScalarField};
    use crate::polynomial::multilinear_polynomial::decimal_to_boolean_vector;
    use crate::polynomial::multilinear_polynomial::eq_poly::EqPolynomial;

    use super::*;

    #[test]
    fn tests_partial_eval() {
        let p = MultilinearPolynomial::<ScalarField, E>::rand(3, &mut thread_rng());
        let r_1 = vec![
            ScalarField::ONE,
            ScalarField::ZERO,
        ];

        let r_2 = vec![
            ScalarField::rand(&mut thread_rng()),
        ];

        let concat = {
            let mut temp = r_1.clone();
            temp.extend(r_2.clone());
            temp
        };

        let p_prime = p.partial_evaluation(r_1.as_slice());
        let val_1 = p_prime.evaluate(r_2.as_slice());
        let val_2 = p.evaluate(concat.as_slice());

        assert_eq!(val_1, val_2);
    }

    fn evaluate_with_LR<E: Pairing>(Z: &[E::ScalarField], r: &[E::ScalarField]) -> E::ScalarField {
        let eq = EqPolynomial::<E::ScalarField>::new(r.to_vec());
        let (L, R) = eq.compute_factored_evals();

        let ell = r.len();
        // ensure ell is even
        assert_eq!(ell % 2, 0);
        // compute n = 2^\ell
        let n = ell.pow2();
        // compute m = sqrt(n) = 2^{\ell/2}
        let m = n.square_root();

        // compute vector-matrix product between L and Z viewed as a matrix
        let LZ = (0..m)
            .map(|i| (0..m).map(|j| L[j] * Z[j * m + i]).sum())
            .collect::<Vec<E::ScalarField>>();

        // compute dot product between LZ and R
        inner_product(&LZ, &R)
    }

    #[test]
    fn check_polynomial_evaluation() {
        check_polynomial_evaluation_helper::<E>()
    }

    fn check_polynomial_evaluation_helper<E: Pairing>() {
        // Z = [1, 2, 1, 4]
        let Z = vec![
            E::ScalarField::one(),
            E::ScalarField::from(2u64),
            E::ScalarField::one(),
            E::ScalarField::from(4u64),
        ];

        // r = [4,3]
        let r = vec![E::ScalarField::from(4u64), E::ScalarField::from(3u64)];

        let eval_with_LR = evaluate_with_LR::<E>(&Z, &r);
        let poly: MultilinearPolynomial<<E as Pairing>::ScalarField, E> = MultilinearPolynomial::new(Z);

        let eval = poly.evaluate(&r);
        assert_eq!(eval, E::ScalarField::from(28u64));
        assert_eq!(eval_with_LR, eval);
    }

    pub fn compute_factored_chis_at_r<F: PrimeField>(r: &[F]) -> (Vec<F>, Vec<F>) {
        let mut L: Vec<F> = Vec::new();
        let mut R: Vec<F> = Vec::new();

        let ell = r.len();
        assert_eq!(ell % 2, 0); // ensure ell is even
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
        check_memoized_chis_helper::<E>()
    }

    fn check_memoized_chis_helper<E: Pairing>() {
        let mut prng = test_rng();

        let s = 10;
        let mut r: Vec<E::ScalarField> = Vec::new();
        for _i in 0..s {
            r.push(E::ScalarField::rand(&mut prng));
        }
        let chis = compute_chis_at_r::<E::ScalarField>(&r);
        let chis_m = EqPolynomial::<E::ScalarField>::new(r).evals();
        assert_eq!(chis, chis_m);
    }

    #[test]
    fn check_factored_chis() {
        check_factored_chis_helper::<ScalarField>()
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
        check_memoized_factored_chis_helper::<ScalarField>()
    }

    fn check_memoized_factored_chis_helper<F: PrimeField>() {
        let mut prng = test_rng();

        let s = 10;
        let mut r: Vec<F> = Vec::new();
        for _i in 0..s {
            r.push(F::rand(&mut prng));
        }
        let (L, R) = compute_factored_chis_at_r(&r);
        let eq = EqPolynomial::new(r);
        let (L2, R2) = eq.compute_factored_evals();
        assert_eq!(L, L2);
        assert_eq!(R, R2);
    }
}
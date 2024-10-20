// mostly borrowed from Arkworks

#![allow(clippy::too_many_arguments)]
use core::ops::{Add, Index};
use std::fmt::Debug;

use ark_ff::{Field, PrimeField};
use ark_serialize::*;
use rand::{Rng, RngCore};

#[cfg(feature = "parallel")]
use crate::polynomial::eq_poly::EqPolynomial;
use crate::math::Math;
use crate::utils::inner_product;

#[derive(Debug, Clone, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
pub struct MultilinearPolynomial<F: PrimeField> {
    /// the number of variables in the multilinear polynomial
    pub(crate) num_variables: usize,
    /// evaluations of the polynomial in all the 2^num_vars Boolean inputs
    pub(crate) evaluation_over_boolean_hypercube: Vec<F>,
    /// length of Z = 2^num_vars
    pub(crate) len: usize,
}

impl<F: PrimeField> MultilinearPolynomial<F> {
    pub fn new(evaluation_over_boolean_hypercube: Vec<F>) -> Self {
        MultilinearPolynomial {
            num_variables: evaluation_over_boolean_hypercube.len().log_2(),
            len: evaluation_over_boolean_hypercube.len(),
            evaluation_over_boolean_hypercube,
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

        let chis = EqPolynomial::new(r.to_vec()).evals();
        assert_eq!(chis.len(), self.evaluation_over_boolean_hypercube.len());

        inner_product(&self.evaluation_over_boolean_hypercube, &chis)
    }

    fn vec(&self) -> &Vec<F> {
        &self.evaluation_over_boolean_hypercube
    }

    pub fn extend(&mut self, other: &MultilinearPolynomial<F>) {
        // TODO: allow extension even when some vars are bound
        assert_eq!(self.evaluation_over_boolean_hypercube.len(), self.len);
        let other_vec = other.vec();
        assert_eq!(other_vec.len(), self.len);
        self.evaluation_over_boolean_hypercube.extend(other_vec);
        self.num_variables += 1;
        self.len *= 2;
        assert_eq!(self.evaluation_over_boolean_hypercube.len(), self.len);
    }

    pub fn merge(polys: &[MultilinearPolynomial<F>]) -> MultilinearPolynomial<F> {
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

    /// Get c(x) = self(x) v other(x)
    pub fn get_bitfield_union_poly(&self, other: &Self) -> Self {
        assert_eq!(self.num_variables, other.num_variables);

        let evaluations: Vec<F> = self.evaluation_over_boolean_hypercube.iter()
            .zip(&other.evaluation_over_boolean_hypercube)
            .map(|(a, b)| {
                *a + *b - *a * *b // Since a, b are either 0 or 1, this is equivalent to a | b
            })
            .collect();

        let len = evaluations.len();

        Self {
            num_variables: self.num_variables,
            evaluation_over_boolean_hypercube: evaluations,
            len: len,
        }
    }

    /// Given f(x), compute r * f(x)
    pub fn scalar_mul(&mut self, r: &F) {
        for f_i in self.evaluation_over_boolean_hypercube.iter_mut() {
            *f_i = *f_i * r;
        }
    }

    /// extend it to some number of variables
    pub fn extend_number_of_variables(&self, num_variables: usize) -> MultilinearPolynomial<F> {
        assert!(self.num_variables <= num_variables);
        let mut temp =self.evaluation_over_boolean_hypercube.clone();
        temp.extend(vec![F::ZERO; (1 << num_variables) - self.evaluation_over_boolean_hypercube.len()]);
        MultilinearPolynomial {
            num_variables,
            evaluation_over_boolean_hypercube: temp,
            len: 1 << num_variables,
        }
    }
}

impl<F: PrimeField> MultilinearPolynomial<F> {
    /// Perform partial evaluation by fixing the first `fixed_vars.len()` variables to `fixed_vars`.
    /// Returns a new MultilinearPolynomial in the remaining variables.
    pub fn partial_evaluation(&self, fixed_vars: &[F]) -> MultilinearPolynomial<F> {
        let mut temp = self.clone();
        for r in fixed_vars {
            temp.bound_poly_var_top(r);
        }
        temp
    }

    pub fn get_partial_evaluation_for_boolean_input(&self, index: usize, n: usize) -> Vec<F> {
        self.evaluation_over_boolean_hypercube[n * index..n * index + n].to_vec()
    }

    pub fn rand<T: RngCore>(num_variables: usize, rng: &mut T) -> MultilinearPolynomial<F> {
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
        }
    }

    /// Return a multilinear poly with evaluations that are either 0 or 1
    pub fn random_binary<T: RngCore>(num_variables: usize, rng: &mut T) -> MultilinearPolynomial<F> {
        let evals_len = 1 << num_variables;

        let evals = (0..evals_len).map(|_| {
            let random_bit = rng.gen_bool(0.5); // Generates a random boolean with equal probability
            if random_bit { F::one() } else { F::zero() }
        }).collect();

        MultilinearPolynomial {
            num_variables,
            evaluation_over_boolean_hypercube: evals,
            len: evals_len,
        }
    }
}

impl<F: PrimeField> Index<usize> for MultilinearPolynomial<F> {
    type Output = F;

    #[inline(always)]
    fn index(&self, _index: usize) -> &F {
        &(self.evaluation_over_boolean_hypercube[_index])
    }
}


impl<F: PrimeField> Add for MultilinearPolynomial<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        // Ensure that both polynomials have the same number of variables
        assert_eq!(self.num_variables, other.num_variables, "Polynomials must have the same number of variables");

        // Ensure both polynomials have evaluation vectors of the same length
        assert_eq!(self.len, other.len, "Evaluation vectors must have the same length");

        // Perform element-wise addition over the evaluation vectors
        let new_evaluation_over_boolean_hypercube: Vec<F> = self.evaluation_over_boolean_hypercube
            .iter()
            .zip(other.evaluation_over_boolean_hypercube.iter())
            .map(|(a, b)| *a + *b)
            .collect();

        // Return a new MultilinearPolynomial with the result
        MultilinearPolynomial {
            num_variables: self.num_variables,
            evaluation_over_boolean_hypercube: new_evaluation_over_boolean_hypercube,
            len: self.len,
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::pairing::Pairing;
use ark_ff::AdditiveGroup;
    use ark_std::One;
    use ark_std::test_rng;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{E, ScalarField};
    use crate::polynomial::eq_poly::EqPolynomial;

    use super::*;

    #[test]
    fn extend_test() {
        // define a polynomial
        let poly = MultilinearPolynomial::<ScalarField>::rand(3, &mut thread_rng());
        // extend it with 2 variables
        let new_poly = poly.extend_number_of_variables(5);
        // define what extension looks like
        let mut temp = poly.evaluation_over_boolean_hypercube.clone();
        temp.extend(vec![ScalarField::ZERO; 32 - 8]);
        // assert the equality
        assert_eq!(temp, new_poly.evaluation_over_boolean_hypercube);
    }

    #[test]
    fn tests_partial_eval() {
        let p = MultilinearPolynomial::<ScalarField>::rand(3, &mut thread_rng());
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
        let poly: MultilinearPolynomial<<E as Pairing>::ScalarField> = MultilinearPolynomial::new(Z);

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

    #[test]
    fn check_extend_function() {
        // random bivariate polynomial
        let polynomial = MultilinearPolynomial::rand(1, &mut thread_rng());

        // random points and evaluation
        let x1 = vec![
            ScalarField::rand(&mut thread_rng()),
        ];

        let z1 = polynomial.evaluate(&x1);

        let polynomial = polynomial.extend_number_of_variables(4);

        // random points and evaluation
        let mut x2 = vec![
            ScalarField::ZERO,
            ScalarField::ZERO,
            ScalarField::ZERO,
        ];
        x2.extend(x1);

        let z2 = polynomial.evaluate(&x2);

        assert_eq!(z1, z2);
    }
}

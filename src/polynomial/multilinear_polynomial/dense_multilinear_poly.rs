// Hossein
// Stolen from: https://github.com/nexus-xyz/nexus-zkvm/blob/main/spartan/src/dense_mlpoly.rs
// The bound_* functions will come useful when you want to partially evaluate something

#![allow(clippy::too_many_arguments)]
use core::ops::Index;
use std::fmt::Debug;

use ark_ec::CurveGroup;
use ark_ec::VariableBaseMSM;
use ark_ff::{Field, PrimeField};
use ark_serialize::*;
use ark_std::Zero;
use itertools::Itertools;
use rand::{Rng, RngCore};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::polynomial::multilinear_polynomial::{boolean_vector_to_decimal, compute_dot_product};
use crate::polynomial::multilinear_polynomial::eq_poly::EqPolynomial;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultilinearPolynomial<F: PrimeField> {
    pub num_variables: usize,

    /// `evaluation_over_boolean_hypercube` represents the evaluation of the multilinear polynomial
    /// over all possible combinations of boolean values for its variables. Specifically, this field
    /// is a vector where each entry corresponds to the polynomial's evaluation at a particular point
    /// in the boolean hypercube.
    ///
    /// For a polynomial with `num_variables` variables, the boolean hypercube consists of all possible
    /// combinations of 0s and 1s for these variables. There are `2^num_variables` such combinations.
    /// The index of an entry in `evaluation_over_boolean_hypercube` corresponds to a particular
    /// combination of boolean variable values, represented as a decimal number. The value at this index
    /// is the result of evaluating the polynomial at the corresponding combination of boolean values.
    ///
    /// For example, consider a polynomial `f(x1, x2, x3) = 10 + x1 + 2 * x2*x3 + 5 * x1*x3`. With 3 variables
    /// (x1, x2, x3), the boolean hypercube has `2^3 = 8` points:
    ///
    /// - `x1 = 0, x2 = 0, x3 = 0` (decimal index 0): `f(0, 0, 0) = 10`
    /// - `x1 = 1, x2 = 0, x3 = 0` (decimal index 1): `f(1, 0, 0) = 11`
    /// - `x1 = 0, x2 = 1, x3 = 0` (decimal index 2): `f(0, 1, 0) = 10`
    /// - `x1 = 1, x2 = 1, x3 = 0` (decimal index 3): `f(1, 1, 0) = 11`
    /// - `x1 = 0, x2 = 0, x3 = 1` (decimal index 4): `f(0, 0, 1) = 10`
    /// - `x1 = 1, x2 = 0, x3 = 1` (decimal index 5): `f(1, 0, 1) = 16`
    /// - `x1 = 0, x2 = 1, x3 = 1` (decimal index 6): `f(0, 1, 1) = 12`
    /// - `x1 = 1, x2 = 1, x3 = 1` (decimal index 7): `f(1, 1, 1) = 18`
    ///
    pub evaluation_over_boolean_hypercube: Vec<F>,
}

impl<F: PrimeField> MultilinearPolynomial<F> {
    pub fn new(num_variables: usize, evaluation_over_boolean_hypercube: Vec<F>) -> Self {
        // Ensure length is 2^num_variables
        assert_eq!(evaluation_over_boolean_hypercube.len(), 1 << num_variables);
        MultilinearPolynomial {
            num_variables,
            evaluation_over_boolean_hypercube,
        }
    }

    pub fn get_num_vars(&self) -> usize {
        self.num_variables
    }

    /// r[i] will represent value of x_i in the polynomial
    pub fn evaluate(&self, r: &[F]) -> F {
        // r must have a value for each variable
        assert_eq!(r.len(), self.get_num_vars());
        let chis = EqPolynomial::evals(r);
        assert_eq!(chis.len(), self.evaluation_over_boolean_hypercube.len());
        compute_dot_product(&self.evaluation_over_boolean_hypercube, &chis)
    }

    pub fn evaluate_binary(&self, r: &[F]) -> F {
        let decimal = boolean_vector_to_decimal(r);
        // Return the corresponding evaluation from the boolean hypercube
        self.evaluation_over_boolean_hypercube[decimal]
    }
}

impl<F: PrimeField> MultilinearPolynomial<F> {
    /// Perform partial evaluation by fixing the first `fixed_vars.len()` variables to `fixed_vars`.
    /// Returns a new MultilinearPolynomial in the remaining variables.
    pub fn partial_evaluation(&self, fixed_vars: &[F]) -> MultilinearPolynomial<F> {
        let fixed_vars_len = fixed_vars.len();
        let remaining_vars_len = self.num_variables - fixed_vars_len;
        assert!(remaining_vars_len > 0);

        // Compute the size of the boolean hypercube for the remaining variables
        let remaining_size = 1 << remaining_vars_len;

        // Initialize the vector for the new evaluations
        let mut new_evaluations = vec![F::zero(); remaining_size];

        // Iterate over all possible combinations of the remaining variables
        for i in 0..remaining_size {
            // Convert `i` to a boolean vector representing the current combination of remaining variables
            let mut remaining_vars_assignment = vec![F::zero(); remaining_vars_len];
            for j in 0..remaining_vars_len {
                remaining_vars_assignment[j] = if (i >> j) & 1 == 1 { F::one() } else { F::zero() };
            }

            // Combine `fixed_vars` and `remaining_vars_assignment` to form the full assignment
            let mut full_assignment = Vec::with_capacity(self.num_variables);
            full_assignment.extend_from_slice(fixed_vars);
            full_assignment.extend_from_slice(&remaining_vars_assignment);

            // Evaluate the polynomial at this full assignment
            new_evaluations[i] = self.evaluate(&full_assignment);
        }

        // Return the new MultilinearPolynomial with the remaining variables
        MultilinearPolynomial::new(remaining_vars_len, new_evaluations)
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
        }
    }
}


#[cfg(test)]
pub(crate) mod tests {
    use ark_ff::{AdditiveGroup, Field};

    use crate::constant_for_curves::ScalarField;

    use super::*;

    type F = ScalarField;

    // Helper function to create the polynomial used in the tests
    pub fn setup_polynomial() -> MultilinearPolynomial<F> {
        // Define the number of variables
        let num_variables = 3;

        // Define the evaluations over the boolean hypercube
        MultilinearPolynomial::new(num_variables, vec![
            F::from(10), // x1 = x2 = x3 = 0 -> 10
            F::from(11), // x1 = 1, x2 = x3 = 0 -> 11
            F::from(10), // x1 = 0, x2 = 1, x3 = 0 -> 10
            F::from(11), // x1 = x2 = 1, x3 = 0 -> 11
            F::from(10), // x1 = x2 = 0, x3 = 1 -> 10
            F::from(16), // x1 = 1, x2 = 0, x3 = 1 -> 16
            F::from(12), // x1 = 0, x2 = x3 = 1 -> 12
            F::from(18), // x1 = x2 = x3 = 1 -> 18
        ])
    }

    #[test]
    fn test_evaluation_over_boolean_hypercube() {
        let poly = setup_polynomial();

        // Test cases for the evaluate_binary function
        let test_cases = vec![
            (vec![F::ZERO, F::ZERO, F::ZERO], F::from(10)),
            (vec![F::ONE, F::ZERO, F::ZERO], F::from(11)),
            (vec![F::ONE, F::ZERO, F::ONE], F::from(16)),
            (vec![F::ZERO, F::ONE, F::ONE], F::from(12)),
        ];

        for (r, eval) in test_cases {
            assert_eq!(poly.evaluate_binary(r.as_slice()), eval);
            assert_eq!(poly.evaluate(r.as_slice()), eval);
        }

        // Test the evaluate function with a non-binary input
        let input = vec![F::from(3), F::from(4), F::from(5)]; // Example input
        let expected = F::from(128); // Expected result based on the polynomial f(x1, x2, x3) = 10 + x1 + 2 * x2x3 + 5 * x1x3
        assert_eq!(poly.evaluate(&input), expected);
    }

    #[test]
    fn test_partial_evaluation() {
        let poly = setup_polynomial();

        let p = poly.partial_evaluation(&[F::from(5)]);
        assert_eq!(
            p.evaluate(&[F::from(3), F::from(6)]),
            poly.evaluate(&[F::from(5), F::from(3), F::from(6)])
        );
    }
}

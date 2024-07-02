/// TODO: Code based on https://github.com/PolyhedraZK/Expander-rs/blob/main/bi-kzg/src/poly.rs
/// Changes:
/// - Adapted to arkworks instead of halo2curves
/// - Simplify it to only consider square matrices (degree_x == degree_y)
/// - Change the coefficient format of the polynomial to what is described below

use ark_ff::{Field, Zero, PrimeField};
use itertools::Itertools;
use rand::RngCore;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::univariate_poly::UnivariatePolynomial;
use crate::utils::compute_powers;

/// A bivariate polynomial in coefficient form:
///   f(X, Y) = sum_{i=0}^{d-1} sum_{j=0}^{d-1} f_{i,j} * X^i * Y^j
///
/// The coefficients of the bivariate polynomial are stored in a flat vector and logically
/// organized in a matrix of dimensions `d x d`. The matrix representation
/// allows us to associate the coefficients with the powers of x and y.
///
/// For example, given a polynomial with 16 coefficients and degree `d = 4`:
///
/// f(X, Y) = f_{0,0} + f_{0,1}Y + f_{0,2}Y^2 + f_{0,3}Y^3 +
///           f_{1,0}X + f_{1,1}XY + f_{1,2}XY^2 + f_{1,3}XY^3 +
///           f_{2,0}X^2 + f_{2,1}X^2Y + f_{2,2}X^2Y^2 + f_{2,3}X^2Y^3 +
///           f_{3,0}X^3 + f_{3,1}X^3Y + f_{3,2}X^3Y^2 + f_{3,3}X^3Y^3
///
/// The coefficients are stored in a vector as:
/// [f_{0,0}, f_{0,1}, f_{0,2}, f_{0,3}, f_{1,0}, f_{1,1}, f_{1,2}, f_{1,3},
///  f_{2,0}, f_{2,1}, f_{2,2}, f_{2,3}, f_{3,0}, f_{3,1}, f_{3,2}, f_{3,3}]
///
/// Each row in the matrix corresponds to increasing powers of Y, and within each row,
/// the columns correspond to increasing powers of X.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BivariatePolynomial<F> {
    pub coeffs: Vec<F>,  // Coefficients of the bivariate polynomial stored in a flat vector
    pub degree: usize,   // Degree of the polynomial in both X and Y
}

impl<F: PrimeField> BivariatePolynomial<F> {
    /// Creates a new bivariate polynomial with given coefficients and degree.
    #[inline]
    pub fn new(coefficients: Vec<F>, degree: usize) -> Self {
        assert_eq!(coefficients.len(), degree * degree); // Ensure the number of coefficients matches the square of the degree
        Self {
            coeffs: coefficients,
            degree,
        }
    }

    /// Test only: Generates a random bivariate polynomial with specified degree.
    pub fn random<T: RngCore>(rng: &mut T, degree: usize) -> Self {
        let coefficients = (0..degree * degree)
            .map(|_| F::rand(rng)) // Randomly generate coefficients
            .collect();
        Self::new(coefficients, degree)
    }

    /// Test only: Generates a random binary bivariate polynomial with specified degree.
    pub fn random_binary<T: RngCore>(_rng: &mut T, _degree: usize) -> Self {
        unimplemented!();
    }

    /// Evaluates the bivariate polynomial at given values (b, c).
    ///
    /// To evaluate the polynomial at (b, c), the function computes:
    /// f(a, b) = sum_{i=0} sum_{j=0} f_{i,j} * b^i * c^j
    pub fn evaluate(&self, b: &F, c: &F) -> F {
        let vec_b_powers = compute_powers(b, self.degree); // Compute powers of b: [b^0, b^1, ..., b^(d-1)]
        let vec_c_powers = compute_powers(c, self.degree); // Compute powers of c: [c^0, c^1, ..., c^(d-1)]

        // Iterate over the powers of b
        vec_b_powers.iter().enumerate().fold(F::ZERO, |acc, (i, x_i)| {
            let col_start = i * self.degree; // Start index of the current column in the flat vector
            let col = &self.coeffs[col_start..col_start + self.degree]; // Slice representing the current column

            // Accumulate the result by multiplying the coefficients with the corresponding powers of c and b
            acc + col.iter().zip(vec_c_powers.iter()).fold(F::ZERO, |acc, (c, y_i)| acc + *c * *y_i) * x_i
        })
    }

    /// Evaluates the polynomial at f(b,Y), returning a univariate polynomial in Y.
    ///
    /// To partially evaluate the polynomial at X = b, the function computes:
    /// f*(Y) = sum_{i=0} (sum_{j=0} f_{i,j} * b^i) * Y^j
    pub fn partially_evaluate_at_x(&self, b: &F) -> UnivariatePolynomial<F> {
        let mut f_y_a = vec![F::zero(); self.degree];
        let vec_b_powers = compute_powers(b, self.degree); // Compute powers of b: [b^0, b^1, ..., b^(d-1)]

        // Iterate over the powers of b
        for i in 0..self.degree {
            let col_start = i * self.degree; // Start index of the current column in the flat vector
            let col = &self.coeffs[col_start..col_start + self.degree]; // Slice representing the current column

            for (j, c) in col.iter().enumerate() {
                f_y_a[j] += *c * vec_b_powers[i]; // Add the result to the corresponding coefficient of the univariate polynomial in Y
            }
        }

        UnivariatePolynomial::new(f_y_a)
    }

    /// Evaluates the polynomial at f(X, y), returning a univariate polynomial in X.
    ///
    /// To partially evaluate the polynomial at y = b, the function computes:
    /// q(X) = sum_{i=0}^{d-1} (sum_{j=0}^{d-1} f_{i,j} * b^j) * X^i
    pub fn partially_evaluate_at_y(&self, y: &F) -> UnivariatePolynomial<F> {
        let mut f_x_b = vec![F::zero(); self.degree]; // Initialize with zero coefficients
        let powers_of_y = compute_powers(y, self.degree); // Compute powers of y: [y^0, y^1, ..., y^(d-1)]

        // Iterate over the chunks of coefficients corresponding to each power of x
        for i in 0..self.degree {
            for j in 0..self.degree {
                f_x_b[i] += self.coeffs[i * self.degree + j] * powers_of_y[j]; // Accumulate the coefficients multiplied by the corresponding power of y
            }
        }

        UnivariatePolynomial::new(f_x_b) // Return the resulting univariate polynomial in X
    }

    /// Compute r(x) = \sum_{j\inH_y} f(X, j)
    ///
    /// Evaluates the polynomial at all roots of unity in the domain and sums the results.
    pub fn sum_partial_evaluations_in_domain(&self) -> UnivariatePolynomial<F> {
        // XXX This could be precomputed and stored somewhere (e.g. on the poly itself)
        let domain = GeneralEvaluationDomain::<F>::new(self.degree).unwrap();

        // This can probably be sped up...
        let mut r_poly = UnivariatePolynomial::new(vec![F::zero(); self.degree]);
        for w_i in domain.elements() {
            r_poly = r_poly + self.partially_evaluate_at_y(&w_i);
        }

        r_poly
    }

    /// Computes the bitfield union of two bivariate polynomials.
    ///
    /// The coefficients of the resulting polynomial are the bitwise OR of the coefficients
    /// of the two input polynomials.
    pub fn bitfield_union(&self, other: &Self) -> Self {
        assert_eq!(self.degree, other.degree, "Polynomials must have the same degree");

        let coefficients: Vec<F> = self.coeffs.iter()
            .zip(&other.coeffs)
            .map(|(a, b)| {
                *a + *b - *a * *b // Since a, b are either 0 or 1, this is equivalent to a | b
            })
            .collect();

        Self::new(coefficients, self.degree)
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use ark_std::test_rng;
    use ark_std::UniformRand;
    use ark_bn254::{Fr, G1Projective};

    #[test]
    fn test_evaluate() {
        let poly = BivariatePolynomial::new(
            vec![
                Fr::from(1u64),  // f_00
                Fr::from(2u64),  // f_01
                Fr::from(3u64),  // f_02
                Fr::from(4u64),  // f_03
                Fr::from(5u64),  // f_10
                Fr::from(6u64),  // f_11
                Fr::from(7u64),  // f_12
                Fr::from(8u64),  // f_13
                Fr::from(9u64),  // f_20
                Fr::from(10u64), // f_21
                Fr::from(11u64), // f_22
                Fr::from(12u64), // f_23
                Fr::from(13u64), // f_30
                Fr::from(14u64), // f_31
                Fr::from(15u64), // f_32
                Fr::from(16u64), // f_33
            ],
            4, // degree
        );

        let x = Fr::from(2u64);
        let y = Fr::from(3u64);
        let result = poly.evaluate(&x, &y);

        // Manually compute the expected result:
        // f(2, 3) = 1 + 2*3 + 3*3^2 + 4*3^3 +
        //           5*2 + 6*2*3 + 7*2*3^2 + 8*2*3^3 +
        //           9*2^2 + 10*2^2*3 + 11*2^2*3^2 + 12*2^2*3^3 +
        //           13*2^3 + 14*2^3*3 + 15*2^3*3^2 + 16*2^3*3^3

        let expected_result = Fr::from(
            1 + 2*3 + 3*9 + 4*27 +
            5*2 + 6*6 + 7*18 + 8*54 +
            9*4 + 10*12 + 11*36 + 12*108 +
            13*8 + 14*24 + 15*72 + 16*216
        );

        assert_eq!(result, expected_result);
    }


    #[test]
    fn test_partial_eval_at_x() {
        let poly = BivariatePolynomial::new(
            vec![
                Fr::from(1u64),  // f_00
                Fr::from(2u64),  // f_01
                Fr::from(3u64),  // f_02
                Fr::from(4u64),  // f_03
                Fr::from(5u64),  // f_10
                Fr::from(6u64),  // f_11
                Fr::from(7u64),  // f_12
                Fr::from(8u64),  // f_13
                Fr::from(9u64),  // f_20
                Fr::from(10u64), // f_21
                Fr::from(11u64), // f_22
                Fr::from(12u64), // f_23
                Fr::from(13u64), // f_30
                Fr::from(14u64), // f_31
                Fr::from(15u64), // f_32
                Fr::from(16u64), // f_33
            ],
            4, // degree
        );
        let eval_at_x = poly.partially_evaluate_at_x(&Fr::from(2u64));

        // Manually compute the expected result:
        // For x = 2, we need to evaluate:
        // f(y) = (1 + 5*2 + 9*2^2 + 13*2^3) + (2 + 6*2 + 10*2^2 + 14*2^3)y +
        //        (3 + 7*2 + 11*2^2 + 15*2^3)y^2 + (4 + 8*2 + 12*2^2 + 16*2^3)y^3

        let expected_coeffs = vec![
            Fr::from(1 + 5*2 + 9*4 + 13*8), // constant term
            Fr::from(2 + 6*2 + 10*4 + 14*8), // coefficient of y
            Fr::from(3 + 7*2 + 11*4 + 15*8), // coefficient of y^2
            Fr::from(4 + 8*2 + 12*4 + 16*8), // coefficient of y^3
        ];

        assert_eq!(eval_at_x.poly.coeffs, expected_coeffs);
    }

    #[test]
    fn test_partial_eval_at_y() {
        let poly = BivariatePolynomial::new(
            vec![
                Fr::from(1u64),  Fr::from(2u64),  Fr::from(3u64),  Fr::from(4u64),
                Fr::from(5u64),  Fr::from(6u64),  Fr::from(7u64),  Fr::from(8u64),
                Fr::from(9u64),  Fr::from(10u64), Fr::from(11u64), Fr::from(12u64),
                Fr::from(13u64), Fr::from(14u64), Fr::from(15u64), Fr::from(16u64),
            ],
            4,
        );

        let y = Fr::from(2u64); // Evaluate at y = 2
        let eval_at_y = poly.partially_evaluate_at_y(&y);

        // Manually compute the expected result:
        // For y = 2, we need to evaluate:
        // f_x_b = [1 + 2*2 + 3*4 + 4*8, 5 + 6*2 + 7*4 + 8*8, 9 + 10*2 + 11*4 + 12*8, 13 + 14*2 + 15*4 + 16*8]
        let expected_coeffs = vec![
            Fr::from(1 + 2*2 + 3*4 + 4*8),
            Fr::from(5 + 6*2 + 7*4 + 8*8),
            Fr::from(9 + 10*2 + 11*4 + 12*8),
            Fr::from(13 + 14*2 + 15*4 + 16*8),
        ];

        assert_eq!(eval_at_y.poly.coeffs, expected_coeffs);
    }

    /// For polynomial f(X,Y) and polynomial f*(b,Y), check that f(b,c) == f*(c)
    #[test]
    fn test_partial_eval_end_to_end() {
        let rng = &mut test_rng();

        let b = &Fr::rand(rng);
        let c = &Fr::rand(rng);

        let bivariate_poly = BivariatePolynomial::random(rng, 16);

        // Compute f*(Y)
        let partial_eval_poly = bivariate_poly.partially_evaluate_at_x(b);
        // Compute f*(c)
        let y = partial_eval_poly.evaluate(c);

        // Compute f(b,c)
        let y_expected = bivariate_poly.evaluate(b, c);

        assert_eq!(y, y_expected);
    }

    #[test]
    fn test_bitfield_union() {
        // Create a prime field element type Fr with elements 0 and 1 representing false and true
        // Here, for illustration purposes, assume Fr::from(0u64) and Fr::from(1u64) represent 0 and 1
        let poly1 = BivariatePolynomial::new(
            vec![
                Fr::from(0u64), Fr::from(0u64), Fr::from(0u64), Fr::from(1u64), Fr::from(1u64),
                Fr::from(0u64), Fr::from(1u64), Fr::from(0u64), Fr::from(1u64), Fr::from(0u64),
                Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(0u64), Fr::from(0u64),
                Fr::from(0u64), Fr::from(0u64), Fr::from(1u64), Fr::from(1u64), Fr::from(0u64),
                Fr::from(1u64), Fr::from(0u64), Fr::from(0u64), Fr::from(1u64), Fr::from(1u64),
            ],
            5,
        );
        let poly2 = BivariatePolynomial::new(
            vec![
                Fr::from(0u64), Fr::from(1u64), Fr::from(1u64), Fr::from(0u64), Fr::from(1u64),
                Fr::from(1u64), Fr::from(0u64), Fr::from(1u64), Fr::from(0u64), Fr::from(1u64),
                Fr::from(1u64), Fr::from(0u64), Fr::from(0u64), Fr::from(1u64), Fr::from(1u64),
                Fr::from(1u64), Fr::from(1u64), Fr::from(0u64), Fr::from(1u64), Fr::from(1u64),
                Fr::from(0u64), Fr::from(1u64), Fr::from(1u64), Fr::from(0u64), Fr::from(1u64),
            ],
            5,
        );

        let union_poly = poly1.bitfield_union(&poly2);

        let expected_coeffs = vec![
            Fr::from(0u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64),
            Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64),
            Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64),
            Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64),
            Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64),
        ];

        assert_eq!(union_poly.coeffs, expected_coeffs);
    }
}

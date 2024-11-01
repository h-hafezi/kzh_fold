use ark_crypto_primitives::sponge::Absorb;
use ark_ff::{PrimeField, Field};
use ark_std::vec::Vec;


/// Struct for interpolating a polynomial from given evaluations at distinct points (x_i, f(x_i)).
#[derive(Debug)]
pub struct PolynomialInterpolator<F: PrimeField> {
    coefficients: Vec<F>,  // This stores the coefficients after interpolation.
}

impl<F: PrimeField> PolynomialInterpolator<F> {
    /// Constructs a new interpolator without coefficients.
    /// `interpolate` method will populate the coefficients.
    pub fn new() -> Self {
        Self {
            coefficients: Vec::new(),
        }
    }

    /// Interpolates a polynomial from given pairs of points (x_i, y_i).
    pub fn interpolate(&mut self, points: &[(F, F)]) {
        let n = points.len();
        self.coefficients = vec![F::zero(); n];

        // Loop over each Lagrange basis polynomial L_j(x)
        for (j, &(x_j, y_j)) in points.iter().enumerate() {
            let mut basis_coefficients = vec![F::one()];

            // Construct the Lagrange basis polynomial L_j(x)
            for (m, &(x_m, _)) in points.iter().enumerate() {
                if m != j {
                    let factor = x_j - x_m;
                    let inverse_factor = factor.inverse().unwrap(); // assuming `factor` is non-zero
                    let linear_term = vec![-x_m * inverse_factor, inverse_factor];

                    basis_coefficients = Self::poly_multiply(&basis_coefficients, &linear_term);
                }
            }

            // Multiply each coefficient of L_j(x) by y_j
            for (k, coeff) in basis_coefficients.iter().enumerate() {
                self.coefficients[k] += *coeff * y_j;
            }
        }
    }

    /// Multiplies two polynomials and returns the resulting coefficients.
    fn poly_multiply(a: &[F], b: &[F]) -> Vec<F> {
        let mut result = vec![F::zero(); a.len() + b.len() - 1];

        for (i, &coeff_a) in a.iter().enumerate() {
            for (j, &coeff_b) in b.iter().enumerate() {
                result[i + j] += coeff_a * coeff_b;
            }
        }

        result
    }

    /// Evaluates the interpolated polynomial at a given point `x`.
    pub fn evaluate(&self, x: F) -> F {
        let mut result = F::zero();
        let mut x_power = F::one();

        for &coeff in &self.coefficients {
            result += coeff * x_power;
            x_power *= x;
        }

        result
    }
}

#[cfg(test)]
mod test {
    use ark_ff::Field;
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::ScalarField;
    use crate::polynomial::univariate::PolynomialInterpolator;

    type F = ScalarField;

    /// Helper function to evaluate the polynomial P(x) = x^7 + 12 * x^6 + 5 * x^3 + 100 * x^2 + x + 9
    fn original_polynomial(x: F) -> F {
        x.pow([7]) +
            F::from(12u64) * x.pow([6]) +
            F::from(5u64) * x.pow([3]) +
            F::from(100u64) * x.pow([2]) +
            x +
            F::from(9u64)
    }

    #[test]
    fn test_polynomial_interpolation_with_random_points() {
        // Step 1: Generate random points (x_i, f(x_i))
        let mut points = Vec::new();
        for _ in 0..8 {
            let x = F::rand(&mut thread_rng());
            let y = original_polynomial(x);
            points.push((x, y));
        }

        // Step 2: Interpolate the polynomial from the points
        let mut interpolator = PolynomialInterpolator::new();
        interpolator.interpolate(&points);

        // Step 3: Verify the interpolated polynomial matches the original polynomial
        // at a few random points
        for _ in 0..5 {
            let x = F::rand(&mut thread_rng());
            let original_value = original_polynomial(x);
            let interpolated_value = interpolator.evaluate(x);
            assert_eq!(original_value, interpolated_value, "Mismatch at random x = {:?}", x);
        }
    }
}

use ark_ff::{Field, PrimeField};

use crate::polynomial::multilinear_polynomial::math::Math;

pub mod multilinear_poly;

pub mod eq_poly;

pub mod math;
mod indentity;

pub fn compute_dot_product<F: PrimeField>(a: &[F], b: &[F]) -> F {
    assert_eq!(a.len(), b.len());
    (0..a.len()).map(|i| a[i] * b[i]).sum()
}

pub fn boolean_vector_to_decimal<F: Field>(r: &[F]) -> usize {
    // Ensure that every entry in the vector is either F::ONE or F::ZERO
    for &x in r {
        assert!(x == F::ONE || x == F::ZERO, "Vector contains non-boolean values.");
    }

    // Convert the boolean vector into a decimal index (big-endian)
    let mut decimal = 0;
    for (i, &bit) in r.iter().rev().enumerate() {
        if bit == F::ONE {
            decimal += 1 << (r.len() - 1 - i);
        }
    }

    decimal
}

pub fn decimal_to_boolean_vector<F: Field>(d: usize, length: usize) -> Vec<F> {
    let bits = d.get_bits(length);
    let mut output = Vec::with_capacity(length);
    for b in bits {
        output.push(if b { F::ONE } else { F::ZERO });
    }
    output
}

#[cfg(test)]
mod tests {
    use ark_ff::{AdditiveGroup, Field};

    use crate::constant_for_curves::ScalarField;
    use crate::polynomial::multilinear_polynomial::{boolean_vector_to_decimal, decimal_to_boolean_vector};

    type F = ScalarField;

    #[test]
    #[should_panic(expected = "Vector contains non-boolean values.")]
    fn test_boolean_vector_to_decimal() {
        assert_eq!(boolean_vector_to_decimal(&vec![F::ONE, F::ZERO, F::ZERO]), 1);
        assert_eq!(boolean_vector_to_decimal(&vec![F::ZERO, F::ZERO, F::ZERO]), 0);
        assert_eq!(boolean_vector_to_decimal(&vec![F::ONE, F::ONE, F::ONE]), 7);
        assert_eq!(boolean_vector_to_decimal(&vec![F::ZERO, F::ONE, F::ZERO]), 2);
        assert_eq!(decimal_to_boolean_vector::<F>(7, 3), vec![F::ONE, F::ONE, F::ONE]);
        assert_eq!(decimal_to_boolean_vector::<F>(0, 3), vec![F::ZERO, F::ZERO, F::ZERO]);
        assert_eq!(decimal_to_boolean_vector::<F>(1, 3), vec![F::ONE, F::ZERO, F::ZERO]);
        // Should panic
        boolean_vector_to_decimal(&vec![F::ONE, F::ONE + F::ONE, F::ZERO]);
    }
}


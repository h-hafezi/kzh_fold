use ark_ff::{Field, PrimeField};

pub mod dense_multilinear_poly;

pub mod eq_poly;

pub mod math;

fn compute_dot_product<F: PrimeField>(a: &[F], b: &[F]) -> F {
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

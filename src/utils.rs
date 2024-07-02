use ark_ff::{Field};
use rayon::prelude::*;

/// return [x^0, x^1, ..., x^n-1]
pub(crate) fn compute_powers<F: Field + Send + Sync>(x: &F, n: usize) -> Vec<F> {
    let sqrt_n = (n as f64).sqrt().ceil() as usize;

    let mut initial_powers = vec![F::ONE];
    let mut cur = *x;
    for _ in 1..sqrt_n {
        initial_powers.push(cur);
        cur *= x;
    }

    let powers: Vec<F> = (0..sqrt_n)
        .into_par_iter()
        .flat_map(|i| {
            let mut block = Vec::with_capacity(sqrt_n);
            if i == 0 {
                block.extend_from_slice(&initial_powers);
            } else {
                let initial_power = initial_powers[sqrt_n - 1] * x;
                let mut base = initial_powers[sqrt_n - 1] * x;
                for _ in 1..i {
                    base *= initial_power;
                }
                let mut cur = base;
                block.push(cur);
                for _ in 1..sqrt_n {
                    cur *= x;
                    block.push(cur);
                }
            }
            block
        })
        .collect();

    powers.into_iter().take(n).collect()
}

pub(crate) fn inner_product<F: Field>(vector1: &[F], vector2: &[F]) -> F {
    // Check if the lengths of the vectors are the same
    assert_eq!(vector1.len(), vector2.len(), "The two vectors must have the same size.");
    // Compute the inner product
    vector1.iter().zip(vector2.iter()).map(|(a, b)| (*a) * (*b)).sum()
}

#[cfg(test)]
mod tests {
    use ark_ff::{AdditiveGroup, Field};
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::ScalarField;
    use super::*;

    pub(crate) fn compute_powers_non_parallel<F: Field>(x: &F, n: usize) -> Vec<F> {
        let mut powers = vec![F::ONE];
        let mut cur = *x;
        for _ in 0..n - 1 {
            powers.push(cur);
            cur *= x;
        }
        powers
    }

    type F = ScalarField;

    #[test]
    fn test_compute_powers() {
        let f = ScalarField::rand(&mut thread_rng());
        let n = 16;

        let result_original = compute_powers_non_parallel(&f, n);
        let result_parallel = compute_powers(&f, n);

        assert_eq!(result_original, result_parallel, "The results of the original and parallel implementations do not match.");

        let n = 32;

        let result_original = compute_powers_non_parallel(&f, n);
        let result_parallel = compute_powers(&f, n);

        assert_eq!(result_original, result_parallel, "The results of the original and parallel implementations do not match.");


        let n = 1234567;

        let result_original = compute_powers_non_parallel(&f, n);
        let result_parallel = compute_powers(&f, n);

        assert_eq!(result_original, result_parallel, "The results of the original and parallel implementations do not match.");


        let n = 1001;

        let result_original = compute_powers_non_parallel(&f, n);
        let result_parallel = compute_powers(&f, n);

        assert_eq!(result_original, result_parallel, "The results of the original and parallel implementations do not match.");
    }

    #[test]
    fn inner_product_test() {
        let vec1 = vec![F::ONE, F::ONE, F::ONE];
        let vec2 = vec![F::ZERO, F::ONE, F::ZERO];
        assert_eq!(inner_product(vec1.as_slice(), vec2.as_slice()), F::ONE);
        let vec1 = vec![F::ONE, F::ONE, F::ONE];
        let vec2 = vec![F::ONE, F::ONE, F::ZERO];
        assert_eq!(inner_product(vec1.as_slice(), vec2.as_slice()), F::from(2u8));
        let vec1 = vec![F::ONE, F::ONE, F::ONE];
        let vec2 = vec![F::ONE, F::ONE, F::ONE];
        assert_eq!(inner_product(vec1.as_slice(), vec2.as_slice()), F::from(3u8));
    }
}


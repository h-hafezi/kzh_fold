use ark_ff::{Field, PrimeField};
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

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use super::*;
    use ark_ff::{Field, PrimeField};
    use ark_std::UniformRand;
    use num::BigInt;
    use rand::{Rng, thread_rng};

    pub(crate) fn compute_powers_non_parallel<F: Field>(x: &F, n: usize) -> Vec<F> {
        let mut powers = vec![F::ONE];
        let mut cur = *x;
        for _ in 0..n - 1 {
            powers.push(cur);
            cur *= x;
        }
        powers
    }

    #[test]
    fn test_compute_powers() {
        let f = Fr::rand(&mut thread_rng());
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
}


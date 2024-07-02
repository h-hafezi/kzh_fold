use ark_ff::{Field, PrimeField};

/// return [x^0, x^1, ..., x^n-1]
pub(crate) fn compute_powers<F: Field>(x: &F, n: usize) -> Vec<F> {
    let mut powers = vec![F::ONE];
    let mut cur = *x;
    for _ in 0..n - 1 {
        powers.push(cur);
        cur *= x;
    }
    powers
}


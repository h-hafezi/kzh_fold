#![allow(warnings)]

use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ff::{AdditiveGroup, Field, PrimeField};
use ark_std::UniformRand;
use rand::Rng;

use crate::math::Math;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;

pub mod kzh2_fold;
mod eq_tree;
pub mod kzh_4_fold;

/// returns a vector of length "degree" of E::G1Affine random elements
fn generate_random_elements<E: Pairing, R: Rng>(degree: usize, rng: &mut R) -> Vec<E::G1Affine> {
    let mut elements = Vec::new();
    for _ in 0..degree {
        elements.push(E::G1Affine::rand(rng));
    }
    elements
}

use rayon::prelude::*; // Ensure you import rayon for parallel iterators

pub fn generic_linear_combination<T, FN>(vec1: &[T], vec2: &[T], combine_fn: FN) -> Vec<T>
where
    T: Send + Sync + Copy, // Traits required for parallel processing and copying
    FN: Fn(T, T) -> T + Sync, // Closure must work on elements of type T
{
    // Ensure both vectors have the same length
    assert_eq!(vec1.len(), vec2.len(), "Vectors must have the same length.");

    // Apply the closure to each pair of entries in the vectors
    vec1.par_iter()
        .zip(vec2.par_iter())
        .map(|(&a, &b)| combine_fn(a, b))
        .collect()
}

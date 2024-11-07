#![allow(warnings)]

use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ff::{AdditiveGroup, Field, PrimeField};
use ark_std::UniformRand;
use rand::Rng;

use crate::math::Math;

pub mod accumulator;
mod eq_tree;

/// returns a vector of length "degree" of E::G1Affine random elements
fn generate_random_elements<E: Pairing, R: Rng>(degree: usize, rng: &mut R) -> Vec<E::G1Affine> {
    let mut elements = Vec::new();
    for _ in 0..degree {
        elements.push(E::G1Affine::rand(rng));
    }
    elements
}

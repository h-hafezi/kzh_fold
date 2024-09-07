#![allow(warnings)]

use ark_ec::AffineRepr;
use ark_ec::pairing::Pairing;
use ark_ff::{AdditiveGroup, Field, PrimeField};
use ark_std::UniformRand;
use rand::Rng;

use crate::gadgets::non_native::util::convert_field_one_to_field_two;
use crate::polynomial::multilinear_polynomial::math::Math;

pub mod accumulator;

fn convert_affine_to_scalars<E: Pairing>(point: E::G1Affine) -> (E::ScalarField, E::ScalarField)
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    if point.is_zero() {
        (E::ScalarField::ONE, E::ScalarField::ZERO)
    } else {
        // Extract x and y coordinates and convert them
        let x = convert_field_one_to_field_two::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(point.x().unwrap());
        let y = convert_field_one_to_field_two::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(point.y().unwrap());
        (x, y)
    }
}

fn generate_random_elements<E: Pairing, R: Rng>(degree: usize, rng: &mut R) -> Vec<E::G1Affine> {
    let mut elements = Vec::new();
    for _ in 0..degree {
        elements.push(E::G1Affine::rand(rng));
    }
    elements
}

use std::marker::PhantomData;

use ark_ec::pairing::Pairing;
use ark_ff::{PrimeField, UniformRand};
use ark_std::test_rng;
use crate::transcript::transcript::Transcript;

pub struct RandomTape<E: Pairing> {
    pub tape: Transcript<E::ScalarField>,
}

impl<E: Pairing> RandomTape<E> {
    pub fn new(name: &'static [u8]) -> Self {
        let tape = {
            let mut prng = test_rng();
            let mut tape = Transcript::new(name);
            Transcript::append_scalar(
                &mut tape,
                b"init_randomness",
                &E::ScalarField::rand(&mut prng),
            );
            tape
        };
        Self {
            tape,
        }
    }

    pub fn random_scalar(&mut self, label: &'static [u8]) -> E::ScalarField {
        Transcript::challenge_scalar(&mut self.tape, label)
    }

    pub fn random_vector(&mut self, label: &'static [u8], len: usize) -> Vec<E::ScalarField> {
        Transcript::challenge_vector(&mut self.tape, label, len)
    }
}
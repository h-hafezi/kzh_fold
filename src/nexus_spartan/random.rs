use std::marker::PhantomData;

use super::transcript::ProofTranscript;
use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use ark_std::test_rng;
use merlin::Transcript;

pub struct RandomTape<E> {
    pub(crate) tape: Transcript,
    phantom: PhantomData<E>,
}

impl<E: Pairing> RandomTape<E> {
    pub fn new(name: &'static [u8]) -> Self {
        let tape = {
            let mut prng = test_rng();
            let mut tape = Transcript::new(name);
            <Transcript as ProofTranscript<E>>::append_scalar(
                &mut tape,
                b"init_randomness",
                &E::ScalarField::rand(&mut prng),
            );
            tape
        };
        Self {
            tape,
            phantom: PhantomData,
        }
    }

    pub fn random_scalar(&mut self, label: &'static [u8]) -> E::ScalarField {
        <Transcript as ProofTranscript<E>>::challenge_scalar(&mut self.tape, label)
    }

    pub fn random_vector(&mut self, label: &'static [u8], len: usize) -> Vec<E::ScalarField> {
        <Transcript as ProofTranscript<E>>::challenge_vector(&mut self.tape, label, len)
    }
}
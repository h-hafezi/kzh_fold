use ark_bn254::g1::G1Affine;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;

pub struct Transcript<F: PrimeField> {
    // This will hold the current state of the transcript
    state: F,
}

impl<F: PrimeField> Transcript<F> {
    pub fn new(label: &'static [u8]) -> Transcript<F> {
        Transcript {
            state: F::ONE,
        }
    }
}

impl<F: PrimeField> Transcript<F> {
    pub fn append_u64(&mut self, label: &'static [u8], n: u64) {}

    pub fn append_message(&mut self, label: &'static [u8], msg: &[u8]) {}

    pub fn append_scalar(&mut self, label: &'static [u8], scalar: &F) {}

    pub fn append_scalars(&mut self, label: &'static [u8], scalars: &[F]) {}

    pub fn challenge_scalar(&mut self, label: &'static [u8]) -> F {
        F::ONE
    }

    pub fn challenge_vector(&mut self, label: &'static [u8], len: usize) -> Vec<F> {
        vec![F::ONE; len]
    }

    pub(crate) fn append_protocol_name(&mut self, protocol_name: &'static [u8]) {}
}

impl<F: PrimeField> Transcript<F> {
    pub fn append_point<E: Pairing<ScalarField=F>>(&mut self, label: &'static [u8], point: &E::G1Affine) {}

    pub fn append_points<E: Pairing<ScalarField=F>>(&mut self, label: &'static [u8], point: &[E::G1Affine]) {}
}

pub trait AppendToTranscript<F: PrimeField> {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>);
}

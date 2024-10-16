use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::PrimeField;

pub struct Transcript<F: PrimeField> {
    // This will hold the current state of the transcript
    state: F,
}

impl<F: PrimeField> Transcript<F> {
    pub fn append_scalar(&mut self, label: &'static [u8], scalar: &F) {

    }

    pub fn append_scalars(&mut self, label: &'static [u8], scalars: &[F]) {

    }

    pub fn challenge_scalar(&mut self, label: &'static [u8]) -> F {
        F::ONE
    }

    pub fn challenge_vector(&mut self, label: &'static [u8], len: usize) -> Vec<F> {
        vec![F::ONE; len]
    }

    fn append_protocol_name(&mut self, protocol_name: &'static [u8]) {

    }
}

impl<F: PrimeField> Transcript<F> {
    pub fn append_point<G: SWCurveConfig<ScalarField=F>>(&mut self, label: &'static [u8], point: &Projective<G>) {

    }

    pub fn append_points<G: SWCurveConfig<ScalarField=F>>(&mut self, label: &'static [u8], point: &Projective<G>) {

    }
}
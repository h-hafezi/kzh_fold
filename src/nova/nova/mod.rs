use ark_ec::AffineRepr;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use crate::constant_for_curves::BaseField;

pub mod verifier_circuit;
pub mod prover;
mod verifier_circuit_var;

#[inline]
fn get_affine_coords<F: PrimeField, G: SWCurveConfig<BaseField=F>>(affine: &Affine<G>) -> (F, F) {
    if affine.is_zero() {
        (F::zero(), F::zero())
    } else {
        let (x, y) = affine.xy().unwrap();
        (x, y)
    }
}

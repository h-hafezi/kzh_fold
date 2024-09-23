#![allow(warnings)]
use ark_ec::AffineRepr;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::{One, Zero};
pub mod instance_circuit;

pub mod verifier_circuit;
pub mod prover;

pub fn affine_to_projective<P: SWCurveConfig>(a: Affine<P>) -> Projective<P> {
    if a.is_zero() {
        // If e.C is the point at infinity (zero in projective coordinates), return the projective zero element
        Projective::zero()
    } else {
        // Otherwise, convert the affine point to a projective point
        Projective::new(
            a.x().unwrap(),
            a.y().unwrap(),
            P::BaseField::one(), // z coordinate set to 1 for standard affine-to-projective conversion
        )
    }
}

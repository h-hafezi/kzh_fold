#![allow(warnings)]
use ark_ec::AffineRepr;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::{One, Zero};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_relations::ns;
use ark_relations::r1cs::ConstraintSystemRef;
use crate::gadgets::non_native::util::cast_field;

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

pub fn randomness_different_formats<E: Pairing>(cs: ConstraintSystemRef<E::ScalarField>, beta: E::ScalarField) -> (
    E::BaseField,
    FpVar<E::ScalarField>,
    NonNativeFieldVar<E::BaseField, E::ScalarField>
) {
    let beta_base = cast_field::<E::ScalarField, E::BaseField>(beta);
    let beta_var = FpVar::new_variable(
        ns!(cs, "beta var"),
        || Ok(beta.clone()),
        AllocationMode::Witness,
    ).unwrap();
    let beta_var_non_native = NonNativeFieldVar::new_variable(
        ns!(cs, "beta var non-native"),
        || Ok(beta_base.clone()),
        AllocationMode::Witness,
    ).unwrap();
    (beta_base, beta_var, beta_var_non_native)
}
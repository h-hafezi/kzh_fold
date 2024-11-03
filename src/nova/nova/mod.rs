use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::AffineVar;
use ark_r1cs_std::select::CondSelectGadget;

pub mod verifier_circuit;
pub mod prover;
mod verifier_circuit_var;


// the following functions have their non-native implementation in gadgets/non_native/non_native_affine_var.rs
// but here we provide the corresponding ones for the native AffineVar

fn get_affine_coords<F, G>(affine: &Affine<G>) -> (F, F)
where
    F: PrimeField,
    G: SWCurveConfig<BaseField=F>,
{
    if affine.is_zero() {
        (F::one(), F::zero())
    } else {
        let (x, y) = affine.xy().unwrap();
        (x, y)
    }
}

fn get_affine_var_coords<F, G>(affine: &AffineVar<G, FpVar<F>>) -> Vec<FpVar<F>>
where
    F: PrimeField,
    G: SWCurveConfig<BaseField=F>,
{
    // Define constants for (1, 0) to represent the coordinates when the point is at infinity
    let one_fpvar = FpVar::constant(F::one());
    let zero_fpvar = FpVar::constant(F::zero());

    // Extract x and y coordinates safely (these will only be used if the point is not at infinity)

    // Conditionally select the output based on whether the point is at infinity
    let x = FpVar::conditionally_select(&affine.infinity, &one_fpvar, &affine.x.clone()).unwrap();
    let y = FpVar::conditionally_select(&affine.infinity, &zero_fpvar, &affine.y.clone()).unwrap();

    vec![x, y]
}

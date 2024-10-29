use std::borrow::Borrow;

use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
use ark_ec::{
    CurveGroup,
    short_weierstrass::{Projective, SWCurveConfig},
};
use ark_ff::{AdditiveGroup, Field, PrimeField};
use ark_r1cs_std::{
    alloc::{AllocationMode, AllocVar},
    boolean::Boolean,
    eq::EqGadget,
    fields::{FieldVar, fp::FpVar, nonnative::NonNativeFieldVar},
    R1CSVar,
    select::CondSelectGadget,
    uint8::UInt8,
};
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use ark_std::Zero;

use crate::gadgets::non_native::util::non_native_to_fpvar;

#[must_use]
#[derive(Debug)]
pub struct NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    pub x: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,
    pub y: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,
    pub infinity: Boolean<G1::ScalarField>,
}

impl<G1> Clone for NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
            infinity: self.infinity.clone(),
        }
    }
}

impl<G1> NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    pub fn into_projective(
        &self,
    ) -> Result<Vec<NonNativeFieldVar<G1::BaseField, G1::ScalarField>>, SynthesisError> {
        let zero_x = NonNativeFieldVar::zero();
        let zero_y = NonNativeFieldVar::one();

        let x = self.infinity.select(&zero_x, &self.x)?;
        let y = self.infinity.select(&zero_y, &self.y)?;
        let z = NonNativeFieldVar::from(self.infinity.not());

        Ok(vec![x, y, z])
    }
}

impl<G1> R1CSVar<G1::ScalarField> for NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    type Value = Projective<G1>;

    fn cs(&self) -> ConstraintSystemRef<G1::ScalarField> {
        self.x.cs().or(self.y.cs()).or(self.infinity.cs())
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        Ok(if self.infinity.value()? {
            Projective::zero()
        } else {
            Projective {
                x: self.x.value()?,
                y: self.y.value()?,
                z: G1::BaseField::ONE,
            }
        })
    }
}

impl<G1> AllocVar<Projective<G1>, G1::ScalarField> for NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    fn new_variable<T: Borrow<Projective<G1>>>(
        cs: impl Into<Namespace<G1::ScalarField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let g = match f() {
            Ok(g) => *g.borrow(),
            Err(_) => Projective::zero(),
        };

        let affine = g.into_affine();

        let x = NonNativeFieldVar::new_variable(
            cs.clone(),
            || Ok(affine.x),
            mode,
        ).unwrap();

        let y = NonNativeFieldVar::new_variable(
            cs.clone(),
            || Ok(affine.y),
            mode,
        ).unwrap();

        let infinity = Boolean::new_variable(
            cs.clone(),
            || Ok(affine.infinity),
            mode
        ).unwrap();

        Ok(Self { x, y, infinity })
    }
}

impl<G1> AbsorbGadget<G1::ScalarField> for NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    fn to_sponge_bytes(&self) -> Result<Vec<UInt8<G1::ScalarField>>, SynthesisError> {
        unreachable!()
    }

    fn to_sponge_field_elements(&self) -> Result<Vec<FpVar<G1::ScalarField>>, SynthesisError> {
        // Convert the non-native fields to FpVars
        let x_fpvar = non_native_to_fpvar(&self.x);
        let y_fpvar = non_native_to_fpvar(&self.y);

        // Define the constants for infinity (FpVar::ONE, FpVar::ZERO)
        let one_fpvar = FpVar::constant(G1::ScalarField::ONE);
        let zero_fpvar = FpVar::constant(G1::ScalarField::ZERO);

        // Conditionally select based on the infinity flag
        let x = FpVar::conditionally_select(&self.infinity, &one_fpvar, &x_fpvar)?;
        let y = FpVar::conditionally_select(&self.infinity, &zero_fpvar, &y_fpvar)?;

        Ok(vec![x, y])
    }
}

impl<G1> CondSelectGadget<G1::ScalarField> for NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    fn conditionally_select(
        cond: &Boolean<G1::ScalarField>,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError> {
        let x = cond.select(&true_value.x, &false_value.x)?;
        let y = cond.select(&true_value.y, &false_value.y)?;
        let infinity = cond.select(&true_value.infinity, &false_value.infinity)?;

        Ok(Self { x, y, infinity })
    }
}

impl<G1> NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    pub fn enforce_equal(
        &self,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        let projective_self = self.into_projective().unwrap();
        let projective_other = other.into_projective().unwrap();
        for i in 0..3 {
            projective_self[i].enforce_equal(&projective_other[i]).unwrap();
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
    use ark_ec::short_weierstrass::Projective;
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_relations::ns;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{G1, ScalarField};
    use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;

    #[test]
    fn constraint_count_test() {
        let cs: ConstraintSystemRef<ScalarField> = ConstraintSystem::new_ref();
        let Q_var: NonNativeAffineVar<G1> = NonNativeAffineVar::new_variable(
            ns!(cs, "Q"),
            || Ok(Projective::rand(&mut thread_rng())),
            AllocationMode::Witness,
        ).unwrap();
        println!("{}", cs.num_constraints());
        let _ = Q_var.to_sponge_field_elements().unwrap();
        println!("constraint count for NonNativeAffineVar into native field: {}", cs.num_constraints());
    }
}


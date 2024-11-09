// borrowed from Nexus: https://github.com/nexus-xyz/nexus-zkvm/blob/f37401c477b680ce5334b2ca523ded8a7273d8c8/nova/src/gadgets/cyclefold/nova/primary.rs
// these structs are used to implement Nova, e.g. running/new instance on the primary curve

use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::gadgets::r1cs::{R1CSInstance, RelaxedR1CSInstance};
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::{Field, PrimeField};
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct R1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    /// Commitment to witness.
    pub commitment_W: NonNativeAffineVar<G1>,
    /// Public input of non-relaxed instance.
    pub X: Vec<FpVar<G1::ScalarField>>,

    _commitment_scheme: PhantomData<C1>,
}


impl<G1, C1> Clone for R1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    fn clone(&self) -> Self {
        Self {
            commitment_W: self.commitment_W.clone(),
            X: self.X.clone(),
            _commitment_scheme: self._commitment_scheme,
        }
    }
}

impl<G1, C1> R1CSVar<G1::ScalarField> for R1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>>,
{
    type Value = R1CSInstance<G1, C1>;

    fn cs(&self) -> ConstraintSystemRef<G1::ScalarField> {
        self.X
            .iter()
            .fold(ConstraintSystemRef::None, |cs, x| cs.or(x.cs()))
            .or(self.commitment_W.cs())
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let commitment_W = self.commitment_W.value()?;
        let X = self.X.value()?;
        Ok(R1CSInstance { commitment_W: commitment_W.into(), X })
    }
}

impl<G1, C1> AllocVar<R1CSInstance<G1, C1>, G1::ScalarField> for R1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>>,
{
    fn new_variable<T: Borrow<R1CSInstance<G1, C1>>>(
        cs: impl Into<Namespace<G1::ScalarField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let r1cs = f()?;
        let X = &r1cs.borrow().X;
        // Only allocate valid instance, which starts with F::ONE.
        assert_eq!(X[0], G1::ScalarField::ONE);

        let commitment_W = NonNativeAffineVar::new_variable(
            cs.clone(),
            || Ok(r1cs.borrow().commitment_W.into()),
            mode,
        )?;
        let alloc_X = X[1..]
            .iter()
            .map(|x| FpVar::<G1::ScalarField>::new_variable(cs.clone(), || Ok(x), mode));
        let X = std::iter::once(Ok(FpVar::constant(G1::ScalarField::ONE)))
            .chain(alloc_X)
            .collect::<Result<_, _>>()?;

        Ok(Self {
            commitment_W,
            X,
            _commitment_scheme: PhantomData,
        })
    }
}

#[derive(Debug)]
pub struct RelaxedR1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    /// Commitment to witness.
    pub commitment_W: NonNativeAffineVar<G1>,
    /// Commitment to error vector.
    pub commitment_E: NonNativeAffineVar<G1>,
    /// Public input of relaxed instance. Expected to start with `u`.
    pub X: Vec<FpVar<G1::ScalarField>>,

    _commitment_scheme: PhantomData<C1>,
}

impl<G1, C1> Clone for RelaxedR1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    fn clone(&self) -> Self {
        Self {
            commitment_W: self.commitment_W.clone(),
            commitment_E: self.commitment_E.clone(),
            X: self.X.clone(),
            _commitment_scheme: self._commitment_scheme,
        }
    }
}

impl<G1, C1> R1CSVar<G1::ScalarField> for RelaxedR1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>>,
{
    type Value = RelaxedR1CSInstance<G1, C1>;

    fn cs(&self) -> ConstraintSystemRef<G1::ScalarField> {
        self.X
            .iter()
            .fold(ConstraintSystemRef::None, |cs, x| cs.or(x.cs()))
            .or(self.commitment_W.cs())
            .or(self.commitment_E.cs())
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let commitment_W = self.commitment_W.value()?;
        let commitment_E = self.commitment_E.value()?;

        let X = self.X.value()?;
        Ok(RelaxedR1CSInstance {
            commitment_W: commitment_W.into(),
            commitment_E: commitment_E.into(),
            X,
        })
    }
}

impl<G1, C1> AllocVar<RelaxedR1CSInstance<G1, C1>, G1::ScalarField>
for RelaxedR1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>>,
{
    fn new_variable<T: Borrow<RelaxedR1CSInstance<G1, C1>>>(
        cs: impl Into<Namespace<G1::ScalarField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let r1cs = f()?;
        let X = &r1cs.borrow().X;

        let commitment_W = NonNativeAffineVar::new_variable(
            cs.clone(),
            || Ok(r1cs.borrow().commitment_W.into()),
            mode,
        )?;
        let commitment_E = NonNativeAffineVar::new_variable(
            cs.clone(),
            || Ok(r1cs.borrow().commitment_E.into()),
            mode,
        )?;

        let X = X
            .iter()
            .map(|x| FpVar::<G1::ScalarField>::new_variable(cs.clone(), || Ok(x), mode))
            .collect::<Result<_, _>>()?;

        Ok(Self {
            commitment_W,
            commitment_E,
            X,
            _commitment_scheme: PhantomData,
        })
    }
}

#[cfg(test)]
mod test {
    use crate::commitment::CommitmentScheme;
    use crate::constant_for_curves::{ScalarField, C1, G1};
    use crate::gadgets::r1cs::conversion::{get_random_r1cs_instance_witness, get_random_relaxed_r1cs_instance_witness};
    use ark_ec::short_weierstrass::Affine;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use crate::gadgets::r1cs::r1cs_var::{R1CSInstanceVar, RelaxedR1CSInstanceVar};

    type F = ScalarField;


    #[test]
    fn test() {
        // simply write a test that allocates R1CSVar and RelaxedR1CSVar and then using .value().unwrap() checks out the result
        let (num_constraints, num_io, num_vars) = (10, 3, 17);
        let pp: Vec<Affine<G1>> = C1::setup(num_vars, b"test", &());
        let cs = ConstraintSystem::new_ref();

        let (_, instance, _) = get_random_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &pp);

        // Allocate value
        let instance_var = R1CSInstanceVar::new_variable(
              cs.clone(),
              || Ok(instance.clone()),
              AllocationMode::Witness,
        ).unwrap();

        // assert equality when converted back to non-zk
        assert_eq!(instance, instance_var.value().unwrap());

        // generate a relaxed instance/witness this time
        let (_, relaxed_instance, _) = get_random_relaxed_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &pp);

        // Allocate value
        let relaxed_instance_var = RelaxedR1CSInstanceVar::new_variable(
            cs.clone(),
            || Ok(relaxed_instance.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // assert equality when converted back to non-zk
        assert_eq!(relaxed_instance, relaxed_instance_var.value().unwrap());
    }
}

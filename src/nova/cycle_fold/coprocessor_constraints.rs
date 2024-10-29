use std::{borrow::Borrow, marker::PhantomData};

use ark_ec::{
    short_weierstrass::{Affine, Projective, SWCurveConfig},
    AffineRepr,
};
use ark_ff::{Field, PrimeField};
use ark_r1cs_std::{
    alloc::{AllocVar, AllocationMode},
    boolean::Boolean,
    fields::{
        fp::FpVar,
        nonnative::{NonNativeFieldMulResultVar, NonNativeFieldVar},
        FieldVar,
    },
    groups::{curves::short_weierstrass::ProjectiveVar, CurveVar},
    select::CondSelectGadget,
    R1CSVar,
    ToBitsGadget,
};
use ark_r1cs_std::eq::EqGadget;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};

use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::gadgets::r1cs::{OvaInstance, RelaxedOvaInstance};
use crate::nova::cycle_fold::coprocessor::SecondaryCircuit as SecondaryCircuit;

#[must_use]
#[derive(Debug)]
/// Struct native on G2::BaseField = G1::ScalarField, hence native on primary curve
pub struct OvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
{
    /// Commitment to witness.
    pub commitment: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    /// Public input of non-relaxed instance.
    pub X: Vec<NonNativeFieldVar<G2::ScalarField, G2::BaseField>>,

    _commitment_scheme: PhantomData<C2>,
}

impl<G2, C2> Clone for OvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
{
    fn clone(&self) -> Self {
        Self {
            commitment: self.commitment.clone(),
            X: self.X.clone(),
            _commitment_scheme: self._commitment_scheme,
        }
    }
}

impl<G2, C2> R1CSVar<G2::BaseField> for OvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    type Value = OvaInstance<G2, C2>;

    fn cs(&self) -> ConstraintSystemRef<G2::BaseField> {
        self.X
            .iter()
            .fold(ConstraintSystemRef::None, |cs, x| cs.or(x.cs()))
            .or(self.commitment.cs())
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let commitment = self.commitment.value()?.into();
        let X = self.X.value()?;
        Ok(OvaInstance { commitment, X })
    }
}

impl<G2, C2> AllocVar<OvaInstance<G2, C2>, G2::BaseField> for OvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    fn new_variable<T: Borrow<OvaInstance<G2, C2>>>(
        cs: impl Into<Namespace<G2::BaseField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let ova_instance = f()?;
        let X = &ova_instance.borrow().X;
        // Only allocate valid instance, which starts with F::ONE.
        assert_eq!(X[0], G2::ScalarField::ONE);

        let commitment = ProjectiveVar::<G2, FpVar<G2::BaseField>>::new_variable(
            cs.clone(),
            || Ok(ova_instance.borrow().commitment.into()),
            mode,
        )?;
        let alloc_X = X[1..]
            .iter()
            .map(|x| NonNativeFieldVar::new_variable(cs.clone(), || Ok(x), mode));

        let X = std::iter::once(Ok(NonNativeFieldVar::constant(G2::ScalarField::ONE)))
            .chain(alloc_X)
            .collect::<Result<_, _>>()?;

        Ok(Self {
            commitment,
            X,
            _commitment_scheme: PhantomData,
        })
    }
}

impl<G2, C2> From<&RelaxedOvaInstanceVar<G2, C2>> for OvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
{
    fn from(U: &RelaxedOvaInstanceVar<G2, C2>) -> Self {
        Self {
            commitment: U.commitment.clone(),
            X: U.X.clone(),
            _commitment_scheme: PhantomData,
        }
    }
}

#[must_use]
#[derive(Debug)]
pub struct RelaxedOvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
{
    /// Commitment to witness.
    pub commitment: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    /// Public input of relaxed instance.
    pub X: Vec<NonNativeFieldVar<G2::ScalarField, G2::BaseField>>,

    _commitment_scheme: PhantomData<C2>,
}

impl<G2, C2> Clone for RelaxedOvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
{
    fn clone(&self) -> Self {
        Self {
            commitment: self.commitment.clone(),
            X: self.X.clone(),
            _commitment_scheme: self._commitment_scheme,
        }
    }
}

impl<G2, C2> R1CSVar<G2::BaseField> for RelaxedOvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    type Value = RelaxedOvaInstance<G2, C2>;

    fn cs(&self) -> ConstraintSystemRef<G2::BaseField> {
        self.X
            .iter()
            .fold(ConstraintSystemRef::None, |cs, x| cs.or(x.cs()))
            .or(self.commitment.cs())
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let commitment = self.commitment.value()?.into();
        let X = self.X.value()?;

        Ok(RelaxedOvaInstance { commitment, X })
    }
}

impl<G2, C2> AllocVar<RelaxedOvaInstance<G2, C2>, G2::BaseField> for RelaxedOvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    fn new_variable<T: Borrow<RelaxedOvaInstance<G2, C2>>>(
        cs: impl Into<Namespace<G2::BaseField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let relaxed_ova_instance = f()?;
        let X = &relaxed_ova_instance.borrow().X;

        let commitment = ProjectiveVar::<G2, FpVar<G2::BaseField>>::new_variable(
            cs.clone(),
            || Ok(relaxed_ova_instance.borrow().commitment.into()),
            mode,
        )?;

        let X = X
            .iter()
            .map(|x| {
                NonNativeFieldVar::<G2::ScalarField, G2::BaseField>::new_variable(
                    cs.clone(),
                    || Ok(x),
                    mode,
                )
            })
            .collect::<Result<_, _>>()?;

        Ok(Self {
            commitment,
            X,
            _commitment_scheme: PhantomData,
        })
    }
}


impl<G2, C2> CondSelectGadget<G2::BaseField> for RelaxedOvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    fn conditionally_select(
        cond: &Boolean<G2::BaseField>,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError> {
        let commitment = cond.select(&true_value.commitment, &false_value.commitment)?;

        let X = true_value
            .X
            .iter()
            .zip(&false_value.X)
            .map(|(x1, x2)| cond.select(x1, x2))
            .collect::<Result<_, _>>()?;

        Ok(Self {
            commitment,
            X,
            _commitment_scheme: PhantomData,
        })
    }
}

impl<G2, C2> RelaxedOvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    pub fn fold(
        &self,
        instances: &[(
            (
                &OvaInstanceVar<G2, C2>,
                Option<&ProjectiveVar<G2, FpVar<G2::BaseField>>>,
            ),
            &'_ ProjectiveVar<G2, FpVar<G2::BaseField>>,
            &'_ NonNativeFieldVar<G2::ScalarField, G2::BaseField>,
            &'_ [Boolean<G2::BaseField>],
        )],
    ) -> Result<Self, SynthesisError> {
        let mut commitment = self.commitment.clone();
        let mut X: Vec<NonNativeFieldMulResultVar<_, _>> = self
            .X
            .iter()
            .map(NonNativeFieldMulResultVar::from)
            .collect();

        for ((U, _), commitment_T, r, r_bits) in instances {
            let res = U.commitment.clone() + *commitment_T;
            commitment += res.scalar_mul_le(r_bits.iter()).unwrap();
            for (x1, x2) in X.iter_mut().zip(&U.X) {
                *x1 += x2.mul_without_reduce(r)?;
            }
        }

        let X = X
            .iter()
            .map(NonNativeFieldMulResultVar::reduce)
            .collect::<Result<_, _>>()?;
        Ok(Self {
            commitment,
            X,
            _commitment_scheme: PhantomData,
        })
    }

    pub fn enforce_equal(&self, other: &Self) -> Result<(), SynthesisError> {
        // Enforce equality for the commitment field
        self.commitment.enforce_equal(&other.commitment)?;

        // Enforce equality for each element in the X vector
        if self.X.len() != other.X.len() {
            return Err(SynthesisError::AssignmentMissing);
        }

        for (x_self, x_other) in self.X.iter().zip(&other.X) {
            x_self.enforce_equal(x_other)?;
        }

        Ok(())
    }
}


macro_rules! parse_projective {
    ($X:ident) => {
        match &$X[..3] {
            [x, y, z, ..] => {
                let zero = Affine::<G1>::zero();
                let zero_x = NonNativeFieldVar::constant(zero.x);
                let zero_y = NonNativeFieldVar::constant(zero.y);
                let infinity = z.is_zero()?;

                let x = infinity.select(&zero_x, x)?;
                let y = infinity.select(&zero_y, y)?;

                let point = NonNativeAffineVar { x, y, infinity };
                $X = &$X[3..];
                point
            }
            _ => return Err(SynthesisError::Unsatisfiable),
        }
    };
}

impl<G2, C2> OvaInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    /// Parses `r, g_out` from the public input of the secondary circuit.
    pub fn parse_secondary_io<G1>(
        &self,
    ) -> Result<
        (
            NonNativeFieldVar<G1::BaseField, G1::ScalarField>,
            NonNativeFieldVar<G1::BaseField, G1::ScalarField>,
            NonNativeAffineVar<G1>,
            NonNativeAffineVar<G1>,
            NonNativeAffineVar<G1>,
        ),
        SynthesisError,
    >
    where
        G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    {
        let r = self.X[SecondaryCircuit::<G1>::NUM_IO - 2].clone();
        let flag = self.X[SecondaryCircuit::<G1>::NUM_IO - 1].clone();

        let mut X = &self.X[1..];
        let g1 = parse_projective!(X);
        let g2 = parse_projective!(X);
        let g_out = parse_projective!(X);

        let _ = X;

        Ok((flag, r, g1, g2, g_out))
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::Zero;
    use ark_pallas::PallasConfig;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_vesta::{Fq, VestaConfig};

    use crate::commitment::*;
    use crate::hash::pederson::PedersenCommitment;
    use crate::nova::cycle_fold::coprocessor::{setup_shape, synthesize};
    use crate::nova::cycle_fold::coprocessor_constraints::OvaInstanceVar;
    use crate::nova::cycle_fold::test::tests::get_random_circuit;

    #[test]
    fn initialisation_test() {
        // define a circuit on PallasConfig, the base are on Pallas (Fq)
        let c = get_random_circuit();

        // let G1 = PallasConfig and G2 = VestaConfig, circuit on G1::BaseField (Fr)
        let shape = setup_shape::<VestaConfig, PallasConfig>().unwrap();

        // pp which is Pedersen's pp, it has to be on opposite curve from the circuit ==> Vesta
        let pp = PedersenCommitment::<ark_vesta::Projective>::setup(shape.num_vars, b"test", &());

        // counting the number of constraints on the second curve
        let cs = ConstraintSystem::<Fq>::new_ref();
        println!("number of constraints on second curve: {}", cs.num_constraints());

        // define a running instance for second curve
        let (U, _W) = synthesize::<
            PallasConfig,
            VestaConfig,
            PedersenCommitment<ark_vesta::Projective>,
        >(c.clone(), &pp).unwrap();

        let U_var = OvaInstanceVar::new_variable(
            cs,
            || Ok(U.clone()),
            AllocationMode::Constant,
        ).unwrap();

        // test if the value is consistent
        let U_prime = U_var.value().unwrap();
        assert_eq!(U_prime, U);

        // parsing IO
        let (flag, r, g1, g2, g_out) = U_var.parse_secondary_io::<PallasConfig>().unwrap();
        assert_eq!(r.value().unwrap(), c.r);
        assert_eq!(!flag.value().unwrap().is_zero(), c.flag);
        assert_eq!(g1.value().unwrap(), c.g1);
        assert_eq!(g2.value().unwrap(), c.g2);
        assert_eq!(g_out.value().unwrap(), c.g_out);
    }
}


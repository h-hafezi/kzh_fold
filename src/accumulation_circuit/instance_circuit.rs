use std::borrow::Borrow;
use std::fmt::Debug;

use ark_ec::{AffineRepr, CurveConfig, CurveGroup};
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::{Field, One, PrimeField, Zero};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::{AffineVar, ProjectiveVar};
use ark_r1cs_std::R1CSVar;
use ark_r1cs_std::uint32::UInt32;
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};

use crate::accumulation::accumulator::AccInstance;
use crate::accumulation_circuit::affine_to_projective;
use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;

#[derive(Clone)]
/// the circuit is defined on scalar of G1
pub struct AccumulatorInstanceVar<G1>
where
    G1: SWCurveConfig + Clone,
    <G1 as CurveConfig>::ScalarField: PrimeField,
    <G1 as CurveConfig>::BaseField: PrimeField,
{
    // group points with base field G1::BaseField
    pub C_var: NonNativeAffineVar<G1>,
    pub T_var: NonNativeAffineVar<G1>,
    pub E_var: NonNativeAffineVar<G1>,
    // the field elements G1::ScalarField
    pub b_var: FpVar<G1::ScalarField>,
    pub c_var: FpVar<G1::ScalarField>,
    pub y_var: FpVar<G1::ScalarField>,
    pub z_b_var: FpVar<G1::ScalarField>,
    pub z_c_var: FpVar<G1::ScalarField>,
}


impl<G1> AccumulatorInstanceVar<G1>
where
    G1: SWCurveConfig + Clone,
    <G1 as CurveConfig>::ScalarField: PrimeField,
    <G1 as CurveConfig>::BaseField: PrimeField,
{
    pub(crate) fn cs(&self) -> ConstraintSystemRef<G1::ScalarField> {
        self.C_var.cs().or(self.T_var.cs())
            .or(self.E_var.cs())
            .or(self.b_var.cs())
            .or(self.c_var.cs())
            .or(self.y_var.cs())
            .or(self.z_b_var.cs())
            .or(self.z_c_var.cs())
    }

    pub(crate) fn value<E>(&self) -> Result<AccInstance<E>, SynthesisError>
    where
        E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField>,
    {
        Ok(AccInstance {
            C: self.C_var.value().unwrap().into(),
            T: self.T_var.value().unwrap().into(),
            E: self.E_var.value().unwrap().into(),
            b: self.b_var.value().unwrap(),
            c: self.c_var.value().unwrap(),
            y: self.y_var.value().unwrap(),
            z_b: self.z_b_var.value().unwrap(),
            z_c: self.z_c_var.value().unwrap(),
        })
    }
}


impl<G1, E> AllocVar<AccInstance<E>, <G1 as CurveConfig>::ScalarField> for AccumulatorInstanceVar<G1>
where
    G1: SWCurveConfig + Clone,
    <G1 as CurveConfig>::ScalarField: PrimeField,
    <G1 as CurveConfig>::BaseField: PrimeField,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField>,
{
    fn new_variable<T: Borrow<AccInstance<E>>>(
        cs: impl Into<Namespace<<G1 as CurveConfig>::ScalarField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        let C_var = NonNativeAffineVar::new_variable(
            ns!(cs, "C"),
            || circuit.map(|e| affine_to_projective(e.C)),
            mode,
        ).unwrap();

        let T_var = NonNativeAffineVar::new_variable(
            ns!(cs, "T"),
            || circuit.map(|e| affine_to_projective(e.T)),
            mode,
        ).unwrap();

        let E_var = NonNativeAffineVar::new_variable(
            ns!(cs, "E"),
            || circuit.map(|e| affine_to_projective(e.E)),
            mode,
        ).unwrap();

        let b_var = FpVar::new_variable(
            ns!(cs, "b"),
            || circuit.map(|e| e.b),
            mode,
        ).unwrap();

        let c_var = FpVar::new_variable(
            ns!(cs, "c"),
            || circuit.map(|e| e.c),
            mode,
        ).unwrap();

        let y_var = FpVar::new_variable(
            ns!(cs, "y"),
            || circuit.map(|e| e.y),
            mode,
        ).unwrap();

        let z_b_var = FpVar::new_variable(
            ns!(cs, "z_b"),
            || circuit.map(|e| e.z_b),
            mode,
        ).unwrap();

        let z_c_var = FpVar::new_variable(
            ns!(cs, "z_c"),
            || circuit.map(|e| e.z_c),
            mode,
        ).unwrap();

        Ok(AccumulatorInstanceVar {
            C_var,
            T_var,
            E_var,
            b_var,
            c_var,
            y_var,
            z_b_var,
            z_c_var,
        })
    }
}


#[cfg(test)]
pub mod tests {
    use std::fmt::Debug;

    use ark_ec::AffineRepr;
    use ark_ff::AdditiveGroup;
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::accumulation::accumulator::AccInstance;
    use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
    use crate::constant_for_curves::{E, ScalarField};

    #[test]
    fn initialisation_test() {
        // build an instance of AccInstanceCircuit
        let instance = AccInstance::<E> {
            C: <E as ark_ec::pairing::Pairing>::G1Affine::rand(&mut thread_rng()),
            T: <E as ark_ec::pairing::Pairing>::G1Affine::rand(&mut thread_rng()),
            E: <E as ark_ec::pairing::Pairing>::G1Affine::rand(&mut thread_rng()),
            b: ScalarField::rand(&mut thread_rng()),
            c: ScalarField::rand(&mut thread_rng()),
            y: ScalarField::rand(&mut thread_rng()),
            z_b: ScalarField::rand(&mut thread_rng()),
            z_c: ScalarField::rand(&mut thread_rng()),
        };

        // a constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // make a circuit_var
        let circuit_var = AccumulatorInstanceVar::new_variable(cs, || Ok(instance.clone()), AllocationMode::Constant).unwrap();
        // get its value and assert its equal to the original instance
        let c = circuit_var.value().unwrap();
        assert_eq!(c, instance, "the value function doesn't work well");
    }
}

use std::borrow::Borrow;
use std::fmt::Debug;

use ark_ec::{AffineRepr, CurveConfig, CurveGroup};
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::{Field, PrimeField, Zero};
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
use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccumulatorInstanceCircuit<G1>
where
    G1: SWCurveConfig + Clone,
    G1::ScalarField: PrimeField,
{
    pub C: Projective<G1>,
    pub T: Projective<G1>,
    // if non-relaxed instance then it should be zero
    pub E: Projective<G1>,
    pub b: G1::ScalarField,
    pub c: G1::ScalarField,
    pub y: G1::ScalarField,
    pub z_b: G1::ScalarField,
    pub z_c: G1::ScalarField,
}

impl<G1> AccumulatorInstanceCircuit<G1>
where
    G1: SWCurveConfig + Clone,
    G1::ScalarField: PrimeField,
{
    #[inline(always)]
    /// Returns if the error term E is zero
    pub fn is_fresh(&self) -> bool {
        self.E.is_zero()
    }

    #[inline(always)]
    /// Returns if the error term E is non-zero
    pub fn is_relaxed(&self) -> bool {
        !self.E.is_zero()
    }
}

#[derive(Clone)]
/// the circuit is defined on scalar of G1
pub struct AccumulatorInstanceCircuitVar<G1>
where
    G1: SWCurveConfig + Clone,
    G1::ScalarField: PrimeField,
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


impl<G1: SWCurveConfig + Clone> AccumulatorInstanceCircuitVar<G1>
where
    G1: SWCurveConfig,
    G1::ScalarField: PrimeField,
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

    pub(crate) fn value(&self) -> Result<AccumulatorInstanceCircuit<G1>, SynthesisError> {
        Ok(AccumulatorInstanceCircuit {
            C: self.C_var.value().unwrap(),
            T: self.T_var.value().unwrap(),
            E: self.E_var.value().unwrap(),
            b: self.b_var.value().unwrap(),
            c: self.c_var.value().unwrap(),
            y: self.y_var.value().unwrap(),
            z_b: self.z_b_var.value().unwrap(),
            z_c: self.z_c_var.value().unwrap(),
        })
    }
}


impl<G1> AllocVar<AccumulatorInstanceCircuit<G1>, G1::ScalarField> for AccumulatorInstanceCircuitVar<G1>
where
    G1: SWCurveConfig + Clone,
    G1::ScalarField: PrimeField,
    <G1 as CurveConfig>::BaseField: PrimeField,
{
    fn new_variable<T: Borrow<AccumulatorInstanceCircuit<G1>>>(
        cs: impl Into<Namespace<G1::ScalarField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        let C_var = NonNativeAffineVar::new_variable(
            ns!(cs, "C"),
            || circuit.map(|e| e.C),
            mode,
        ).unwrap();

        let T_var = NonNativeAffineVar::new_variable(
            ns!(cs, "T"),
            || circuit.map(|e| e.T),
            mode,
        ).unwrap();

        let E_var = NonNativeAffineVar::new_variable(
            ns!(cs, "E"),
            || circuit.map(|e| e.E),
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

        Ok(AccumulatorInstanceCircuitVar {
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

    use ark_ec::short_weierstrass::Projective;
    use ark_ff::{AdditiveGroup, Zero};
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::accumulation::accumulator::AccInstance;
    use crate::accumulation_circuit::acc_instance_circuit::{AccumulatorInstanceCircuit, AccumulatorInstanceCircuitVar};
    use crate::constant_for_curves::{E, G1, ScalarField};

    pub fn accumulator_instance_to_circuit(acc_instance: AccInstance<E>) -> AccumulatorInstanceCircuit<G1> {
        AccumulatorInstanceCircuit {
            C: acc_instance.C.into(),
            T: acc_instance.T.into(),
            E: acc_instance.E.into(),
            b: acc_instance.b,
            c: acc_instance.c,
            y: acc_instance.y,
            z_b: acc_instance.z_b,
            z_c: acc_instance.z_c,
        }
    }

    pub fn circuit_to_accumulator_instance(circuit: AccumulatorInstanceCircuit<G1>) -> AccInstance<E> {
        AccInstance {
            C: circuit.C.into(),
            T: circuit.T.into(),
            E: circuit.E.into(),
            b: circuit.b,
            c: circuit.c,
            y: circuit.y,
            z_b: circuit.z_b,
            z_c: circuit.z_c,
        }
    }


    #[test]
    fn initialisation_test() {
        // build an instance of AccInstanceCircuit
        let instance = AccumulatorInstanceCircuit::<G1> {
            C: Projective::zero(),
            T: Projective::zero(),
            E: Projective::zero(),
            b: ScalarField::ZERO,
            c: ScalarField::ZERO,
            y: ScalarField::ZERO,
            z_b: ScalarField::ZERO,
            z_c: ScalarField::ZERO,
        };

        // a constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // make a circuit_var
        let circuit_var = AccumulatorInstanceCircuitVar::new_variable(cs, || Ok(instance.clone()), AllocationMode::Constant).unwrap();
        // get its value and assert its equal to the original instance
        let c = circuit_var.value().unwrap();
        assert!(c == instance, "the value function doesn't work well");
    }

    #[test]
    fn test_conversion_functions() {
        // build an instance of AccInstanceCircuit
        let circuit = AccumulatorInstanceCircuit::<G1> {
            C: Projective::rand(&mut thread_rng()),
            T: Projective::rand(&mut thread_rng()),
            E: Projective::rand(&mut thread_rng()),
            b: ScalarField::rand(&mut thread_rng()),
            c: ScalarField::rand(&mut thread_rng()),
            y: ScalarField::rand(&mut thread_rng()),
            z_b: ScalarField::rand(&mut thread_rng()),
            z_c: ScalarField::rand(&mut thread_rng()),
        };

        let acc_instance = circuit_to_accumulator_instance(circuit.clone());

        let circuit_new = accumulator_instance_to_circuit(acc_instance);

        assert!(circuit_new == circuit, "equality problem: circuits are not equal");
    }
}

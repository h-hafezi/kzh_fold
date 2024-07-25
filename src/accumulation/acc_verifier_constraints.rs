use std::borrow::Borrow;
use ark_ec::CurveConfig;
use ark_ec::short_weierstrass::SWCurveConfig;
use ark_ff::{Field, PrimeField};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_relations::ns;
use ark_relations::r1cs::{Namespace, SynthesisError};
use crate::accumulation::acc_instance_constraints::{AccumulatorInstance, AccumulatorInstanceVar};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccumulatorVerifier<G1>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    FpVar<
        <G1 as CurveConfig>::BaseField
    >: FieldVar<
        <G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField
    >,
{
    instance: AccumulatorInstance<G1>,
    acc: AccumulatorInstance<G1>,
}

#[derive(Clone)]
pub struct AccumulatorVerifierVar<G1>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    FpVar<
        <G1 as CurveConfig>::BaseField
    >: FieldVar<
        <G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField
    >,
{
    instance: AccumulatorInstanceVar<G1>,
    acc: AccumulatorInstanceVar<G1>,
}

impl<G1> AllocVar<AccumulatorVerifier<G1>, G1::BaseField> for AccumulatorVerifierVar<G1>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    FpVar<
        <G1 as CurveConfig>::BaseField
    >: FieldVar<
        <G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField
    >,
{
    fn new_variable<T: Borrow<AccumulatorVerifier<G1>>>(cs: impl Into<Namespace<G1::BaseField>>, f: impl FnOnce() -> Result<T, SynthesisError>, mode: AllocationMode) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        let instance = AccumulatorInstanceVar::new_variable(
            ns!(cs, "instance"),
            || circuit.map(|e| e.instance.clone()),
            mode,
        ).unwrap();

        let acc = AccumulatorInstanceVar::new_variable(
            ns!(cs, "acc"),
            || circuit.map(|e| e.acc.clone()),
            mode,
        ).unwrap();

        Ok(AccumulatorVerifierVar {
            instance,
            acc,
        })
    }
}

impl<G1: SWCurveConfig> AccumulatorVerifierVar<G1>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    FpVar<
        <G1 as CurveConfig>::BaseField
    >: FieldVar<
        <G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField
    >,
{
    pub fn accumulate(&self) {
        // Poseidon hash
        // scalar multiplication: linear combination of C
        // scalar multiplication: linear combination of T
        // scalar multiplication: linear combination of E
        // Non-native field operation: linear combination of b
        // Non-native field operation: linear combination of c
        // Non-native field operation: linear combination of y
        // Non-native field operation: linear combination of z_b
        // Non-native field operation: linear combination of z_c
        // equality assertion that z_b = b^n-1 for the first instance
        // equality assertion that z_c = c^m-1 for the first instance
    }
}
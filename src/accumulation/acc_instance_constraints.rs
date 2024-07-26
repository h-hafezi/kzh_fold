use std::borrow::Borrow;
use std::fmt::Debug;
use ark_ec::CurveConfig;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::{Field, PrimeField};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::groups::curves::short_weierstrass::{AffineVar, ProjectiveVar};
use ark_r1cs_std::R1CSVar;
use ark_r1cs_std::uint32::UInt32;
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccumulatorInstance<G1>
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
    pub C: Projective<G1>,
    pub T: Projective<G1>,
    pub E: Projective<G1>,
    pub b: G1::BaseField,
    pub c: G1::BaseField,
    pub y: G1::BaseField,
    pub z_b: G1::BaseField,
    pub z_c: G1::BaseField,
    // these are constant values
    pub n: u32,
    pub m: u32,
}

#[derive(Clone)]
pub struct AccumulatorInstanceVar<G1>
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
    pub C_var: ProjectiveVar<G1, FpVar<G1::BaseField>>,
    pub T_var: ProjectiveVar<G1, FpVar<G1::BaseField>>,
    pub E_var: ProjectiveVar<G1, FpVar<G1::BaseField>>,
    pub b_var: FpVar<G1::BaseField>,
    pub c_var: FpVar<G1::BaseField>,
    pub y_var: FpVar<G1::BaseField>,
    pub z_b_var: FpVar<G1::BaseField>,
    pub z_c_var: FpVar<G1::BaseField>,
    // these are constant values
    pub n: u32,
    pub m: u32,
}


impl<G1: SWCurveConfig + Clone> AccumulatorInstanceVar<G1>
where
    FpVar<<G1 as CurveConfig>::BaseField>: FieldVar<<G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField>,
    <G1 as CurveConfig>::BaseField: PrimeField,
{
    pub(crate) fn cs(&self) -> ConstraintSystemRef<G1::BaseField> {
        self.C_var.cs().or(self.T_var.cs())
            .or(self.b_var.cs())
            .or(self.c_var.cs())
            .or(self.y_var.cs())
            .or(self.z_b_var.cs())
            .or(self.z_c_var.cs())
            .or(self.E_var.cs())
    }

    pub(crate) fn value(&self) -> Result<AccumulatorInstance<G1>, SynthesisError> {
        Ok(AccumulatorInstance {
            C: self.C_var.value().unwrap(),
            T: self.T_var.value().unwrap(),
            b: self.b_var.value().unwrap(),
            c: self.c_var.value().unwrap(),
            y: self.y_var.value().unwrap(),
            z_b: self.z_b_var.value().unwrap(),
            z_c: self.z_c_var.value().unwrap(),
            n: self.n,
            E: self.E_var.value().unwrap(),
            m: self.m,
        })
    }
}

impl<G1> AllocVar<AccumulatorInstance<G1>, G1::BaseField> for AccumulatorInstanceVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: Field<BasePrimeField=<G1 as CurveConfig>::BaseField> + PrimeField,
    FpVar<G1::BaseField>: FieldVar<G1::BaseField, <G1::BaseField as Field>::BasePrimeField>,
    G1: Clone,
{
    fn new_variable<T: Borrow<AccumulatorInstance<G1>>>(
        cs: impl Into<Namespace<G1::BaseField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        let C_var = ProjectiveVar::new_variable(
            ns!(cs, "C"),
            || circuit.map(|e| e.C),
            mode,
        ).unwrap();

        let T_var = ProjectiveVar::new_variable(
            ns!(cs, "T"),
            || circuit.map(|e| e.T),
            mode,
        ).unwrap();

        let E_var = ProjectiveVar::new_variable(
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

        Ok(AccumulatorInstanceVar {
            C_var,
            T_var,
            E_var,
            b_var,
            c_var,
            y_var,
            z_b_var,
            z_c_var,
            n: circuit.unwrap().n,
            m: circuit.unwrap().m,
        })
    }
}


#[cfg(test)]
mod tests {
    use std::fmt::{Debug, Formatter};
    use ark_ff::{AdditiveGroup, Field, PrimeField, Zero};
    use ark_ec::short_weierstrass::{SWCurveConfig, Projective};
    use ark_std::UniformRand;
    use ark_bn254::g1::Config;
    use ark_bn254::{Fr, Fq};
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use itertools::equal;
    use crate::accumulation::acc_instance_constraints::{AccumulatorInstance, AccumulatorInstanceVar};

    #[test]
    fn initialisation_test() {
        // build an instance of AccInstanceCircuit
        let instance = AccumulatorInstance::<Config> {
            C: Projective::zero(),
            T: Projective::zero(),
            E: Projective::zero(),
            b: Fq::ZERO,
            c: Fq::ZERO,
            y: Fq::ZERO,
            z_b: Fq::ZERO,
            z_c: Fq::ZERO,
            n: 0u32,
            m: 0u32,
        };
        // a constraint system
        let cs = ConstraintSystem::<Fq>::new_ref();
        // make a circuit_var
        let circuit_var = AccumulatorInstanceVar::new_variable(cs, || Ok(instance.clone()), AllocationMode::Constant).unwrap();
        // get its value and assert its equal to the original instance
        let c = circuit_var.value().unwrap();
        if !(c.T == instance.T && c.C == instance.C && c.E == instance.E && c.b == instance.b && c.c == instance.c
            && c.y == instance.y && c.z_b == instance.z_b && c.z_c == instance.z_c && c.m == instance.m && c.n == instance.n) {
            panic!("the value function doesn't work well")
        }
    }
}

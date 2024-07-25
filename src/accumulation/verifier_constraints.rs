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
pub struct AccInstanceCircuit<G1: SWCurveConfig>
where
    FpVar<<G1 as CurveConfig>::BaseField>: FieldVar<<G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField>,
    <G1 as CurveConfig>::BaseField: PrimeField,
{
    C: Projective<G1>,
    T: Projective<G1>,
    E: Projective<G1>,
    b: G1::BaseField,
    c: G1::BaseField,
    y: G1::BaseField,
    z_b: G1::BaseField,
    z_c: G1::BaseField,
    beta: G1::BaseField,
    n: u32,
    m: u32,
}

#[derive(Clone)]
pub struct AccInstanceCircuitVar<G1: SWCurveConfig>
where
    FpVar<<G1 as CurveConfig>::BaseField>: FieldVar<<G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField>,
    <G1 as CurveConfig>::BaseField: PrimeField,
{
    C_var: ProjectiveVar<G1, FpVar<G1::BaseField>>,
    T_var: ProjectiveVar<G1, FpVar<G1::BaseField>>,
    E_var: ProjectiveVar<G1, FpVar<G1::BaseField>>,
    b_var: FpVar<G1::BaseField>,
    c_var: FpVar<G1::BaseField>,
    y_var: FpVar<G1::BaseField>,
    z_b_var: FpVar<G1::BaseField>,
    z_c_var: FpVar<G1::BaseField>,
    beta_var: FpVar<G1::BaseField>,
    n_var: UInt32<G1::BaseField>,
    m_var: UInt32<G1::BaseField>,
}


impl<G1: SWCurveConfig + Clone + Eq> AccInstanceCircuitVar<G1>
where
    FpVar<<G1 as CurveConfig>::BaseField>: FieldVar<<G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField>,
    <G1 as CurveConfig>::BaseField: PrimeField,
{
    fn cs(&self) -> ConstraintSystemRef<G1::BaseField> {
        self.C_var.cs().or(self.T_var.cs())
            .or(self.b_var.cs())
            .or(self.c_var.cs())
            .or(self.y_var.cs())
            .or(self.n_var.cs())
            .or(self.n_var.cs())
            .or(self.m_var.cs())
            .or(self.z_b_var.cs())
            .or(self.z_c_var.cs())
            .or(self.E_var.cs())
            .or(self.beta_var.cs())
    }

    fn value(&self) -> Result<AccInstanceCircuit<G1>, SynthesisError> {
        Ok(AccInstanceCircuit {
            C: self.C_var.value().unwrap(),
            T: self.T_var.value().unwrap(),
            b: self.b_var.value().unwrap(),
            c: self.c_var.value().unwrap(),
            y: self.y_var.value().unwrap(),
            n: self.n_var.value().unwrap(),
            m: self.m_var.value().unwrap(),
            z_b: self.z_b_var.value().unwrap(),
            z_c: self.z_c_var.value().unwrap(),
            E: self.E_var.value().unwrap(),
            beta: self.beta_var.value().unwrap(),
        })
    }
}

impl<G1> AllocVar<AccInstanceCircuit<G1>, G1::BaseField> for AccInstanceCircuitVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: Field<BasePrimeField=<G1 as CurveConfig>::BaseField> + PrimeField,
    FpVar<G1::BaseField>: FieldVar<G1::BaseField, <G1::BaseField as Field>::BasePrimeField>,
{
    fn new_variable<T: Borrow<AccInstanceCircuit<G1>>>(
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

        let n_var = UInt32::new_variable(
            ns!(cs, "n"),
            || circuit.map(|e| e.n),
            mode,
        ).unwrap();

        let m_var = UInt32::new_variable(
            ns!(cs, "m"),
            || circuit.map(|e| e.m),
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

        let beta_var = FpVar::new_variable(
            ns!(cs, "beta"),
            || circuit.map(|e| e.beta),
            mode,
        ).unwrap();

        Ok(AccInstanceCircuitVar {
            C_var,
            T_var,
            E_var,
            b_var,
            c_var,
            y_var,
            z_b_var,
            z_c_var,
            beta_var,
            n_var,
            m_var,
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
    use crate::accumulation::verifier_constraints::{AccInstanceCircuit, AccInstanceCircuitVar};

    #[test]
    fn initialisation_test() {
        // build an instance of AccInstanceCircuit
        let instance = AccInstanceCircuit::<Config> {
            C: Projective::zero(),
            T: Projective::zero(),
            E: Projective::zero(),
            b: Fq::ZERO,
            c: Fq::ZERO,
            y: Fq::ZERO,
            z_b: Fq::ZERO,
            z_c: Fq::ZERO,
            beta: Fq::ZERO,
            n: 0u32,
            m: 0u32,
        };
        // a constraint system
        let cs = ConstraintSystem::<Fq>::new_ref();
        // make a circuit_var
        let circuit_var = AccInstanceCircuitVar::new_variable(cs, || Ok(instance.clone()), AllocationMode::Constant).unwrap();
        // get its value and assert its equal to the original instance
        let c = circuit_var.value().unwrap();
        if !(c.T == instance.T && c.C == instance.C && c.E == instance.E && c.b == instance.b && c.c == instance.c && c.n == instance.n
            && c.y == instance.y && c.z_b == instance.z_b && c.z_c == instance.z_c && c.beta == instance.beta && c.m == instance.m) {
            panic!("the value function doesn't work well")
        }
    }
}

use std::borrow::Borrow;
use std::ops::Add;
use ark_ec::CurveConfig;
use ark_ec::short_weierstrass::SWCurveConfig;
use ark_ff::{BigInteger64, Field, PrimeField};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::{R1CSVar, ToBitsGadget};
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::groups::CurveVar;
use ark_relations::ns;
use ark_relations::r1cs::{Namespace, SynthesisError};
use ark_std::UniformRand;
use rand::thread_rng;
use crate::accumulation_circuit::acc_instance_constraints::{AccumulatorInstance, AccumulatorInstanceVar};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccumulatorVerifier<G1>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    FpVar<
        <G1 as CurveConfig>::BaseField
    >: FieldVar<
        <G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField
    >,
{
    pub instance: AccumulatorInstance<G1>,
    pub acc: AccumulatorInstance<G1>,
}

#[derive(Clone)]
pub struct AccumulatorVerifierVar<G1>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    FpVar<
        <G1 as CurveConfig>::BaseField
    >: FieldVar<
        <G1 as CurveConfig>::BaseField,
        <<G1 as CurveConfig>::BaseField as Field>::BasePrimeField
    >,
{
    pub instance: AccumulatorInstanceVar<G1>,
    pub acc: AccumulatorInstanceVar<G1>,
}

impl<G1> AllocVar<AccumulatorVerifier<G1>, G1::BaseField> for AccumulatorVerifierVar<G1>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
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
        let beta_fr = NonNativeFieldVar::new_variable(
            ns!(self.acc.cs(), "beta"),
            || Ok(G1::ScalarField::rand(&mut thread_rng())),
            AllocationMode::Input,
        ).unwrap();
        let beta_bits = beta_fr.to_bits_le().unwrap();

        // Non-native scalar multiplication: linear combination of C
        let _ = &self.acc.C_var + &self.instance.C_var.scalar_mul_le(beta_bits.iter()).unwrap();

        // Non-native scalar multiplication: linear combination of T
        let _ = &self.acc.T_var + &self.instance.T_var.scalar_mul_le(beta_bits.iter()).unwrap();

        // Non-native scalar multiplication: linear combination of E
        let _ = &self.acc.E_var + &self.instance.E_var.scalar_mul_le(beta_bits.iter()).unwrap();

        // Native field operation: linear combination of b
        let _ = &self.acc.b_var + &beta_fr * &self.instance.b_var;

        // Native field operation: linear combination of c
        let _ = &self.acc.c_var + &beta_fr * &self.instance.c_var;

        // Native field operation: linear combination of y
        let _ = &self.acc.y_var + &beta_fr * &self.instance.y_var;

        // Native field operation: linear combination of z_b
        let _ = &self.acc.z_b_var + &beta_fr * &self.instance.z_b_var;

        // Native field operation: linear combination of z_c
        let _ = &self.acc.z_c_var + &beta_fr * &self.instance.z_c_var;

        // Native field operation: equality assertion that z_b = b^n-1 for the first instance
        let n = BigInteger64::from(self.acc.n);
        let _ = self.acc.b_var.pow_by_constant(n.as_ref()).expect("TODO: panic message");

        // Native field operation: equality assertion that z_c = c^m-1 for the first instance
        let m = BigInteger64::from(self.acc.m);
        let _ = self.acc.c_var.pow_by_constant(m.as_ref()).expect("TODO: panic message");
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
    use rand::thread_rng;
    use crate::accumulation_circuit::acc_instance_constraints::{AccumulatorInstance, AccumulatorInstanceVar};
    use crate::accumulation_circuit::acc_verifier_constraints::AccumulatorVerifierVar;

    #[test]
    fn initialisation_test() {
        // build an instance of AccInstanceCircuit
        let instance = AccumulatorInstance::<Config> {
            C: Projective::rand(&mut thread_rng()),
            T: Projective::rand(&mut thread_rng()),
            E: Projective::rand(&mut thread_rng()),
            b: Fr::rand(&mut thread_rng()),
            c: Fr::rand(&mut thread_rng()),
            y: Fr::rand(&mut thread_rng()),
            z_b: Fr::rand(&mut thread_rng()),
            z_c: Fr::rand(&mut thread_rng()),
            n: 1000u32,
            m: 1000u32,
        };

        let acc = AccumulatorInstance::<Config> {
            C: Projective::rand(&mut thread_rng()),
            T: Projective::rand(&mut thread_rng()),
            E: Projective::rand(&mut thread_rng()),
            b: Fr::rand(&mut thread_rng()),
            c: Fr::rand(&mut thread_rng()),
            y: Fr::rand(&mut thread_rng()),
            z_b: Fr::rand(&mut thread_rng()),
            z_c: Fr::rand(&mut thread_rng()),
            n: 1000u32,
            m: 1000u32,
        };
        // a constraint system
        let cs = ConstraintSystem::<Fq>::new_ref();

        // make a circuit_var
        let instance_var = AccumulatorInstanceVar::new_variable(
            cs.clone(),
            || Ok(instance.clone()),
            AllocationMode::Witness
        ).unwrap();

        let acc_var = AccumulatorInstanceVar::new_variable(cs.clone(),
                                                           || Ok(acc.clone()),
                                                           AllocationMode::Witness
        ).unwrap();

        let verifier = AccumulatorVerifierVar { instance: instance_var, acc: acc_var };
        println!("before: {}", cs.num_constraints());
        verifier.accumulate();
        println!("after: {}", cs.num_constraints());
    }
}

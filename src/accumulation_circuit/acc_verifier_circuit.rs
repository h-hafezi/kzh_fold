use std::borrow::Borrow;
use std::ops::Add;

use ark_ec::CurveConfig;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::{BigInteger64, Field, PrimeField};
use ark_r1cs_std::{R1CSVar, ToBitsGadget};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::groups::CurveVar;
use ark_relations::ns;
use ark_relations::r1cs::{Namespace, SynthesisError};
use ark_std::UniformRand;
use rand::thread_rng;

use crate::accumulation_circuit::acc_instance_circuit::{AccumulatorInstance, AccumulatorInstanceVar};
use crate::gadgets::r1cs::R1CSInstance;
use crate::nova::commitment::CommitmentScheme;
use crate::nova::cycle_fold::coprocessor::Circuit as SecondCircuit;
use crate::nova::cycle_fold::coprocessor_constraints::R1CSInstanceVar;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccumulatorVerifier<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    /// the randomness used for taking linear combination
    pub r: G1::BaseField,
    /// auxiliary input which helps to have w without scalar multiplication
    // todo: make sure using G2, C2 is correct
    pub auxiliary_input_w: R1CSInstance<G2, C2>,
    /// the instance to be folded
    pub instance: AccumulatorInstance<G1>,
    /// the running accumulator
    pub acc: AccumulatorInstance<G1>,
}


#[derive(Clone)]
pub struct AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
{
    /// The auxiliary input is on base field of G2 which is equal to scalar field of G1
    pub auxiliary_input_w: R1CSInstanceVar<G2, C2>,
    /// the randomness used for taking linear combination
    pub r: FpVar<G1::BaseField>,
    pub instance: AccumulatorInstanceVar<G1>,
    pub acc: AccumulatorInstanceVar<G1>,
}

impl<G1, G2, C2> AllocVar<AccumulatorVerifier<G1, G2, C2>, G1::ScalarField> for AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>
{
    fn new_variable<T: Borrow<AccumulatorVerifier<G1, G2, C2>>>(cs: impl Into<Namespace<G1::ScalarField>>, f: impl FnOnce() -> Result<T, SynthesisError>, mode: AllocationMode) -> Result<Self, SynthesisError> {
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

        let auxiliary_input_w = R1CSInstanceVar::new_variable(
            ns!(cs, "acc"),
            || circuit.map(|e| e.auxiliary_input_w.clone()),
            mode,
        ).unwrap();

        let r = FpVar::new_variable(
            ns!(cs, "acc"),
            || circuit.map(|e| e.r.clone()),
            mode,
        ).unwrap();


        Ok(AccumulatorVerifierVar {
            auxiliary_input_w,
            r,
            instance,
            acc,
        })
    }
}

impl<G1: SWCurveConfig, G2: SWCurveConfig, C2> AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>
{
    pub fn accumulate(&self) {
        // Poseidon hash
        let beta_fr = FpVar::new_variable(
            ns!(self.acc.cs(), "beta"),
            || Ok(G1::ScalarField::rand(&mut thread_rng())),
            AllocationMode::Input,
        ).unwrap();

        // Non-native scalar multiplication: linear combination of C
        //let _ = &self.acc.C_var + &self.instance.C_var.scalar_mul_le(beta_bits.iter()).unwrap();

        // Non-native scalar multiplication: linear combination of T
        //let _ = &self.acc.T_var + &self.instance.T_var.scalar_mul_le(beta_bits.iter()).unwrap();

        // Non-native scalar multiplication: linear combination of E
        //let _ = &self.acc.E_var + &self.instance.E_var.scalar_mul_le(beta_bits.iter()).unwrap();

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
    use std::fmt::Debug;

    use ark_bn254::g1::Config;
    use ark_pallas::{PallasConfig, Projective};
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use ark_vesta::{Fr, VestaConfig};
    use rand::thread_rng;

    use crate::accumulation_circuit::acc_instance_circuit::{AccumulatorInstance, AccumulatorInstanceVar};
    use crate::accumulation_circuit::acc_verifier_circuit::AccumulatorVerifierVar;
    use crate::hash::pederson::PedersenCommitment;
    use crate::nova::commitment::CommitmentScheme;
    use crate::nova::cycle_fold::coprocessor::{Circuit, setup_shape, synthesize};

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
        let cs = ConstraintSystem::<Fr>::new_ref();

        // make a circuit_var
        let instance_var = AccumulatorInstanceVar::new_variable(
            cs.clone(),
            || Ok(instance.clone()),
            AllocationMode::Witness,
        ).unwrap();

        let acc_var = AccumulatorInstanceVar::new_variable(
            cs.clone(),
            || Ok(acc.clone()),
            AllocationMode::Witness,
        ).unwrap();

        let r = Fr::from(2u128);
        let r_var = FpVar::new_variable(
            cs.clone(),
            || Ok(r.clone()),
            AllocationMode::Witness,
        ).unwrap();


        let c = {
            let g1 = instance.C.clone();
            let g2 = acc.C.clone();
            let g_out = g1 * r + g2;
            Circuit {
                g1: instance.C,
                g2: acc.C,
                g_out,
                r,
            }
        };

        let auxiliary_input_w = {
            let shape = setup_shape::<PallasConfig, VestaConfig>().unwrap();
            let pp = PedersenCommitment::<ark_vesta::Projective>::setup(shape.num_vars, b"test", &());
            let (u, _) = synthesize::<
                PallasConfig,
                VestaConfig,
                PedersenCommitment<ark_vesta::Projective>,
            >(c, &pp).unwrap();
        };

        let verifier = AccumulatorVerifierVar { auxiliary_input_w, r: r_var, instance: instance_var, acc: acc_var };
        println!("before: {}", cs.num_constraints());
        verifier.accumulate();
        println!("after: {}", cs.num_constraints());
    }
}


use std::borrow::Borrow;
use std::marker::PhantomData;
use std::ops::Add;

use ark_ec::CurveConfig;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::{BigInteger64, Field, PrimeField};
use ark_r1cs_std::{R1CSVar, ToBitsGadget};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::groups::CurveVar;
use ark_relations::ns;
use ark_relations::r1cs::{Namespace, SynthesisError};
use ark_std::UniformRand;
use rand::thread_rng;

use crate::accumulation_circuit::acc_instance_circuit::{AccumulatorInstanceCircuit, AccumulatorInstanceCircuitVar};
use crate::constant_for_curves::{BaseField, ScalarField};
use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;
use crate::gadgets::non_native::util::{convert_field_one_to_field_two, non_native_to_fpvar};
use crate::gadgets::r1cs::R1CSInstance;
use crate::nova::commitment::CommitmentScheme;
use crate::nova::cycle_fold::coprocessor::{SecondaryCircuit as SecondaryCircuit, synthesize};
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
// this condition is needed for cycle of curves
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
{
    /// the randomness used for taking linear combination
    pub beta: G1::ScalarField,
    /// auxiliary input which helps to have C'' = (1-beta) * C + beta * C' without scalar multiplication
    pub auxiliary_input_C: R1CSInstance<G2, C2>,
    /// auxiliary input which helps to have T'' = (1-beta) * T + beta * T' without scalar multiplication
    pub auxiliary_input_T: R1CSInstance<G2, C2>,
    /// auxiliary input which helps to have E_{temp} = (1-beta) * E + beta * E' without scalar multiplication
    pub auxiliary_input_E_1: R1CSInstance<G2, C2>,
    /// auxiliary input which helps to have E'' = E_{temp} + beta * (1-beta) * Q without scalar multiplication
    pub auxiliary_input_E_2: R1CSInstance<G2, C2>,
    /// accumulation proof
    pub Q: Projective<G1>,
    /// the instance to be folded
    pub instance: AccumulatorInstanceCircuit<G1>,
    /// the running accumulator
    pub acc: AccumulatorInstanceCircuit<G1>,

    // these are constant values
    pub n: u32,
    pub m: u32,
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
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
{
    /// auxiliary input which helps to have C'' = (1-beta) * C + beta * C' without scalar multiplication
    pub auxiliary_input_C: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have T'' = (1-beta) * T + beta * T' without scalar multiplication
    pub auxiliary_input_T: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E_{temp} = (1-beta) * E + beta * E' without scalar multiplication
    pub auxiliary_input_E_1: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E'' = E_{temp} + beta * (1-beta) * Q without scalar multiplication
    pub auxiliary_input_E_2: R1CSInstanceVar<G2, C2>,

    /// the randomness used for taking linear combination and its non-native counterpart
    pub beta: FpVar<G1::ScalarField>,
    pub beta_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// accumulation proof
    pub Q: NonNativeAffineVar<G1>,

    pub instance: AccumulatorInstanceCircuitVar<G1>,
    pub acc: AccumulatorInstanceCircuitVar<G1>,

    // these are constant values
    pub n: u32,
    pub m: u32,
}

impl<G1, G2, C2> AllocVar<AccumulatorVerifier<G1, G2, C2>, G1::ScalarField> for AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
{
    fn new_variable<T: Borrow<AccumulatorVerifier<G1, G2, C2>>>(cs: impl Into<Namespace<G1::ScalarField>>, f: impl FnOnce() -> Result<T, SynthesisError>, mode: AllocationMode) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        let auxiliary_input_C = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_C"),
            || Ok(circuit.map(|e| e.auxiliary_input_C.clone()).unwrap()),
            mode,
        ).unwrap();

        let auxiliary_input_T = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_T"),
            || Ok(circuit.map(|e| e.auxiliary_input_T.clone()).unwrap()),
            mode,
        ).unwrap();

        let auxiliary_input_E_1 = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_1"),
            || Ok(circuit.map(|e| e.auxiliary_input_E_1.clone()).unwrap()),
            mode,
        ).unwrap();

        let auxiliary_input_E_2 = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_2"),
            || Ok(circuit.map(|e| e.auxiliary_input_E_2.clone()).unwrap()),
            mode,
        ).unwrap();


        let instance = AccumulatorInstanceCircuitVar::new_variable(
            ns!(cs, "instance"),
            || circuit.map(|e| e.instance.clone()),
            mode,
        ).unwrap();

        let acc = AccumulatorInstanceCircuitVar::new_variable(
            ns!(cs, "acc"),
            || circuit.map(|e| e.acc.clone()),
            mode,
        ).unwrap();

        let beta = FpVar::new_variable(
            ns!(cs, "beta"),
            || circuit.map(|e| e.beta.clone()),
            mode,
        ).unwrap();

        let beta_non_native = NonNativeFieldVar::new_variable(
            ns!(cs, "non native beta"),
            || circuit.map(|e| convert_field_one_to_field_two::<G1::ScalarField, G1::BaseField>(e.beta.clone())),
            mode,
        ).unwrap();

        let Q = NonNativeAffineVar::new_variable(
            ns!(cs, "Q"),
            || circuit.map(|e| e.Q),
            mode,
        ).unwrap();

        Ok(AccumulatorVerifierVar {
            auxiliary_input_C,
            auxiliary_input_T,
            auxiliary_input_E_1,
            auxiliary_input_E_2,
            beta,
            beta_non_native,
            Q,
            instance,
            acc,
            n: 0,
            m: 0,
        })
    }
}

/// Here we assume instance to be A.X and acc to be A.X'
impl<G1: SWCurveConfig, G2: SWCurveConfig, C2> AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
{
    pub fn accumulate(&self) {
        // checking beta and non_native beta are consistent
        let beta_ = non_native_to_fpvar(&self.beta_non_native);
        self.beta.enforce_equal(&beta_).expect("error while enforcing equality");

        // compute Poseidon hash and make sure it's consistent with input beta

        // Non-native scalar multiplication: linear combination of C
        let (flag,
            r,
            g1,
            g2,
            new_C
        ) = self.auxiliary_input_C.parse_secondary_io().unwrap();
        // g1 == acc.C
        self.acc.C_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.C
        self.instance.C_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_non_native).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination of T
        let (flag,
            r,
            g1,
            g2,
            new_T
        ) = self.auxiliary_input_T.parse_secondary_io().unwrap();
        // g1 == acc.T
        self.acc.T_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.C
        self.instance.T_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_non_native).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination E_temp = (instance.E * (1-beta) + acc.E * beta)
        let (flag,
            r,
            g1,
            g2,
            E_temp
        ) = self.auxiliary_input_E_1.parse_secondary_io().unwrap();
        // g1 == acc.E
        self.acc.E_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.E
        self.instance.E_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_non_native).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination E'' = E_{temp} + (1-beta) * beta * Q
        let (flag,
            r,
            g1,
            g2,
            new_E
        ) = self.auxiliary_input_E_2.parse_secondary_io().unwrap();
        // g1 == Q
        g1.enforce_equal(&self.Q).expect("error while enforcing equality");
        // g2 == E_temp
        g2.enforce_equal(&E_temp).expect("error while enforcing equality");
        // enforce flag to be true
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        // check r to be equal to beta
        let beta_times_beta_minus_one = self.beta_non_native.clone() * (NonNativeFieldVar::one() - self.beta_non_native.clone());
        r.enforce_equal(&beta_times_beta_minus_one).expect("error while enforcing equality");


        // Native field operation: linear combination of b
        let _ = &self.acc.b_var + &self.beta * &self.instance.b_var;

        // Native field operation: linear combination of c
        let _ = &self.acc.c_var + &self.beta * &self.instance.c_var;

        // Native field operation: linear combination of y
        let _ = &self.acc.y_var + &self.beta * &self.instance.y_var;

        // Native field operation: linear combination of z_b
        let _ = &self.acc.z_b_var + &self.beta * &self.instance.z_b_var;

        // Native field operation: linear combination of z_c
        let _ = &self.acc.z_c_var + &self.beta * &self.instance.z_c_var;

        // Native field operation: equality assertion that z_b = b^n-1 for the first instance
        let n = BigInteger64::from(self.n);
        let z_b_ = self.acc.b_var.pow_by_constant(n.as_ref()).expect("error while enforcing equality");

        // Native field operation: equality assertion that z_c = c^m-1 for the first instance
        let m = BigInteger64::from(self.m);
        let z_c_ = self.acc.c_var.pow_by_constant(m.as_ref()).expect("error while enforcing equality");

        // Conditional check: if instance.E_var == 0, then enforce z_b_ == instance.z_b_var and z_c_ == instance.z_c_var
        let is_E_zero = &self.instance.E_var.infinity.clone();
        z_b_.conditional_enforce_equal(&self.instance.z_b_var, is_E_zero).expect("error while enforcing z_b equality under condition");
        z_c_.conditional_enforce_equal(&self.instance.z_c_var, is_E_zero).expect("error while enforcing z_c equality under condition");
    }
}


#[cfg(test)]
mod tests {
    use std::fmt::Debug;
    use ark_ec::short_weierstrass::Projective;
    use ark_ff::Field;
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
    use ark_relations::ns;
    use ark_relations::r1cs::ConstraintSystem;

    use crate::accumulation::accumulator::{Accumulator, AccumulatorTrait};
    use crate::accumulation::accumulator::tests::{get_satisfying_accumulator, get_srs};
    use crate::accumulation_circuit::acc_instance_circuit::AccumulatorInstanceCircuitVar;
    use crate::accumulation_circuit::acc_instance_circuit::tests::accumulator_instance_to_circuit;
    use crate::accumulation_circuit::acc_verifier_circuit::AccumulatorVerifierVar;
    use crate::constant_for_curves::{BaseField, G1, G2, ScalarField};
    use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;
    use crate::gadgets::non_native::util::convert_field_one_to_field_two;
    use crate::hash::pederson::PedersenCommitment;
    use crate::nova::commitment::CommitmentScheme;
    use crate::nova::cycle_fold::coprocessor::{SecondaryCircuit, setup_shape, synthesize};
    use crate::nova::cycle_fold::coprocessor_constraints::R1CSInstanceVar;

    type GrumpkinCurveGroup = ark_grumpkin::Projective;

    #[test]
    fn shitting() {}

    #[test]
    fn initialisation_test() {
        // specifying degrees of polynomials
        let n = 16;
        let m = 16;

        // get a random srs
        let srs = get_srs(n, m);

        // build an instance of AccInstanceCircuit
        let instance = get_satisfying_accumulator(&srs);

        let acc = get_satisfying_accumulator(&srs);

        // accumulate proof
        let (_, _, Q) = Accumulator::prove(&srs, &acc, &instance);
        let Q = Projective::new(Q.x, Q.y, BaseField::ONE);

        // build instance/acc circuits
        let instance = accumulator_instance_to_circuit(instance.instance);
        let acc = accumulator_instance_to_circuit(acc.instance);

        // a constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // make a circuit_var
        let instance_var = AccumulatorInstanceCircuitVar::new_variable(
            cs.clone(),
            || Ok(instance.clone()),
            AllocationMode::Witness,
        ).unwrap();

        let acc_var = AccumulatorInstanceCircuitVar::new_variable(
            cs.clone(),
            || Ok(acc.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // the randomness in different formats
        let beta_scalar = ScalarField::from(2u128);
        let beta_base = convert_field_one_to_field_two::<ScalarField, BaseField>(beta_scalar);
        let beta_var = FpVar::new_variable(
            cs.clone(),
            || Ok(beta_scalar.clone()),
            AllocationMode::Witness,
        ).unwrap();
        let beta_non_native = NonNativeFieldVar::new_variable(
            cs.clone(),
            || Ok(beta_base.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // the shape of the R1CS instance
        let shape = setup_shape::<G1, G2>().unwrap();

        // public parameters of Pedersen
        let pp = PedersenCommitment::<GrumpkinCurveGroup>::setup(shape.num_vars, b"test", &());

        let auxiliary_input_C = {
            // the circuit
            let c = {
                let g1 = acc.C.clone();
                let g2 = instance.C.clone();
                // C'' = beta * acc.C + (1 - beta) * instance.C
                let g_out = (g1 * beta_scalar) + (g2 * (ScalarField::ONE - beta_scalar));
                SecondaryCircuit {
                    g1: instance.C,
                    g2: acc.C,
                    g_out,
                    r: beta_base,
                    flag: false,
                }
            };

            let r1cs_instance = synthesize::<
                G1,
                G2,
                PedersenCommitment<GrumpkinCurveGroup>,
            >(c, &pp).unwrap().0;

            R1CSInstanceVar::new_variable(
                cs.clone(),
                || Ok(r1cs_instance),
                AllocationMode::Witness,
            ).unwrap()
        };

        let auxiliary_input_T = {
            // the circuit
            let c = {
                let g1 = acc.T.clone();
                let g2 = instance.T.clone();
                // T'' = beta * acc.T + (1 - beta) * instance.T
                let g_out = (g1 * beta_scalar) + (g2 * (ScalarField::ONE - beta_scalar));
                SecondaryCircuit {
                    g1: acc.T,
                    g2: instance.T,
                    g_out,
                    r: beta_base,
                    flag: false,
                }
            };

            let r1cs_instance = synthesize::<
                G1,
                G2,
                PedersenCommitment<GrumpkinCurveGroup>,
            >(c, &pp).unwrap().0;

            R1CSInstanceVar::new_variable(
                cs.clone(),
                || Ok(r1cs_instance),
                AllocationMode::Witness,
            ).unwrap()
        };

        let auxiliary_input_E_1 = {
            // the circuit
            let c = {
                let g1 = instance.E.clone();
                let g2 = acc.E.clone();
                // E_temp = beta * acc.T + (1 - beta) * instance.T
                let g_out = (g1 * beta_scalar) + (g2 * (ScalarField::ONE - beta_scalar));
                SecondaryCircuit {
                    g1: acc.E,
                    g2: instance.E,
                    g_out,
                    r: beta_base,
                    flag: false,
                }
            };

            let r1cs_instance = synthesize::<
                G1,
                G2,
                PedersenCommitment<GrumpkinCurveGroup>,
            >(c, &pp).unwrap().0;

            R1CSInstanceVar::new_variable(
                cs.clone(),
                || Ok(r1cs_instance),
                AllocationMode::Witness,
            ).unwrap()
        };

        let auxiliary_input_E_2 = {
            // the circuit
            let c = {
                let g1 = instance.E.clone();
                let g2 = acc.E.clone();
                // E = E_temp + (beta * (1- beta)) * Q
                let E_temp = (instance.E.clone() * beta_scalar) + (acc.E.clone() * (ScalarField::ONE - beta_scalar));
                let g_out = E_temp + Q * (beta_scalar * (ScalarField::ONE - beta_scalar)) ;
                SecondaryCircuit {
                    g1: Q,
                    g2: acc.E,
                    g_out,
                    r: convert_field_one_to_field_two::<ScalarField, BaseField>(beta_scalar * (ScalarField::ONE - beta_scalar)),
                    flag: true,
                }
            };

            let r1cs_instance = synthesize::<
                G1,
                G2,
                PedersenCommitment<GrumpkinCurveGroup>,
            >(c, &pp).unwrap().0;

            R1CSInstanceVar::new_variable(
                cs.clone(),
                || Ok(r1cs_instance),
                AllocationMode::Witness,
            ).unwrap()
        };

        let Q_var = NonNativeAffineVar::new_variable(
            ns!(cs, "Q_var"),
            || Ok(Q),
            AllocationMode::Witness,
        ).unwrap();


        let verifier = AccumulatorVerifierVar {
            auxiliary_input_C,
            auxiliary_input_T,
            auxiliary_input_E_1,
            auxiliary_input_E_2,
            beta: beta_var,
            beta_non_native,
            Q: Q_var,
            instance: instance_var,
            acc: acc_var,
            n: n as u32,
            m: m as u32,
        };

        println!("number of constraint for initialisation: {}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap())
    }
}



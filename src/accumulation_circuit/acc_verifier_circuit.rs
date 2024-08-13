use std::borrow::Borrow;
use std::marker::PhantomData;
use std::ops::Add;

use ark_ec::CurveConfig;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::{BigInteger64, Field, PrimeField};
use ark_r1cs_std::{R1CSVar, ToBitsGadget};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
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
use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;
use crate::gadgets::r1cs::R1CSInstance;
use crate::nova::commitment::CommitmentScheme;
use crate::nova::cycle_fold::coprocessor::{Circuit as SecondaryCircuit, synthesize};
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
    /// the randomness used for taking linear combination
    pub beta: FpVar<G1::ScalarField>,
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
            Q,
            instance,
            acc,
            n: 0,
            m: 0,
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
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
{
    pub fn accumulate(&self) {
        // Poseidon hash
        let beta = FpVar::new_variable(
            ns!(self.acc.cs(), "beta"),
            || Ok(G1::ScalarField::rand(&mut thread_rng())),
            AllocationMode::Input,
        ).unwrap();

        // Non-native scalar multiplication: linear combination of C
        {
            // parse inputs
            let (flag,
                r,
                g1,
                g2,
                g_out
            ) = self.auxiliary_input_C.parse_secondary_io().unwrap();
            // check equality tests
            self.instance.C_var.enforce_equal(&g1).expect("error while enforcing equality");
            self.acc.C_var.enforce_equal(&g2).expect("error while enforcing equality");
            // enforce correctness of the challenge
        }


        // Non-native scalar multiplication: linear combination of T
        //let _ = &self.acc.T_var + &self.instance.T_var.scalar_mul_le(beta_bits.iter()).unwrap();

        // Non-native scalar multiplication: linear combination of E
        //let _ = &self.acc.E_var + &self.instance.E_var.scalar_mul_le(beta_bits.iter()).unwrap();

        // Native field operation: linear combination of b
        let _ = &self.acc.b_var + &beta * &self.instance.b_var;

        // Native field operation: linear combination of c
        let _ = &self.acc.c_var + &beta * &self.instance.c_var;

        // Native field operation: linear combination of y
        let _ = &self.acc.y_var + &beta * &self.instance.y_var;

        // Native field operation: linear combination of z_b
        let _ = &self.acc.z_b_var + &beta * &self.instance.z_b_var;

        // Native field operation: linear combination of z_c
        let _ = &self.acc.z_c_var + &beta * &self.instance.z_c_var;

        // Native field operation: equality assertion that z_b = b^n-1 for the first instance
        let n = BigInteger64::from(self.n);
        let _ = self.acc.b_var.pow_by_constant(n.as_ref()).expect("TODO: panic message");

        // Native field operation: equality assertion that z_c = c^m-1 for the first instance
        let m = BigInteger64::from(self.m);
        let _ = self.acc.c_var.pow_by_constant(m.as_ref()).expect("TODO: panic message");
    }
}


#[cfg(test)]
mod tests {
    use std::fmt::Debug;
    use ark_ec::short_weierstrass::Projective;
    use ark_ff::Field;
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_relations::ns;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::accumulation_circuit::acc_instance_circuit::{AccumulatorInstanceCircuit, AccumulatorInstanceCircuitVar};
    use crate::accumulation_circuit::acc_verifier_circuit::{AccumulatorVerifier, AccumulatorVerifierVar};
    use crate::constant_for_curves::{G1, ScalarField, BaseField, G2};
    use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;
    use crate::hash::pederson::PedersenCommitment;
    use crate::nova::commitment::CommitmentScheme;
    use crate::nova::cycle_fold::coprocessor::{Circuit, setup_shape, synthesize};
    use crate::nova::cycle_fold::coprocessor_constraints::R1CSInstanceVar;
    use crate::utils::cast_field_element;

    /// TODO: change this to actual instance where I can write decider for it
    fn random_instance(n: u32, m: u32) -> AccumulatorInstanceCircuit::<G1> {
        AccumulatorInstanceCircuit::<G1> {
            C: Projective::rand(&mut thread_rng()),
            T: Projective::rand(&mut thread_rng()),
            E: Projective::rand(&mut thread_rng()),
            b: ScalarField::rand(&mut thread_rng()),
            c: ScalarField::rand(&mut thread_rng()),
            y: ScalarField::rand(&mut thread_rng()),
            z_b: ScalarField::rand(&mut thread_rng()),
            z_c: ScalarField::rand(&mut thread_rng()),
        }
    }



    #[test]
    fn initialisation_test() {
        // build an instance of AccInstanceCircuit
        let instance = random_instance(1000u32, 1000u32);

        let acc = random_instance(1000u32, 1000u32);

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

        let r = ScalarField::from(2u128);
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
                g1: instance.C.clone(),
                g2: acc.C,
                g_out,
                r: unsafe { cast_field_element::<ScalarField, BaseField>(&r) },
                flag: true,
            }
        };

        let auxiliary_input = {
            let shape = setup_shape::<G1, G2>().unwrap();
            let pp = PedersenCommitment::<ark_grumpkin::Projective>::setup(shape.num_vars, b"test", &());
            let (u, _) = synthesize::<
                G1,
                G2,
                PedersenCommitment<ark_grumpkin::Projective>,
            >(c, &pp).unwrap();

            R1CSInstanceVar::new_variable(
                cs.clone(),
                || Ok(u.clone()),
                AllocationMode::Constant,
            ).unwrap()
        };

        let Q = NonNativeAffineVar::new_variable(
            ns!(cs, "C"),
            || Ok(instance.C),
            AllocationMode::Witness,
        ).unwrap();


        let verifier = AccumulatorVerifierVar {
            auxiliary_input_C: auxiliary_input.clone(),
            auxiliary_input_T: auxiliary_input.clone(),
            auxiliary_input_E_1: auxiliary_input.clone(),
            auxiliary_input_E_2: auxiliary_input.clone(),
            beta: r_var,
            Q,
            instance: instance_var,
            acc: acc_var,
            n: 0,
            m: 0,
        };

        println!("number of constraint for initialisation: {}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap())
    }

    /*
    #[test]
    fn initialisation_test() {
        // build an instance of AccInstanceCircuit
        let instance = random_instance(1000u32, 1000u32);

        let acc = random_instance(1000u32, 1000u32);


        // a constraint system
        let cs = ConstraintSystem::<Fr>::new_ref();

        // randomness
        let beta = Fr::from(2u128);

        // the shape of the R1CS instance
        let shape = setup_shape::<Config, GrumpkinConfig>().unwrap();

        // public parameters of Pedersen
        let pp = PedersenCommitment::<ark_vesta::Projective>::setup(shape.num_vars, b"test", &());

        let auxiliary_input_C = {
            // the circuit
            let c = {
                let g1 = instance.C.clone();
                let g2 = acc.C.clone();
                // C'' = (1 - beta) * acc.C + beta * instance.C
                let g_out = (g1 * beta) + (g2 * (Fr::ONE - beta));
                Circuit {
                    g1: instance.C,
                    g2: acc.C,
                    g_out,
                    r: unsafe { cast_field_element::<Fr, Fq>(&beta) },
                    flag: false,
                }
            };

            synthesize::<
                Config,
                GrumpkinConfig,
                PedersenCommitment<ark_vesta::Projective>,
            >(c, &pp).unwrap().0
        };

        let auxiliary_input_T = {
            // the circuit
            let c = {
                let g1 = instance.T.clone();
                let g2 = acc.T.clone();
                // C'' = (1 - beta) * acc.T + beta * instance.T
                let g_out = (g1 * beta) + (g2 * (Fr::ONE - beta));
                Circuit {
                    g1: instance.T,
                    g2: acc.T,
                    g_out,
                    r: unsafe { cast_field_element::<Fr, Fq>(&beta) },
                    flag: false,
                }
            };

            synthesize::<
                Config,
                GrumpkinConfig,
                PedersenCommitment<ark_vesta::Projective>,
            >(c, &pp).unwrap().0
        };

        let auxiliary_input_E_1 = {
            // the circuit
            let c = {
                let g1 = instance.E.clone();
                let g2 = acc.E.clone();
                // C'' = (1 - beta) * acc.T + beta * instance.T
                let g_out = (g1 * beta) + (g2 * (Fr::ONE - beta));
                Circuit {
                    g1: instance.E,
                    g2: acc.E,
                    g_out,
                    r: unsafe { cast_field_element::<Fr, Fq>(&beta) },
                    flag: false,
                }
            };

            synthesize::<
                Config,
                GrumpkinConfig,
                PedersenCommitment<ark_vesta::Projective>,
            >(c, &pp).unwrap().0
        };


        let accumulator_verifier = AccumulatorVerifier {
            beta,
            auxiliary_input_C: R1CSInstance {},
            auxiliary_input_T: R1CSInstance {},
            auxiliary_input_E_1: R1CSInstance {},
            auxiliary_input_E_2: R1CSInstance {},
            instance: AccumulatorInstance {},
            acc: AccumulatorInstance {},
        };
    }
     */
}



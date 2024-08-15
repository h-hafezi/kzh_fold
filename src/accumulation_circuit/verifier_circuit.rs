use std::borrow::Borrow;
use std::marker::PhantomData;
use std::ops::Add;

use ark_ec::CurveConfig;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
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

use crate::accumulation::accumulator::AccInstance;
use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use crate::constant_for_curves::{BaseField, ScalarField};
use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;
use crate::gadgets::non_native::util::{convert_field_one_to_field_two, non_native_to_fpvar};
use crate::gadgets::r1cs::{R1CSInstance, RelaxedR1CSInstance};
use crate::nova::commitment::CommitmentScheme;
use crate::nova::cycle_fold::coprocessor::{SecondaryCircuit as SecondaryCircuit, synthesize};
use crate::nova::cycle_fold::coprocessor_constraints::{R1CSInstanceVar, RelaxedR1CSInstanceVar};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccumulatorVerifier<G1, G2, C2, E>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField>,
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
    pub acc_instance: AccInstance<E>,
    /// the running accumulator
    pub acc_running: AccInstance<E>,
    /// the result accumulator
    pub acc_result: AccInstance<E>,

    /// running cycle fold instance
    pub cycle_fold_running_instance: RelaxedR1CSInstance<G2, C2>,

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
    pub auxiliary_input_C_var: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have T'' = (1-beta) * T + beta * T' without scalar multiplication
    pub auxiliary_input_T_var: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E_{temp} = (1-beta) * E + beta * E' without scalar multiplication
    pub auxiliary_input_E_1_var: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E'' = E_{temp} + beta * (1-beta) * Q without scalar multiplication
    pub auxiliary_input_E_2_var: R1CSInstanceVar<G2, C2>,

    /// the randomness used for taking linear combination and its non-native counterpart
    pub beta_var: FpVar<G1::ScalarField>,
    pub beta_var_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// accumulation proof
    pub Q_var: NonNativeAffineVar<G1>,

    pub acc_instance_var: AccumulatorInstanceVar<G1>,
    pub acc_running_var: AccumulatorInstanceVar<G1>,
    pub acc_result_var: AccumulatorInstanceVar<G1>,

    pub cycle_fold_running_instance_var: RelaxedR1CSInstanceVar<G2, C2>,

    // these are constant values
    pub n: u32,
    pub m: u32,
}

impl<G1, G2, C2, E> AllocVar<AccumulatorVerifier<G1, G2, C2, E>, G1::ScalarField> for AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField>,
{
    fn new_variable<T: Borrow<AccumulatorVerifier<G1, G2, C2, E>>>(cs: impl Into<Namespace<G1::ScalarField>>, f: impl FnOnce() -> Result<T, SynthesisError>, mode: AllocationMode) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        let auxiliary_input_C_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_C"),
            || Ok(circuit.map(|e| e.auxiliary_input_C.clone()).unwrap()),
            mode,
        ).unwrap();

        let auxiliary_input_T_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_T"),
            || Ok(circuit.map(|e| e.auxiliary_input_T.clone()).unwrap()),
            mode,
        ).unwrap();

        let auxiliary_input_E_1_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_1"),
            || Ok(circuit.map(|e| e.auxiliary_input_E_1.clone()).unwrap()),
            mode,
        ).unwrap();

        let auxiliary_input_E_2_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_2"),
            || Ok(circuit.map(|e| e.auxiliary_input_E_2.clone()).unwrap()),
            mode,
        ).unwrap();


        let acc_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "instance"),
            || circuit.map(|e| e.acc_instance.clone()),
            mode,
        ).unwrap();

        let acc_running_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "acc"),
            || circuit.map(|e| e.acc_running.clone()),
            mode,
        ).unwrap();

        let acc_result_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "result acc"),
            || circuit.map(|e| e.acc_result.clone()),
            mode,
        ).unwrap();


        let cycle_fold_running_instance_var = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.cycle_fold_running_instance.clone()),
            mode,
        ).unwrap();


        let beta_var = FpVar::new_variable(
            ns!(cs, "beta"),
            || circuit.map(|e| e.beta.clone()),
            mode,
        ).unwrap();

        let beta_var_non_native = NonNativeFieldVar::new_variable(
            ns!(cs, "non native beta"),
            || circuit.map(|e| convert_field_one_to_field_two::<G1::ScalarField, G1::BaseField>(e.beta.clone())),
            mode,
        ).unwrap();


        let Q_var = NonNativeAffineVar::new_variable(
            ns!(cs, "Q"),
            || circuit.map(|e| e.Q),
            mode,
        ).unwrap();


        Ok(AccumulatorVerifierVar {
            auxiliary_input_C_var,
            auxiliary_input_T_var,
            auxiliary_input_E_1_var,
            auxiliary_input_E_2_var,
            beta_var,
            beta_var_non_native,
            Q_var,
            acc_instance_var,
            acc_running_var,
            acc_result_var,
            cycle_fold_running_instance_var,
            n: circuit.map(|e| e.n).unwrap(),
            m: circuit.map(|e| e.m).unwrap(),
        })
    }
}

/// Here we assume instance to be A.X and acc to be A.X' ==> beta * instance + (1-beta) * instance
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
        let beta_ = non_native_to_fpvar(&self.beta_var_non_native);
        self.beta_var.enforce_equal(&beta_).expect("error while enforcing equality");

        // todo: compute Poseidon hash and make sure it's consistent with input beta


        // Non-native scalar multiplication: linear combination of C
        let (flag,
            r,
            g1,
            g2,
            C_var
        ) = self.auxiliary_input_C_var.parse_secondary_io().unwrap();
        // g1 == acc.C
        self.acc_running_var.C_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.C
        self.acc_instance_var.C_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).expect("error while enforcing equality");
        // check out the result C_var is consistent with result_acc
        C_var.enforce_equal(&self.acc_result_var.C_var).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination of T
        let (flag,
            r,
            g1,
            g2,
            T_var
        ) = self.auxiliary_input_T_var.parse_secondary_io().unwrap();
        // g1 == acc.T
        self.acc_running_var.T_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.C
        self.acc_instance_var.T_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).expect("error while enforcing equality");
        // check out the result T_var is consistent with result_acc
        T_var.enforce_equal(&self.acc_result_var.T_var).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination E_temp = (instance.E * (1-beta) + acc.E * beta)
        let (flag,
            r,
            g1,
            g2,
            E_temp
        ) = self.auxiliary_input_E_1_var.parse_secondary_io().unwrap();
        // g1 == acc.E
        self.acc_running_var.E_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.E
        self.acc_instance_var.E_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination E'' = E_{temp} + (1-beta) * beta * Q
        let (flag,
            r,
            g1,
            g2,
            E_var
        ) = self.auxiliary_input_E_2_var.parse_secondary_io().unwrap();
        // g1 == Q
        g1.enforce_equal(&self.Q_var).expect("error while enforcing equality");
        // g2 == E_temp
        g2.enforce_equal(&E_temp).expect("error while enforcing equality");
        // enforce flag to be true
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        // check r to be equal to beta
        let beta_times_beta_minus_one = self.beta_var_non_native.clone() * (NonNativeFieldVar::one() - self.beta_var_non_native.clone());
        r.enforce_equal(&beta_times_beta_minus_one).expect("error while enforcing equality");
        // check out the result E_var is consistent with result_acc
        E_var.enforce_equal(&self.acc_result_var.E_var).expect("error while enforcing equality");


        let beta_minus_one = FpVar::<G1::ScalarField>::one() - &self.beta_var;

        // Native field operation: linear combination of b
        let b_var = &self.beta_var * &self.acc_running_var.b_var + &beta_minus_one * &self.acc_instance_var.b_var;
        // check out the result b_var is consistent with result_acc
        b_var.enforce_equal(&self.acc_result_var.b_var).expect("error while enforcing equality");


        // Native field operation: linear combination of c
        let c_var = &self.beta_var * &self.acc_running_var.c_var + &beta_minus_one * &self.acc_instance_var.c_var;
        // check out the result c_var is consistent with result_acc
        c_var.enforce_equal(&self.acc_result_var.c_var).expect("error while enforcing equality");


        // Native field operation: linear combination of y
        let y_var = &self.beta_var * &self.acc_running_var.y_var + &beta_minus_one * &self.acc_instance_var.y_var;
        // check out the result y_var is consistent with result_acc
        y_var.enforce_equal(&self.acc_result_var.y_var).expect("error while enforcing equality");


        // Native field operation: linear combination of z_b
        let z_b_var = &self.beta_var * &self.acc_running_var.z_b_var + &beta_minus_one * &self.acc_instance_var.z_b_var;
        // check out the result z_b_var is consistent with result_acc
        z_b_var.enforce_equal(&self.acc_result_var.z_b_var).expect("error while enforcing equality");


        // Native field operation: linear combination of z_c
        let z_c_var = &self.beta_var * &self.acc_running_var.z_c_var + &beta_minus_one * &self.acc_instance_var.z_c_var;
        // check out the result z_c_var is consistent with result_acc
        z_c_var.enforce_equal(&self.acc_result_var.z_c_var).expect("error while enforcing equality");


        // Native field operation: equality assertion that z_b = b^n-1 for the first instance
        let n = BigInteger64::from(self.n);
        let z_b_ = self.acc_running_var.b_var.pow_by_constant(n.as_ref()).expect("error while enforcing equality");


        // Native field operation: equality assertion that z_c = c^m-1 for the first instance
        let m = BigInteger64::from(self.m);
        let z_c_ = self.acc_running_var.c_var.pow_by_constant(m.as_ref()).expect("error while enforcing equality");


        // Conditional check: if instance.E_var == 0, then enforce z_b_ == instance.z_b_var and z_c_ == instance.z_c_var
        let is_E_zero = &self.acc_instance_var.E_var.infinity.clone();
        z_b_.conditional_enforce_equal(&self.acc_instance_var.z_b_var, is_E_zero).expect("error while enforcing z_b equality under condition");
        z_c_.conditional_enforce_equal(&self.acc_instance_var.z_c_var, is_E_zero).expect("error while enforcing z_c equality under condition");
    }
}


#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use ark_ec::short_weierstrass::Projective;
    use ark_ff::{Field, Fp};
    use ark_grumpkin::GrumpkinConfig;
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
    use ark_relations::ns;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};

    use crate::accumulation::accumulator::{Accumulator, AccumulatorTrait};
    use crate::accumulation::accumulator::tests::{get_satisfying_accumulator, get_srs};
    use crate::accumulation_circuit::affine_to_projective;
    use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
    use crate::accumulation_circuit::verifier_circuit::AccumulatorVerifierVar;
    use crate::constant_for_curves::{BaseField, G1, G2, ScalarField};
    use crate::gadgets::non_native::short_weierstrass::NonNativeAffineVar;
    use crate::gadgets::non_native::util::convert_field_one_to_field_two;
    use crate::gadgets::r1cs::{R1CSInstance, R1CSWitness, RelaxedR1CSInstance};
    use crate::gadgets::r1cs::r1cs::{commit_T, RelaxedR1CSWitness};
    use crate::hash::pederson::PedersenCommitment;
    use crate::nova::commitment::CommitmentScheme;
    use crate::nova::cycle_fold::coprocessor::{SecondaryCircuit, setup_shape, synthesize};
    use crate::nova::cycle_fold::coprocessor_constraints::{R1CSInstanceVar, RelaxedR1CSInstanceVar};

    type GrumpkinCurveGroup = ark_grumpkin::Projective;
    type C2 = PedersenCommitment<GrumpkinCurveGroup>;

    pub fn secondary_circuit_to_r1cs(
        pp: &Vec<<Projective<GrumpkinConfig> as ark_ec::ScalarMul>::MulBase>,
        cs: ConstraintSystemRef<ScalarField>,
        c: SecondaryCircuit<G1>,
    ) -> (R1CSInstance<G2, C2>, R1CSWitness<G2>, R1CSInstanceVar<G2, C2>) {
        let (r1cs_instance, r1cs_witness) = synthesize::<
            G1,
            G2,
            PedersenCommitment<GrumpkinCurveGroup>,
        >(c, pp).unwrap();

        (r1cs_instance.clone(),
         r1cs_witness,
         R1CSInstanceVar::new_variable(
             cs.clone(),
             || Ok(r1cs_instance),
             AllocationMode::Witness,
         ).unwrap())
    }

    pub fn randomness_different_formats(cs: ConstraintSystemRef<ScalarField>,
                                        beta: ScalarField) -> (
        BaseField,
        FpVar<ScalarField>,
        NonNativeFieldVar<BaseField, ScalarField>
    ) {
        let beta_base = convert_field_one_to_field_two::<ScalarField, BaseField>(beta);
        let beta_var = FpVar::new_variable(
            cs.clone(),
            || Ok(beta.clone()),
            AllocationMode::Witness,
        ).unwrap();
        let beta_var_non_native = NonNativeFieldVar::new_variable(
            cs.clone(),
            || Ok(beta_base.clone()),
            AllocationMode::Witness,
        ).unwrap();
        (beta_base, beta_var, beta_var_non_native)
    }

    #[test]
    fn initialisation_test() {
        // specifying degrees of polynomials
        let n = 16;
        let m = 16;

        // get a random srs
        let srs = get_srs(n, m);

        // a constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // build an instance of AccInstanceCircuit
        let acc_instance = get_satisfying_accumulator(&srs);
        let acc_running = get_satisfying_accumulator(&srs);

        // the randomness in different formats
        let beta_scalar = ScalarField::from(2u128);
        let (beta_base, beta_var, beta_var_non_native) = randomness_different_formats(cs.clone(), beta_scalar);


        // accumulate proof
        let (result_instance, _, Q) = Accumulator::prove(&srs, &beta_scalar, &acc_running, &acc_instance);
        let Q = Projective::new(Q.x, Q.y, BaseField::ONE);


        // make a circuit_var
        let instance_var = AccumulatorInstanceVar::new_variable(
            cs.clone(),
            || Ok(acc_instance.instance.clone()),
            AllocationMode::Witness,
        ).unwrap();

        let acc_var = AccumulatorInstanceVar::new_variable(
            cs.clone(),
            || Ok(acc_running.instance.clone()),
            AllocationMode::Witness,
        ).unwrap();

        let result_acc_circuit_var = AccumulatorInstanceVar::new_variable(
            cs.clone(),
            || Ok(result_instance.clone()),
            AllocationMode::Witness,
        ).unwrap();


        // the shape of the R1CS instance
        let shape = setup_shape::<G1, G2>().unwrap();

        // public parameters of Pedersen
        let pp = PedersenCommitment::<GrumpkinCurveGroup>::setup(shape.num_vars, b"test", &());


        let (auxiliary_input_C, auxiliary_input_C_witness, auxiliary_input_C_var) = {
            secondary_circuit_to_r1cs(&pp, cs.clone(), {
                let g1 = affine_to_projective(acc_running.instance.C.clone());
                let g2 = affine_to_projective(acc_instance.instance.C.clone());
                // C'' = beta * acc.C + (1 - beta) * instance.C
                let g_out = (g1 * beta_scalar) + (g2 * (ScalarField::ONE - beta_scalar));
                SecondaryCircuit {
                    g1,
                    g2,
                    g_out,
                    r: beta_base,
                    flag: false,
                }
            })
        };

        let (auxiliary_input_T, auxiliary_input_T_witness, auxiliary_input_T_var) = {
            secondary_circuit_to_r1cs(&pp, cs.clone(), {
                let g1 = affine_to_projective(acc_running.instance.T.clone());
                let g2 = affine_to_projective(acc_instance.instance.T.clone());
                // T'' = beta * acc.T + (1 - beta) * instance.T
                let g_out = (g1 * beta_scalar) + (g2 * (ScalarField::ONE - beta_scalar));
                SecondaryCircuit {
                    g1,
                    g2,
                    g_out,
                    r: beta_base,
                    flag: false,
                }
            })
        };

        let (auxiliary_input_E_1, auxiliary_input_E_1_witness, auxiliary_input_E_1_var) = {
            secondary_circuit_to_r1cs(&pp, cs.clone(), {
                let g1 = affine_to_projective(acc_instance.instance.E.clone());
                let g2 = affine_to_projective(acc_running.instance.E.clone());
                // E_temp = beta * acc.E + (1 - beta) * instance.E
                let g_out = (g1 * beta_scalar) + (g2 * (ScalarField::ONE - beta_scalar));
                SecondaryCircuit {
                    g1,
                    g2,
                    g_out,
                    r: beta_base,
                    flag: false,
                }
            })
        };

        let (auxiliary_input_E_2, auxiliary_input_E_2_witness, auxiliary_input_E_2_var) = {
            // the circuit
            secondary_circuit_to_r1cs(&pp, cs.clone(), {
                let g1 = affine_to_projective(acc_instance.instance.E.clone());
                let g2 = affine_to_projective(acc_running.instance.E.clone());
                // E = E_temp + (beta * (1- beta)) * Q
                let E_temp = (g1 * beta_scalar) + (g2 * (ScalarField::ONE - beta_scalar));
                let g_out = E_temp + Q * (beta_scalar * (ScalarField::ONE - beta_scalar));
                SecondaryCircuit {
                    g1,
                    g2,
                    g_out,
                    r: convert_field_one_to_field_two::<ScalarField, BaseField>(beta_scalar * (ScalarField::ONE - beta_scalar)),
                    flag: true,
                }
            })
        };


        let Q_var = NonNativeAffineVar::new_variable(
            ns!(cs, "Q var"),
            || Ok(Q),
            AllocationMode::Witness,
        ).unwrap();

        // make the trivial running instance
        let trivial_instance = RelaxedR1CSInstance::new(&shape);
        let trivial_witness = RelaxedR1CSWitness::zero(&shape);
        shape.is_relaxed_satisfied(&trivial_instance, &trivial_witness, &pp).expect("panic!");


        let trivial_cycle_fold_running_instance = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "trivial cycle fold running instance"),
            || Ok(trivial_instance),
            AllocationMode::Witness,
        ).unwrap();


        let verifier = AccumulatorVerifierVar {
            auxiliary_input_C_var,
            auxiliary_input_T_var,
            auxiliary_input_E_1_var,
            auxiliary_input_E_2_var,
            beta_var,
            beta_var_non_native,
            Q_var,
            acc_instance_var: instance_var,
            acc_running_var: acc_var,
            acc_result_var: result_acc_circuit_var,
            cycle_fold_running_instance_var: trivial_cycle_fold_running_instance,
            n: n as u32,
            m: m as u32,
        };

        println!("number of constraint for initialisation: {}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap())
    }
}


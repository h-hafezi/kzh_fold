use crate::commitment;
use std::borrow::Borrow;
use std::marker::PhantomData;
use std::ops::Add;

use ark_crypto_primitives::sponge::Absorb;
use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
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
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use ark_std::UniformRand;
use rand::thread_rng;

use crate::accumulation::accumulator::{AccInstance, AccSRS};
use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use crate::accumulation_circuit::prover::AccumulatorVerifierCircuitProver;
use crate::accumulation_circuit::randomness_different_formats;
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::gadgets::non_native::util::{convert_field_one_to_field_two, non_native_to_fpvar};
use crate::gadgets::r1cs::{R1CSInstance, RelaxedR1CSInstance};
use crate::hash::poseidon::{PoseidonHashVar, PoseidonHashVarTrait};
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

    /// accumulation proof for accumulators
    pub Q: Projective<G1>,

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub com_C: Projective<G2>,
    pub com_T: Projective<G2>,
    pub com_E_1: Projective<G2>,
    pub com_E_2: Projective<G2>,

    /// the instance to be folded
    pub current_accumulator_instance: AccInstance<E>,
    /// the running accumulator
    pub running_accumulator_instance: AccInstance<E>,
    /// the result accumulator
    pub final_accumulator_instance: AccInstance<E>,

    /// running cycle fold instance
    pub running_cycle_fold_instance: RelaxedR1CSInstance<G2, C2>,
    pub final_cycle_fold_instance: RelaxedR1CSInstance<G2, C2>,

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

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub com_C_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub com_T_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub com_E_1_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub com_E_2_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,

    pub current_accumulator_instance_var: AccumulatorInstanceVar<G1>,
    pub running_accumulator_instance_var: AccumulatorInstanceVar<G1>,
    pub final_accumulator_instance_var: AccumulatorInstanceVar<G1>,

    pub running_cycle_fold_instance_var: RelaxedR1CSInstanceVar<G2, C2>,
    pub final_cycle_fold_instance_var: RelaxedR1CSInstanceVar<G2, C2>,

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

        // auxiliary inputs
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


        // accumulator instances
        let current_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "instance"),
            || circuit.map(|e| e.current_accumulator_instance.clone()),
            mode,
        ).unwrap();

        let running_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "acc"),
            || circuit.map(|e| e.running_accumulator_instance.clone()),
            mode,
        ).unwrap();

        let final_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "result acc"),
            || circuit.map(|e| e.final_accumulator_instance.clone()),
            mode,
        ).unwrap();


        // cycle fold instances
        let running_cycle_fold_instance_var = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.running_cycle_fold_instance.clone()),
            mode,
        ).unwrap();

        let final_cycle_fold_instance_var = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.final_cycle_fold_instance.clone()),
            mode,
        ).unwrap();


        // randomness variables
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


        // folding proofs for cycle fold and accumulator
        let Q_var = NonNativeAffineVar::new_variable(
            ns!(cs, "Q"),
            || circuit.map(|e| e.Q),
            mode,
        ).unwrap();

        let com_C_var = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.com_C.clone()),
            mode,
        ).unwrap();

        let com_T_var = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.com_T.clone()),
            mode,
        ).unwrap();

        let com_E_1_var = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.com_E_1.clone()),
            mode,
        ).unwrap();

        let com_E_2_var = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.com_E_2.clone()),
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
            com_C_var,
            com_T_var,
            com_E_1_var,
            com_E_2_var,
            current_accumulator_instance_var,
            running_accumulator_instance_var,
            final_accumulator_instance_var,
            running_cycle_fold_instance_var,
            final_cycle_fold_instance_var,
            n: circuit.map(|e| e.n).unwrap(),
            m: circuit.map(|e| e.m).unwrap(),
        })
    }
}

/// Here we assume current_acc to be A.X and running_acc to be A.X' ==> beta * running_acc + (1-beta) * current_acc
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
    pub fn accumulate(&self)
    where
        <G2 as CurveConfig>::BaseField: Absorb,
    {
        // checking beta and non_native beta are consistent
        let beta_bits = self.beta_var_non_native.to_bits_le().unwrap();
        let beta_ = Boolean::le_bits_to_fp_var(beta_bits.as_slice()).unwrap();
        self.beta_var.enforce_equal(&beta_).expect("error while enforcing equality");

        // compute Poseidon hash and make sure it's consistent with input beta
        let mut hash_object = PoseidonHashVar::new(self.current_accumulator_instance_var.cs());
        let mut sponge = Vec::new();
        sponge.extend(self.current_accumulator_instance_var.to_sponge_field_elements().unwrap());
        sponge.extend(self.running_accumulator_instance_var.to_sponge_field_elements().unwrap());
        sponge.extend(self.Q_var.to_sponge_field_elements().unwrap());
        hash_object.update_sponge(sponge);
        hash_object.output().enforce_equal(&self.beta_var).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination of C
        let (flag,
            r,
            g1,
            g2,
            C_var
        ) = self.auxiliary_input_C_var.parse_secondary_io::<G1>().unwrap();
        // g1 == acc.C
        self.running_accumulator_instance_var.C_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.C
        self.current_accumulator_instance_var.C_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).expect("error while enforcing equality");
        // check out the result C_var is consistent with result_acc
        C_var.enforce_equal(&self.final_accumulator_instance_var.C_var).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination of T
        let (flag,
            r,
            g1,
            g2,
            T_var
        ) = self.auxiliary_input_T_var.parse_secondary_io::<G1>().unwrap();
        // g1 == acc.T
        self.running_accumulator_instance_var.T_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.C
        self.current_accumulator_instance_var.T_var.enforce_equal(&g2).expect("error while enforcing equality");
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).expect("error while enforcing equality");
        // check out the result T_var is consistent with result_acc
        T_var.enforce_equal(&self.final_accumulator_instance_var.T_var).expect("error while enforcing equality");


        // Non-native scalar multiplication: linear combination E_temp = (instance.E * (1-beta) + acc.E * beta)
        let (flag,
            r,
            g1,
            g2,
            E_temp
        ) = self.auxiliary_input_E_1_var.parse_secondary_io::<G1>().unwrap();
        // g1 == acc.E
        self.running_accumulator_instance_var.E_var.enforce_equal(&g1).expect("error while enforcing equality");
        // g2 == instance.E
        self.current_accumulator_instance_var.E_var.enforce_equal(&g2).expect("error while enforcing equality");
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
        ) = self.auxiliary_input_E_2_var.parse_secondary_io::<G1>().unwrap();
        // g1 == Q
        g1.enforce_equal(&self.Q_var).expect("error while enforcing equality");
        // g2 == E_temp
        g2.enforce_equal(&E_temp).expect("error while enforcing equality");
        // enforce flag to be true
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        // check r to be equal to beta
        let beta_times_beta_minus_one = self.beta_var_non_native.clone() - self.beta_var_non_native.square().unwrap();
        //r.enforce_equal(&beta_times_beta_minus_one).expect("error while enforcing equality");
        // check out the result E_var is consistent with result_acc
        E_var.enforce_equal(&self.final_accumulator_instance_var.E_var).expect("error while enforcing equality");


        let beta_minus_one = FpVar::<G1::ScalarField>::one() - &self.beta_var;

        // Native field operation: linear combination of x
        for i in 0..self.running_accumulator_instance_var.x_var.len() {
            let x_var = &self.beta_var * &self.running_accumulator_instance_var.x_var[i] + &beta_minus_one * &self.current_accumulator_instance_var.x_var[i];
            // check out the result b_var is consistent with result_acc
            x_var.enforce_equal(&self.final_accumulator_instance_var.x_var[i]).expect("error while enforcing equality");
        }

        // Native field operation: linear combination of x
        for i in 0..self.running_accumulator_instance_var.y_var.len() {
            let y_var = &self.beta_var * &self.running_accumulator_instance_var.y_var[i] + &beta_minus_one * &self.current_accumulator_instance_var.y_var[i];
            // check out the result b_var is consistent with result_acc
            y_var.enforce_equal(&self.final_accumulator_instance_var.y_var[i]).expect("error while enforcing equality");
        }

        // Native field operation: linear combination of z_c
        let z_var = &self.beta_var * &self.running_accumulator_instance_var.z_var + &beta_minus_one * &self.current_accumulator_instance_var.z_var;
        // check out the result z_c_var is consistent with result_acc
        z_var.enforce_equal(&self.final_accumulator_instance_var.z_var).expect("error while enforcing equality");

        let final_instance = self.running_cycle_fold_instance_var.fold(
            &[((&self.auxiliary_input_C_var, None), &self.com_C_var, &self.beta_var_non_native, &beta_bits),
                ((&self.auxiliary_input_T_var, None), &self.com_T_var, &self.beta_var_non_native, &beta_bits),
                ((&self.auxiliary_input_E_1_var, None), &self.com_E_1_var, &self.beta_var_non_native, &beta_bits),
                ((&self.auxiliary_input_E_2_var, None), &self.com_E_2_var, &self.beta_var_non_native, &beta_bits),
            ]
        ).unwrap();

        self.final_cycle_fold_instance_var.X.enforce_equal(&final_instance.X).expect("XXX: panic message");
        self.final_cycle_fold_instance_var.commitment_E.enforce_equal(&final_instance.commitment_E).expect("XXX: panic message");
        self.final_cycle_fold_instance_var.commitment_W.enforce_equal(&final_instance.commitment_W).expect("XXX: panic message");
    }
}

impl<G1, G2, C2> AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP = Vec<Affine<G2>>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    ProjectiveVar<G2, FpVar<<G2 as CurveConfig>::BaseField>>: AllocVar<<C2 as CommitmentScheme<Projective<G2>>>::Commitment, <G2 as CurveConfig>::BaseField>
{
    pub fn rand<E: Pairing>(srs: &AccSRS<E>, cs: ConstraintSystemRef<G1::ScalarField>) -> AccumulatorVerifierVar<G1, G2, C2>
    where
        E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField, BaseField=<G1 as CurveConfig>::BaseField>,
        <G2 as CurveConfig>::BaseField: Absorb,
        <G2 as CurveConfig>::ScalarField: Absorb
    {
        // get the prover
        let prover = AccumulatorVerifierCircuitProver::rand(&srs);

        // the randomness in different formats
        let beta_scalar = prover.beta.clone();
        let (_, beta_var, beta_var_non_native) = randomness_different_formats::<E>(cs.clone(), beta_scalar);

        // initialise accumulator variables
        let current_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "current accumulator instance var"),
            || Ok(prover.get_current_acc_instance().clone()),
            AllocationMode::Witness,
        ).unwrap();

        let running_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "running accumulator instance var"),
            || Ok(prover.get_running_acc_instance().clone()),
            AllocationMode::Witness,
        ).unwrap();

        let final_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "final accumulator instance var"),
            || Ok(prover.compute_result_accumulator_instance()),
            AllocationMode::Witness,
        ).unwrap();


        // initialise auxiliary input variables
        let auxiliary_input_C_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary input C var"),
            || Ok(prover.compute_auxiliary_input_C().0),
            AllocationMode::Witness,
        ).unwrap();

        let auxiliary_input_T_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary input T var"),
            || Ok(prover.compute_auxiliary_input_T().0),
            AllocationMode::Witness,
        ).unwrap();

        let auxiliary_input_E_1_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary input E_1 var"),
            || Ok(prover.compute_auxiliary_input_E_1().0),
            AllocationMode::Witness,
        ).unwrap();

        let auxiliary_input_E_2_var = R1CSInstanceVar::new_variable(
            ns!(cs, "auxiliary input E_2 var"),
            || Ok(prover.compute_auxiliary_input_E_2().0),
            AllocationMode::Witness,
        ).unwrap();


        // initialise Q variables
        let Q_var = NonNativeAffineVar::new_variable(
            ns!(cs, "Q var"),
            || Ok(prover.compute_proof_Q()),
            AllocationMode::Witness,
        ).unwrap();

        let cycle_fold_proof = prover.compute_cycle_fold_proofs_and_final_instance();

        let com_C_var = ProjectiveVar::new_variable(
            ns!(cs, "com_C_var"),
            || Ok(cycle_fold_proof.0),
            AllocationMode::Witness,
        ).unwrap();

        let com_T_var = ProjectiveVar::new_variable(
            ns!(cs, "com_T_var"),
            || Ok(cycle_fold_proof.1),
            AllocationMode::Witness,
        ).unwrap();

        let com_E_1_var = ProjectiveVar::new_variable(
            ns!(cs, "com_E_1_var"),
            || Ok(cycle_fold_proof.2),
            AllocationMode::Witness,
        ).unwrap();

        let com_E_2_var = ProjectiveVar::new_variable(
            ns!(cs, "com_E_2_var"),
            || Ok(cycle_fold_proof.3),
            AllocationMode::Witness,
        ).unwrap();


        // initialise cycle fold running instance var
        let running_cycle_fold_instance_var = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "running cycle fold instance var"),
            || Ok(prover.running_cycle_fold_instance),
            AllocationMode::Witness,
        ).unwrap();

        // initialise cycle fold running instance var
        let final_cycle_fold_instance_var = RelaxedR1CSInstanceVar::new_variable(
            ns!(cs, "final cycle fold instance var"),
            || Ok(cycle_fold_proof.4),
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
            com_C_var,
            com_T_var,
            com_E_1_var,
            com_E_2_var,
            current_accumulator_instance_var,
            running_accumulator_instance_var,
            final_accumulator_instance_var,
            running_cycle_fold_instance_var,
            final_cycle_fold_instance_var,
            n: prover.n,
            m: prover.m,
        };

        verifier
    }
}


#[cfg(test)]
pub mod tests {
    use std::fmt::Debug;
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;

    use ark_ec::short_weierstrass::Projective;
    use ark_ff::PrimeField;
    use ark_relations::r1cs::ConstraintSystem;
    use rand::thread_rng;

    use crate::accumulation::accumulator::Accumulator;
    use crate::accumulation_circuit::verifier_circuit::AccumulatorVerifierVar;
    use crate::commitment::CommitmentScheme;
    use crate::constant_for_curves::{E, G1, G2, ScalarField};
    use crate::hash::pederson::PedersenCommitment;
    use crate::pcs::multilinear_pcs::{PolyCommit, SRS};
    use crate::pcs::multilinear_pcs::PolyCommitTrait;

    type C2 = PedersenCommitment<Projective<G2>>;

    #[test]
    fn initialisation_test() {
        // specifying degrees of polynomials
        let n = 4;
        let m = 4;

        // get a random srs
        let srs = {
            let srs_pcs: SRS<E> = PolyCommit::<E>::setup(n, m, &mut thread_rng());
            Accumulator::setup(srs_pcs.clone(), &mut thread_rng())
        };

        // a constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        let verifier = AccumulatorVerifierVar::<G1, G2, C2>::rand(&srs, cs.clone());

        println!("number of constraint for initialisation: {}", cs.num_constraints());
        verifier.accumulate();
        println!("number of constraint for initialisation: {}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap());

        // write the witness into a file
        let cs_borrow = cs.borrow().unwrap();
        let witness = cs_borrow.witness_assignment.clone();
        let pub_io = cs_borrow.instance_assignment.clone();

        let file = File::create("witness.txt").unwrap();
        let mut writer = BufWriter::new(file);
        for scalar in witness.iter() {
            let big_int = scalar.into_bigint();
            writeln!(writer, "{}", big_int.to_string()).expect("error writing the witness into a file");
        }
    }
}

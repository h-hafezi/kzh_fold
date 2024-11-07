use crate::commitment;
use std::borrow::Borrow;
use std::marker::PhantomData;
use std::ops::Add;

use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::CurveConfig;
use ark_ff::{BigInteger64, Field, PrimeField};
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::groups::CurveVar;
use ark_r1cs_std::{R1CSVar, ToBitsGadget};
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use ark_std::UniformRand;
use ark_std::{end_timer, start_timer};
use rand::thread_rng;

use crate::accumulation::accumulator::{AccInstance, AccSRS};
use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use crate::accumulation_circuit::prover::get_random_prover;
use crate::accumulation_circuit::prover::AccumulatorVerifierCircuitProver;
use crate::accumulation_circuit::randomness_different_formats;
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::gadgets::non_native::util::{cast_field, non_native_to_fpvar};
use crate::gadgets::r1cs::{OvaInstance, RelaxedOvaInstance};
use crate::hash::poseidon::PoseidonHashVar;
use crate::nova::cycle_fold::coprocessor::{synthesize, SecondaryCircuit as SecondaryCircuit};
use crate::nova::cycle_fold::coprocessor_constraints::{OvaInstanceVar, RelaxedOvaInstanceVar};
use crate::transcript::transcript_var::TranscriptVar;

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
    pub ova_auxiliary_input_C: OvaInstance<G2, C2>,
    /// auxiliary input which helps to have T'' = (1-beta) * T + beta * T' without scalar multiplication
    pub ova_auxiliary_input_T: OvaInstance<G2, C2>,
    /// auxiliary input which helps to have E_{temp} = (1-beta) * E + beta * E' without scalar multiplication
    pub ova_auxiliary_input_E_1: OvaInstance<G2, C2>,
    /// auxiliary input which helps to have E'' = E_{temp} + beta * (1-beta) * Q without scalar multiplication
    pub ova_auxiliary_input_E_2: OvaInstance<G2, C2>,

    /// accumulation proof for accumulators
    pub cross_term_error_commitment_Q: Projective<G1>,

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub ova_cross_term_error_commitment_C: Projective<G2>,
    pub ova_cross_term_error_commitment_T: Projective<G2>,
    pub ova_cross_term_error_commitment_E_1: Projective<G2>,
    pub ova_cross_term_error_commitment_E_2: Projective<G2>,

    /// the instance to be folded
    pub current_accumulator_instance: AccInstance<E>,
    /// the running accumulator
    pub running_accumulator_instance: AccInstance<E>,
    /// the result accumulator
    pub final_accumulator_instance: AccInstance<E>,

    /// running cycle fold instance
    pub ova_running_instance: RelaxedOvaInstance<G2, C2>,

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
    pub ova_auxiliary_input_C: OvaInstanceVar<G2, C2>,
    /// auxiliary input which helps to have T'' = (1-beta) * T + beta * T' without scalar multiplication
    pub ova_auxiliary_input_T: OvaInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E_{temp} = (1-beta) * E + beta * E' without scalar multiplication
    pub ova_auxiliary_input_E_1: OvaInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E'' = E_{temp} + beta * (1-beta) * Q without scalar multiplication
    pub ova_auxiliary_input_E_2: OvaInstanceVar<G2, C2>,

    /// the randomness used for taking linear combination and its non-native counterpart
    pub beta_var: FpVar<G1::ScalarField>,
    pub beta_var_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// accumulation proof
    pub cross_term_error_commitment_Q: NonNativeAffineVar<G1>,

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub ova_cross_term_error_commitment_C: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub ova_cross_term_error_commitment_T: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub ova_cross_term_error_commitment_E_1: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub ova_cross_term_error_commitment_E_2: ProjectiveVar<G2, FpVar<G2::BaseField>>,

    pub current_accumulator_instance_var: AccumulatorInstanceVar<G1>,
    pub running_accumulator_instance_var: AccumulatorInstanceVar<G1>,
    pub final_accumulator_instance_var: AccumulatorInstanceVar<G1>,

    pub ova_running_instance: RelaxedOvaInstanceVar<G2, C2>,

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
        let ova_auxiliary_input_C = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_C"),
            || Ok(circuit.map(|e| e.ova_auxiliary_input_C.clone()).unwrap()),
            mode,
        ).unwrap();

        let ova_auxiliary_input_T = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_T"),
            || Ok(circuit.map(|e| e.ova_auxiliary_input_T.clone()).unwrap()),
            mode,
        ).unwrap();

        let ova_auxiliary_input_E_1 = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_1"),
            || Ok(circuit.map(|e| e.ova_auxiliary_input_E_1.clone()).unwrap()),
            mode,
        ).unwrap();

        let ova_auxiliary_input_E_2 = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary_input_E_2"),
            || Ok(circuit.map(|e| e.ova_auxiliary_input_E_2.clone()).unwrap()),
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
        let ova_running_instance = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.ova_running_instance.clone()),
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
            || circuit.map(|e| cast_field::<G1::ScalarField, G1::BaseField>(e.beta.clone())),
            mode,
        ).unwrap();

        // folding proofs for cycle fold and accumulator
        let cross_term_error_commitment_Q = NonNativeAffineVar::new_variable(
            ns!(cs, "Q"),
            || circuit.map(|e| e.cross_term_error_commitment_Q),
            mode,
        ).unwrap();

        let ova_cross_term_error_commitment_C = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.ova_cross_term_error_commitment_C.clone()),
            mode,
        ).unwrap();

        let ova_cross_term_error_commitment_T = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.ova_cross_term_error_commitment_T.clone()),
            mode,
        ).unwrap();

        let ova_cross_term_error_commitment_E_1 = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.ova_cross_term_error_commitment_E_1.clone()),
            mode,
        ).unwrap();

        let ova_cross_term_error_commitment_E_2 = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || circuit.map(|e| e.ova_cross_term_error_commitment_E_2.clone()),
            mode,
        ).unwrap();

        Ok(AccumulatorVerifierVar {
            ova_auxiliary_input_C,
            ova_auxiliary_input_T,
            ova_auxiliary_input_E_1,
            ova_auxiliary_input_E_2,
            beta_var,
            beta_var_non_native,
            cross_term_error_commitment_Q,
            ova_cross_term_error_commitment_C,
            ova_cross_term_error_commitment_T,
            ova_cross_term_error_commitment_E_1,
            ova_cross_term_error_commitment_E_2,
            current_accumulator_instance_var,
            running_accumulator_instance_var,
            final_accumulator_instance_var,
            ova_running_instance,
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
    pub fn accumulate(&self, transcript_var: &mut TranscriptVar<G1::ScalarField>) -> (RelaxedOvaInstanceVar<G2, C2>, &AccumulatorInstanceVar<G1>)
    where
        <G2 as CurveConfig>::BaseField: Absorb,
    {
        // checking beta and non_native beta are consistent
        let beta_bits = self.beta_var_non_native.to_bits_le().unwrap();
        self.beta_var.enforce_equal(&Boolean::le_bits_to_fp_var(beta_bits.as_slice()).unwrap()).unwrap();

        // compute hash and make sure it's consistent with input beta
        transcript_var.append_scalars(b"instance 1", self.current_accumulator_instance_var.to_sponge_field_elements().unwrap().as_slice());
        transcript_var.append_scalars(b"instance 2", self.running_accumulator_instance_var.to_sponge_field_elements().unwrap().as_slice());
        transcript_var.append_scalars(b"Q", self.cross_term_error_commitment_Q.to_sponge_field_elements().unwrap().as_slice());
        transcript_var.challenge_scalar(b"challenge scalar").enforce_equal(&self.beta_var).unwrap();

        // Non-native scalar multiplication: linear combination of C
        let (flag,
            r,
            g1,
            g2,
            C_var
        ) = self.ova_auxiliary_input_C.parse_secondary_io::<G1>().unwrap();
        // g1 == acc.C
        self.running_accumulator_instance_var.C_var.enforce_equal(&g1).unwrap();
        // g2 == instance.C
        self.current_accumulator_instance_var.C_var.enforce_equal(&g2).unwrap();
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).unwrap();
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).unwrap();
        // check out the result C_var is consistent with result_acc
        C_var.enforce_equal(&self.final_accumulator_instance_var.C_var).unwrap();


        // Non-native scalar multiplication: linear combination of T
        let (flag,
            r,
            g1,
            g2,
            T_var
        ) = self.ova_auxiliary_input_T.parse_secondary_io::<G1>().unwrap();
        // g1 == acc.T
        self.running_accumulator_instance_var.T_var.enforce_equal(&g1).unwrap();
        // g2 == instance.C
        self.current_accumulator_instance_var.T_var.enforce_equal(&g2).unwrap();
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).unwrap();
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).unwrap();
        // check out the result T_var is consistent with result_acc
        T_var.enforce_equal(&self.final_accumulator_instance_var.T_var).unwrap();


        // Non-native scalar multiplication: linear combination E_temp = (instance.E * (1-beta) + acc.E * beta)
        let (flag,
            r,
            g1,
            g2,
            E_temp
        ) = self.ova_auxiliary_input_E_1.parse_secondary_io::<G1>().unwrap();
        // g1 == acc.E
        self.running_accumulator_instance_var.E_var.enforce_equal(&g1).unwrap();
        // g2 == instance.E
        self.current_accumulator_instance_var.E_var.enforce_equal(&g2).unwrap();
        // enforce flag to be false
        flag.enforce_equal(&NonNativeFieldVar::zero()).unwrap();
        // check r to be equal to beta
        r.enforce_equal(&self.beta_var_non_native).unwrap();


        // Non-native scalar multiplication: linear combination E'' = E_{temp} + (1-beta) * beta * Q
        let (flag,
            _r,
            g1,
            g2,
            E_var
        ) = self.ova_auxiliary_input_E_2.parse_secondary_io::<G1>().unwrap();
        // g1 == Q
        g1.enforce_equal(&self.cross_term_error_commitment_Q).unwrap();
        // g2 == E_temp
        g2.enforce_equal(&E_temp).unwrap();
        // enforce flag to be true
        flag.enforce_equal(&NonNativeFieldVar::one()).unwrap();
        // check r to be equal to beta
        let _beta_times_beta_minus_one = self.beta_var_non_native.clone() - self.beta_var_non_native.square().unwrap();
        // check out the result E_var is consistent with result_acc
        E_var.enforce_equal(&self.final_accumulator_instance_var.E_var).unwrap();


        let beta_minus_one = FpVar::<G1::ScalarField>::one() - &self.beta_var;

        // Native field operation: linear combination of x
        for i in 0..self.running_accumulator_instance_var.x_var.len() {
            let x_var = &self.beta_var * &self.running_accumulator_instance_var.x_var[i] +
                &beta_minus_one * &self.current_accumulator_instance_var.x_var[i];
            // check out the result b_var is consistent with result_acc
            x_var.enforce_equal(&self.final_accumulator_instance_var.x_var[i]).unwrap();
        }

        // Native field operation: linear combination of x
        for i in 0..self.running_accumulator_instance_var.y_var.len() {
            let y_var = &self.beta_var * &self.running_accumulator_instance_var.y_var[i] +
                &beta_minus_one * &self.current_accumulator_instance_var.y_var[i];
            // check out the result b_var is consistent with result_acc
            y_var.enforce_equal(&self.final_accumulator_instance_var.y_var[i]).unwrap();
        }

        // check out the result z_c_var is consistent with result_acc
        self.final_accumulator_instance_var.z_var.enforce_equal(
            &(&self.beta_var * &self.running_accumulator_instance_var.z_var +
                &beta_minus_one * &self.current_accumulator_instance_var.z_var)
        ).unwrap();

        transcript_var.append_scalars_non_native(b"label", self.ova_auxiliary_input_C.X.as_slice());
        transcript_var.append_scalars_non_native(b"label", self.ova_auxiliary_input_T.X.as_slice());
        transcript_var.append_scalars_non_native(b"label", self.ova_auxiliary_input_E_1.X.as_slice());
        transcript_var.append_scalars_non_native(b"label", self.ova_auxiliary_input_E_2.X.as_slice());
        transcript_var.append_scalars(
            b"label",
            &[
                self.ova_auxiliary_input_C.commitment.x.clone(),
                self.ova_auxiliary_input_C.commitment.x.clone(),
                self.ova_auxiliary_input_C.commitment.x.clone(),
                self.ova_cross_term_error_commitment_C.x.clone(),
                self.ova_cross_term_error_commitment_C.y.clone(),
                self.ova_cross_term_error_commitment_C.z.clone(),
            ],
        );
        transcript_var.append_scalars(
            b"label",
            &[
                self.ova_auxiliary_input_T.commitment.x.clone(),
                self.ova_auxiliary_input_T.commitment.x.clone(),
                self.ova_auxiliary_input_T.commitment.x.clone(),
                self.ova_cross_term_error_commitment_T.x.clone(),
                self.ova_cross_term_error_commitment_T.y.clone(),
                self.ova_cross_term_error_commitment_T.z.clone(),
            ],
        );
        transcript_var.append_scalars(
            b"label",
            &[
                self.ova_auxiliary_input_E_1.commitment.x.clone(),
                self.ova_auxiliary_input_E_1.commitment.x.clone(),
                self.ova_auxiliary_input_E_1.commitment.x.clone(),
                self.ova_cross_term_error_commitment_E_1.x.clone(),
                self.ova_cross_term_error_commitment_E_1.y.clone(),
                self.ova_cross_term_error_commitment_E_1.z.clone(),
            ],
        );
        transcript_var.append_scalars(
            b"label",
            &[
                self.ova_auxiliary_input_E_2.commitment.x.clone(),
                self.ova_auxiliary_input_E_2.commitment.x.clone(),
                self.ova_auxiliary_input_E_2.commitment.x.clone(),
                self.ova_cross_term_error_commitment_E_2.x.clone(),
                self.ova_cross_term_error_commitment_E_2.y.clone(),
                self.ova_cross_term_error_commitment_E_2.z.clone(),
            ],
        );

        let beta_2_non_native = &self.beta_var_non_native * &self.beta_var_non_native;
        let beta_3_non_native = &self.beta_var_non_native * &beta_2_non_native;
        let beta_4_non_native = &self.beta_var_non_native * &beta_3_non_native;

        let final_instance = self.ova_running_instance.fold(
            &[
                (
                    (&self.ova_auxiliary_input_C, None),
                    &self.ova_cross_term_error_commitment_C,
                    &self.beta_var_non_native,
                    &beta_bits
                ),
                (
                    (&self.ova_auxiliary_input_T, None),
                    &self.ova_cross_term_error_commitment_T,
                    &beta_2_non_native,
                    &beta_2_non_native.to_bits_le().unwrap(),
                ),
                (
                    (&self.ova_auxiliary_input_E_1, None),
                    &self.ova_cross_term_error_commitment_E_1,
                    &beta_3_non_native,
                    &beta_3_non_native.to_bits_le().unwrap(),
                ),
                (
                    (&self.ova_auxiliary_input_E_2, None),
                    &self.ova_cross_term_error_commitment_E_2,
                    &beta_4_non_native,
                    &beta_4_non_native.to_bits_le().unwrap(),
                ),
            ]
        ).unwrap();


        // return result of accumulation and final cycle fold instance
        (final_instance, &self.final_accumulator_instance_var)
    }
}

impl<G1, G2, C2> AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    ProjectiveVar<G2, FpVar<<G2 as CurveConfig>::BaseField>>: AllocVar<<C2 as CommitmentScheme<Projective<G2>>>::Commitment, <G2 as CurveConfig>::BaseField>,
{
    pub fn new<E: Pairing>(cs: ConstraintSystemRef<G1::ScalarField>, prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, E::ScalarField>) -> AccumulatorVerifierVar<G1, G2, C2>
    where
        E: Pairing<G1Affine=Affine<G1>, ScalarField=<G1 as CurveConfig>::ScalarField, BaseField=<G1 as CurveConfig>::BaseField>,
        <G2 as CurveConfig>::BaseField: Absorb,
        <G2 as CurveConfig>::ScalarField: Absorb,
    {
        // the randomness in different formats
        let beta_scalar = prover.beta.clone();
        let (_, beta_var, beta_var_non_native) = randomness_different_formats::<E>(cs.clone(), beta_scalar);

        // initialise accumulator variables
        let current_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "current accumulator instance var"),
            || Ok(prover.get_current_acc_instance().clone()),
            AllocationMode::Input,
        ).unwrap();

        let running_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "running accumulator instance var"),
            || Ok(prover.get_running_acc_instance().clone()),
            AllocationMode::Input,
        ).unwrap();

        let final_accumulator_instance_var = AccumulatorInstanceVar::new_variable(
            ns!(cs, "final accumulator instance var"),
            || Ok(prover.compute_result_accumulator_instance()),
            AllocationMode::Input,
        ).unwrap();

        // initialise auxiliary input variables
        let auxiliary_input_C_var = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary input C var"),
            || Ok(prover.compute_auxiliary_input_C().0),
            AllocationMode::Input,
        ).unwrap();

        let auxiliary_input_T_var = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary input T var"),
            || Ok(prover.compute_auxiliary_input_T().0),
            AllocationMode::Input,
        ).unwrap();

        let auxiliary_input_E_1_var = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary input E_1 var"),
            || Ok(prover.compute_auxiliary_input_E_1().0),
            AllocationMode::Input,
        ).unwrap();

        let auxiliary_input_E_2_var = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary input E_2 var"),
            || Ok(prover.compute_auxiliary_input_E_2().0),
            AllocationMode::Input,
        ).unwrap();


        // initialise Q variables
        let Q_var = NonNativeAffineVar::new_variable(
            ns!(cs, "Q var"),
            || Ok(prover.compute_proof_Q()),
            AllocationMode::Input,
        ).unwrap();

        let cycle_fold_proof = prover.compute_cycle_fold_proofs_and_final_instance();

        let com_C_var = ProjectiveVar::new_variable(
            ns!(cs, "com_C_var"),
            || Ok(cycle_fold_proof.0),
            AllocationMode::Input,
        ).unwrap();

        let com_T_var = ProjectiveVar::new_variable(
            ns!(cs, "com_T_var"),
            || Ok(cycle_fold_proof.1),
            AllocationMode::Input,
        ).unwrap();

        let com_E_1_var = ProjectiveVar::new_variable(
            ns!(cs, "com_E_1_var"),
            || Ok(cycle_fold_proof.2),
            AllocationMode::Input,
        ).unwrap();

        let com_E_2_var = ProjectiveVar::new_variable(
            ns!(cs, "com_E_2_var"),
            || Ok(cycle_fold_proof.3),
            AllocationMode::Input,
        ).unwrap();


        // initialise cycle fold running instance var
        let running_cycle_fold_instance_var = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "running cycle fold instance var"),
            || Ok(prover.cycle_fold_running_instance),
            AllocationMode::Input,
        ).unwrap();


        let verifier = AccumulatorVerifierVar {
            ova_auxiliary_input_C: auxiliary_input_C_var,
            ova_auxiliary_input_T: auxiliary_input_T_var,
            ova_auxiliary_input_E_1: auxiliary_input_E_1_var,
            ova_auxiliary_input_E_2: auxiliary_input_E_2_var,
            beta_var,
            beta_var_non_native,
            cross_term_error_commitment_Q: Q_var,
            ova_cross_term_error_commitment_C: com_C_var,
            ova_cross_term_error_commitment_T: com_T_var,
            ova_cross_term_error_commitment_E_1: com_E_1_var,
            ova_cross_term_error_commitment_E_2: com_E_2_var,
            current_accumulator_instance_var,
            running_accumulator_instance_var,
            final_accumulator_instance_var,
            ova_running_instance: running_cycle_fold_instance_var,
            n: prover.n,
            m: prover.m,
        };

        verifier
    }
}


#[cfg(test)]
pub mod tests {
    use std::fmt::Debug;

    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef, SynthesisMode};
    use rand::thread_rng;

    use super::*;
    use crate::constant_for_curves::{ScalarField, C2, E, G1, G2};
    use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
    use crate::pcs::multilinear_pcs::PolynomialCommitmentSRS;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

    // Test helper
    pub fn get_random_acc_verifier_cs() -> ConstraintSystemRef<ScalarField> {
        // a constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // initialise the accumulate verifier circuit
        let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, ScalarField> = get_random_prover();
        let verifier = AccumulatorVerifierVar::<G1, G2, C2>::new::<E>(
            cs.clone(),
            prover.clone(),
        );

        println!("number of constraint for initialisation: {}", cs.num_constraints());

        let mut transcript_var = TranscriptVar::from_transcript(
            cs.clone(),
            prover.initial_transcript.clone(),
        );

        // run the accumulation
        let _ = verifier.accumulate(&mut transcript_var);

        println!("number of constraint after accumulation: {}", cs.num_constraints());

        // assert the constraint system is satisfied
        assert!(cs.is_satisfied().unwrap());

        // these are required to called CRR1CSShape::convert
        cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
        cs.finalize();

        cs
    }

    #[test]
    fn kzh_acc_verifier_circuit_end_to_end_test() {
        let SRS: PolynomialCommitmentSRS<E> = MultilinearPolynomial::setup(18, &mut thread_rng()).unwrap();

        let cs = get_random_acc_verifier_cs();

        println!("number of constraint random cs: {}", cs.num_constraints());
        println!("number of instance random cs: {}", cs.num_instance_variables());
        println!("number of witness random cs: {}", cs.num_witness_variables());

        println!("number of constraint before shape convert: {}", cs.num_constraints());

        println!("number of constraint random cs: {}", cs.num_constraints());
        println!("number of instance random cs: {}", cs.num_instance_variables());
        println!("number of witness random cs: {}", cs.num_witness_variables());

        println!("number of constraint before shape convert: {}", cs.num_constraints());
        /*
        // convert to the corresponding Spartan types
        let shape = CRR1CSShape::<ScalarField>::convert::<G1>(cs.clone());
        let key: CRR1CSKey<E, MultilinearPolynomial<ScalarField>> = CRR1CSKey::new(&SRS, shape.get_num_cons(), shape.get_num_vars());
        // Commitment to w(x) happens here
        let instance: CRR1CSInstance<E, MultilinearPolynomial<ScalarField>> = CRR1CSInstance::convert(cs.clone(), &key.keys.ck);

        let witness = CRR1CSWitness::<ScalarField>::convert(cs.clone());

        // check that the Spartan instance-witness pair is still satisfying
        assert!(is_sat(&shape, &instance, &witness, &key).unwrap());

        //////////// Prover: ////////////////

        let mut prover_transcript = Transcript::new(b"example");

        let (proof, rx, ry) = CRR1CSProof::prove(
            &shape,
            &instance,
            witness,
            &key,
            &mut prover_transcript,
        );

        // evaluate matrices A B C
        let inst_evals = shape.inst.inst.evaluate(&rx, &ry);

        //////////// Verifier ////////////////

        let mut verifier_transcript = Transcript::new(b"example");
        let (num_cons, num_vars, _num_inputs) = (
            shape.get_num_cons(),
            shape.get_num_vars(),
            shape.get_num_inputs(),
        );
        assert!(proof
            .verify(
                num_vars,
                num_cons,
                &instance,
                &inst_evals,
                &mut verifier_transcript,
            )
            .is_ok());
         */
    }
}
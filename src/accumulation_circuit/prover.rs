use crate::accumulation::accumulator::{AccInstance, AccSRS, Accumulator};
use crate::accumulation_circuit::affine_to_projective;
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::util::cast_field;
use crate::gadgets::r1cs::ova::commit_T;
use crate::gadgets::r1cs::{OvaInstance, OvaWitness, R1CSShape, RelaxedOvaInstance, RelaxedOvaWitness};
use crate::hash::pederson::PedersenCommitment;
use crate::nova::cycle_fold::coprocessor::{setup_shape, synthesize, SecondaryCircuit};
use crate::kzh::kzh2::{KZH2, KZH2SRS};
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::{CurveConfig, CurveGroup};
use ark_ff::Field;
use ark_ff::PrimeField;
use rand::thread_rng;
use crate::kzh::KZH;
use crate::math::Math;

#[derive(Clone)]
pub struct AccumulatorVerifierCircuitProver<G1, G2, C2, E, F>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    /// the randomness used for taking linear combination, between two accumulations
    /// it should be input from Accumulator::compute_fiat_shamir_challenge()
    pub beta: F,

    pub initial_transcript: Transcript<F>,
    pub final_transcript: Transcript<F>,

    /// srs for the accumulation
    pub srs: AccSRS<E>,

    /// the instance to be folded
    pub current_accumulator: Accumulator<E>,
    /// the running accumulator
    pub running_accumulator: Accumulator<E>,

    /// running cycle fold instance
    pub shape: R1CSShape<G2>,
    pub commitment_pp: <C2 as CommitmentScheme<Projective<G2>>>::PP,
    pub cycle_fold_running_instance: RelaxedOvaInstance<G2, C2>,
    pub cycle_fold_running_witness: RelaxedOvaWitness<G2>,

    // these are constant values
    pub n: u32,
    pub m: u32,
}

impl<G1, G2, C2, E, F> AccumulatorVerifierCircuitProver<G1, G2, C2, E, F>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    #[inline(always)]
    pub fn get_current_acc_instance(&self) -> &AccInstance<E> {
        &self.current_accumulator.instance
    }

    #[inline(always)]
    pub fn get_running_acc_instance(&self) -> &AccInstance<E> {
        &self.running_accumulator.instance
    }
}

impl<G1, G2, C2, E, F> AccumulatorVerifierCircuitProver<G1, G2, C2, E, F>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    pub fn is_satisfied(&self)
    where
        <G2 as CurveConfig>::ScalarField: Absorb,
    {
        assert!(Accumulator::decide(&self.srs, &self.running_accumulator));
        assert!(Accumulator::decide(&self.srs, &self.current_accumulator));
        self.shape.is_relaxed_ova_satisfied(
            &self.cycle_fold_running_instance,
            &self.cycle_fold_running_witness,
            &self.commitment_pp,
        ).expect("panic!");
    }

    pub fn compute_auxiliary_input_C(&self) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let g1 = affine_to_projective(self.running_accumulator.instance.C.clone());
        let g2 = affine_to_projective(self.current_accumulator.instance.C.clone());

        // C'' = beta * acc_running.instance.C + (1 - beta) * acc_instance.instance.C
        let g_out = (g1 * self.beta) + (g2 * (G1::ScalarField::ONE - self.beta));

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1,
            g2,
            g_out,
            r: cast_field::<G1::ScalarField, G1::BaseField>(self.beta),
            flag: false,
        }, &self.commitment_pp[0..self.shape.num_vars].to_vec()).unwrap()
    }

    pub fn compute_auxiliary_input_T(&self) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let g1 = affine_to_projective(self.running_accumulator.instance.T.clone());
        let g2 = affine_to_projective(self.current_accumulator.instance.T.clone());

        // T'' = beta * acc_running.instance.T + (1 - beta) * acc_instance.instance.T
        let g_out = (g1 * self.beta) + (g2 * (F::ONE - self.beta));

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1,
            g2,
            g_out,
            r: cast_field::<G1::ScalarField, G1::BaseField>(self.beta),
            flag: false,
        }, &self.commitment_pp[0..self.shape.num_vars].to_vec()).unwrap()
    }

    pub fn compute_auxiliary_input_E_1(&self) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let g1 = affine_to_projective(self.running_accumulator.instance.E.clone());
        let g2 = affine_to_projective(self.current_accumulator.instance.E.clone());

        // E_temp = beta * acc_running.instance.E + (1 - beta) * acc_instance.instance.E
        let g_out = (g1 * self.beta) + (g2 * (F::ONE - self.beta));

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1,
            g2,
            g_out,
            r: cast_field::<G1::ScalarField, G1::BaseField>(self.beta),
            flag: false,
        }, &self.commitment_pp[0..self.shape.num_vars].to_vec()).unwrap()
    }

    pub fn compute_auxiliary_input_E_2(&self) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let e1 = affine_to_projective(self.running_accumulator.instance.E.clone());
        let e2 = affine_to_projective(self.current_accumulator.instance.E.clone());

        // E_temp = beta * e1 + (1 - beta) * e2
        let E_temp = (e1 * self.beta) + (e2 * (F::ONE - self.beta));
        let Q = self.compute_proof_Q();
        let g_out = E_temp + Q * (self.beta * (F::ONE - self.beta));

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1: Q,
            g2: E_temp,
            g_out,
            r: cast_field::<G1::ScalarField, G1::BaseField>(self.beta * (F::ONE - self.beta)),
            flag: true,
        }, &self.commitment_pp[0..self.shape.num_vars].to_vec()).unwrap()
    }

    pub fn compute_proof_Q(&self) -> Projective<G1>
    {
        // since acc_instance takes (1- beta) then it should be first in the function argument
        affine_to_projective(Accumulator::helper_function_Q(&self.srs,
                                                            &self.current_accumulator,
                                                            &self.running_accumulator)
        )
    }

    pub fn compute_result_accumulator_instance(&self) -> AccInstance<E>
    {
        let mut transcript = self.initial_transcript.clone();
        Accumulator::prove(&self.srs,
                           &self.current_accumulator,
                           &self.running_accumulator,
                           &mut transcript,
        ).0
    }

    pub fn compute_cycle_fold_proofs_and_final_instance(&self) -> (
        C2::Commitment,
        C2::Commitment,
        C2::Commitment,
        C2::Commitment,
        RelaxedOvaInstance<G2, C2>
    )
    where
        <G2 as CurveConfig>::ScalarField: Absorb,
    {
        let compute_commit_and_fold =
            |running_witness: &RelaxedOvaWitness<G2>,
             running_instance: &RelaxedOvaInstance<G2, C2>,
             witness: &OvaWitness<G2>,
             instance: &OvaInstance<G2, C2>,
             beta: &G2::ScalarField,
            | -> (C2::Commitment, RelaxedOvaWitness<G2>, RelaxedOvaInstance<G2, C2>) {
                let (T, com_T) = commit_T(
                    &self.shape,
                    &self.commitment_pp[self.shape.num_vars..].to_vec(),
                    running_instance,
                    running_witness,
                    instance,
                    witness,
                ).unwrap();

                // Fold the running instance and witness with the first proof
                let new_running_instance = running_instance.fold(instance, &com_T, beta).unwrap();
                let new_running_witness = running_witness.fold(witness, &T, beta).unwrap();

                (com_T, new_running_witness, new_running_instance)
            };

        let beta_non_native = cast_field::<G1::ScalarField, G1::BaseField>(self.beta);

        // first fold auxiliary_input_C with the running instance
        let (instance_C, witness_C) = self.compute_auxiliary_input_C();
        let (com_C, new_running_witness, new_running_instance) = compute_commit_and_fold(
            &self.cycle_fold_running_witness,
            &self.cycle_fold_running_instance,
            &witness_C,
            &instance_C,
            &beta_non_native,
        );

        self.shape.is_ova_satisfied(&instance_C, &witness_C, &self.commitment_pp).unwrap();
        self.shape.is_relaxed_ova_satisfied(&new_running_instance, &new_running_witness, &self.commitment_pp).unwrap();


        // first fold auxiliary_input_T with the running instance
        let beta_2: <G2 as CurveConfig>::ScalarField = beta_non_native * beta_non_native;
        let (instance_T, witness_T) = self.compute_auxiliary_input_T();
        let (com_T, new_running_witness, new_running_instance) = compute_commit_and_fold(
            &new_running_witness,
            &new_running_instance,
            &witness_T,
            &instance_T,
            &beta_2,
        );

        self.shape.is_ova_satisfied(&instance_T, &witness_T, &self.commitment_pp).unwrap();
        self.shape.is_relaxed_ova_satisfied(&new_running_instance, &new_running_witness, &self.commitment_pp).unwrap();

        // first fold auxiliary_input_E_1 with the running instance
        let (instance_E_1, witness_E_1) = self.compute_auxiliary_input_E_1();
        let beta_3: <G2 as CurveConfig>::ScalarField = beta_2 * beta_non_native;
        let (com_E_1, new_running_witness, new_running_instance) = compute_commit_and_fold(
            &new_running_witness,
            &new_running_instance,
            &witness_E_1,
            &instance_E_1,
            &beta_3,
        );

        self.shape.is_ova_satisfied(&instance_E_1, &witness_E_1, &self.commitment_pp).unwrap();
        self.shape.is_relaxed_ova_satisfied(&new_running_instance, &new_running_witness, &self.commitment_pp).unwrap();

        // first fold auxiliary_input_E_1 with the running instance
        let (instance_E_2, witness_E_2) = self.compute_auxiliary_input_E_2();
        let beta_4: <G2 as CurveConfig>::ScalarField = beta_3 * beta_non_native;
        let (com_E_2, new_running_witness, new_running_instance) = compute_commit_and_fold(
            &new_running_witness,
            &new_running_instance,
            &witness_E_2,
            &instance_E_2,
            &beta_4,
        );

        self.shape.is_ova_satisfied(&instance_E_2, &witness_E_2, &self.commitment_pp).unwrap();
        self.shape.is_relaxed_ova_satisfied(&new_running_instance, &new_running_witness, &self.commitment_pp).unwrap();

        (com_C, com_T, com_E_1, com_E_2, new_running_instance)
    }

    pub fn new(srs: &AccSRS<E>,
               commitment_pp: C2::PP,
               running_accumulator: Accumulator<E>,
               current_accumulator: Accumulator<E>,
               cycle_fold_running_instance: RelaxedOvaInstance<G2, C2>,
               cycle_fold_running_witness: RelaxedOvaWitness<G2>,
               initial_transcript: Transcript<F>,
    ) -> AccumulatorVerifierCircuitProver<G1, G2, C2, E, F> {
        // assert accumulators are satisfied
        debug_assert!(Accumulator::decide(&srs, &running_accumulator));
        debug_assert!(Accumulator::decide(&srs, &current_accumulator));

        // the shape of the R1CS instance
        let shape = setup_shape::<G1, G2>().unwrap();

        // assert the length of the commitment__pp is well formatted
        assert_eq!(commitment_pp.len(), shape.num_vars + shape.num_constraints);

        // assert the cycle fold instance/witness is satisfied
        shape.is_relaxed_ova_satisfied(
            &cycle_fold_running_instance,
            &cycle_fold_running_witness,
            &commitment_pp,
        ).expect("error doing equality assertion");

        // compute Q
        let Q = Accumulator::helper_function_Q(&srs, &current_accumulator, &running_accumulator);

        // clone the transcript so it doesn't change
        let mut transcript = initial_transcript.clone();
        let beta = Accumulator::compute_fiat_shamir_challenge(
            &mut transcript,
            &current_accumulator.instance,
            &running_accumulator.instance,
            Q,
        );

        AccumulatorVerifierCircuitProver {
            beta,
            initial_transcript,
            final_transcript: transcript,
            srs: srs.clone(),
            current_accumulator,
            running_accumulator,
            shape,
            commitment_pp,
            cycle_fold_running_instance,
            cycle_fold_running_witness,
            n: srs.pc_srs.degree_x as u32,
            m: srs.pc_srs.degree_y as u32,
        }
    }

    pub fn get_trivial_cycle_fold_running_instance_witness(shape: &R1CSShape<G2>) -> (RelaxedOvaInstance<G2, C2>, RelaxedOvaWitness<G2>) {
        let cycle_fold_running_instance = RelaxedOvaInstance::new(&shape);
        let cycle_fold_running_witness = RelaxedOvaWitness::zero(&shape);

        (cycle_fold_running_instance, cycle_fold_running_witness)
    }

    pub fn get_commitment_pp(shape: &R1CSShape<G2>) -> C2::PP {
        // public parameters of Pedersen
        let commitment_pp: Vec<Affine<G2>> = PedersenCommitment::<Projective<G2>>::setup(shape.num_vars + shape.num_constraints, b"test", &());

        commitment_pp
    }
}

pub fn get_random_prover<G1, G2, C2, E, F>() -> AccumulatorVerifierCircuitProver<G1, G2, C2, E, F>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    // specifying degrees of polynomials
    let (n, m) = (4, 4);

    // get a random srs
    let srs = {
        let srs_pcs: KZH2SRS<E> = KZH2::setup((n * m as usize).log_2(), &mut thread_rng());
        Accumulator::setup(srs_pcs.clone(), &mut thread_rng())
    };

    // the shape of the R1CS instance
    let shape = setup_shape::<G1, G2>().unwrap();

    // get trivial running instance
    let (cycle_fold_running_instance, cycle_fold_running_witness) = AccumulatorVerifierCircuitProver::<G1, G2, C2, E, F>::get_trivial_cycle_fold_running_instance_witness(&shape);

    // get commitment_pp
    let commitment_pp = AccumulatorVerifierCircuitProver::<G1, G2, C2, E, F>::get_commitment_pp(&shape);

    // get two random accumulators
    let current_accumulator = Accumulator::rand(&srs, &mut thread_rng());
    let running_accumulator = Accumulator::rand(&srs, &mut thread_rng());


    let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, F> = AccumulatorVerifierCircuitProver::new(
        &srs,
        commitment_pp,
        running_accumulator,
        current_accumulator,
        cycle_fold_running_instance,
        cycle_fold_running_witness,
        Transcript::new(b"initial_transcript"),
    );

    // assert it's formated correctly
    prover.is_satisfied();

    // return the prover
    prover
}


#[cfg(test)]
pub mod tests {
    use ark_ec::pairing::Pairing;
    use ark_ec::CurveConfig;
    use ark_ff::Field;

    use crate::accumulation::accumulator::Accumulator;
    use crate::accumulation_circuit::prover::{get_random_prover, AccumulatorVerifierCircuitProver};
    use crate::constant_for_curves::{BaseField, ScalarField, E, G1, G2};
    use crate::gadgets::non_native::util::cast_field;
    use crate::hash::pederson::PedersenCommitment;

    type GrumpkinCurveGroup = ark_grumpkin::Projective;
    type C2 = PedersenCommitment<GrumpkinCurveGroup>;
    type F = ScalarField;
    type Q = BaseField;

    #[test]
    pub fn auxiliary_input_C_correctness() {
        let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, F> = get_random_prover();

        let (r1cs_instance, _) = prover.compute_auxiliary_input_C();
        let secondary_circuit = r1cs_instance.parse_secondary_io().unwrap();

        // get the accumulated result
        let new_acc_instance = Accumulator::prove(
            &prover.srs,
            &prover.current_accumulator,
            &prover.running_accumulator,
            &mut prover.initial_transcript.clone(),
        ).0;

        assert_eq!(secondary_circuit.r, cast_field::<F, BaseField>(prover.beta));
        assert_eq!(secondary_circuit.flag, false);
        assert_eq!(secondary_circuit.g1, prover.running_accumulator.instance.C);
        assert_eq!(secondary_circuit.g2, prover.current_accumulator.instance.C);
        assert_eq!(secondary_circuit.g_out, new_acc_instance.C);
    }

    #[test]
    pub fn auxiliary_input_T_correctness() {
        let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, F> = get_random_prover();

        let (r1cs_instance, _) = prover.compute_auxiliary_input_T();
        let secondary_circuit = r1cs_instance.parse_secondary_io().unwrap();

        // get the accumulated result
        let new_acc_instance = Accumulator::prove(
            &prover.srs,
            &prover.current_accumulator,
            &prover.running_accumulator,
            &mut prover.initial_transcript.clone(),
        ).0;

        assert_eq!(secondary_circuit.r, cast_field::<F, BaseField>(prover.beta));
        assert_eq!(secondary_circuit.flag, false);
        assert_eq!(secondary_circuit.g1, prover.running_accumulator.instance.T);
        assert_eq!(secondary_circuit.g2, prover.current_accumulator.instance.T);
        assert_eq!(secondary_circuit.g_out, new_acc_instance.T);
    }


    #[test]
    pub fn auxiliary_input_E_correctness() {
        let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, F> = get_random_prover();

        let (r1cs_instance, _) = prover.compute_auxiliary_input_E_1();
        let secondary_circuit_E_1 = r1cs_instance.parse_secondary_io().unwrap();

        let (r1cs_instance, _) = prover.compute_auxiliary_input_E_2();
        let secondary_circuit_E_2 = r1cs_instance.parse_secondary_io().unwrap();

        let Q = prover.compute_proof_Q();

        // get the accumulated result
        let new_acc_instance = Accumulator::prove(
            &prover.srs,
            &prover.current_accumulator,
            &prover.running_accumulator,
            &mut prover.initial_transcript.clone(),
        ).0;

        // checking correctness of flags
        assert_eq!(secondary_circuit_E_1.flag, false);
        assert_eq!(secondary_circuit_E_2.flag, true);

        // checking correctness of randomness
        assert_eq!(secondary_circuit_E_1.r, cast_field::<F, Q>(prover.beta));
        assert_eq!(secondary_circuit_E_2.r, cast_field::<F, Q>(prover.beta * (F::ONE - prover.beta)));

        // check E_temp is present in two circuits
        assert_eq!(secondary_circuit_E_1.g_out, secondary_circuit_E_2.g2);

        // check input to the first circuit is correct
        assert_eq!(secondary_circuit_E_1.g1, prover.running_accumulator.instance.E);
        assert_eq!(secondary_circuit_E_1.g2, prover.current_accumulator.instance.E);

        // check input to the first circuit is correct
        assert_eq!(secondary_circuit_E_2.g1, Q);
    }

    #[test]
    pub fn compute_cycle_fold_proofs_correctness() {
        let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, F> = get_random_prover();
        let _ = prover.compute_cycle_fold_proofs_and_final_instance();
    }
}

use crate::accumulation_circuit::affine_to_projective;
use crate::commitment::{Commitment, CommitmentScheme};
use crate::gadgets::non_native::util::cast_field;
use crate::gadgets::r1cs::conversion::{convert_constraint_system_into_instance_witness, get_random_r1cs_instance_witness, get_random_relaxed_r1cs_instance_witness};
use crate::gadgets::r1cs::ova::commit_T as Ova_commit_T;
use crate::gadgets::r1cs::r1cs::commit_T as R1CS_commit_T;
use crate::gadgets::r1cs::{OvaInstance, OvaWitness, R1CSInstance, R1CSShape, R1CSWitness, RelaxedOvaInstance, RelaxedOvaWitness, RelaxedR1CSInstance, RelaxedR1CSWitness};
use crate::nova::cycle_fold::coprocessor::{setup_shape, synthesize, SecondaryCircuit};
use crate::nova::nova::get_affine_coords;
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::{CurveConfig, CurveGroup};
use ark_ff::PrimeField;
use ark_relations::r1cs::ConstraintSystemRef;
use rand::thread_rng;

pub struct NovaProver<F, G1, G2, C1, C2>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    /// shape of the main curve
    pub shape: R1CSShape<G1>,

    /// srs for the kzh_fold
    pub commitment_pp: <C1 as CommitmentScheme<Projective<G1>>>::PP,

    /// the instance to be folded
    pub current_accumulator: (R1CSInstance<G1, C1>, R1CSWitness<G1>),

    /// the running accumulator
    pub running_accumulator: (RelaxedR1CSInstance<G1, C1>, RelaxedR1CSWitness<G1>),

    /// running cycle fold instance
    pub ova_shape: R1CSShape<G2>,
    pub ova_commitment_pp: <C2 as CommitmentScheme<Projective<G2>>>::PP,
    pub ova_running_instance: RelaxedOvaInstance<G2, C2>,
    pub ova_running_witness: RelaxedOvaWitness<G2>,
}

impl<F, G1, G2, C1, C2> NovaProver<F, G1, G2, C1, C2>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>, Commitment=Projective<G1>, SetupAux=()>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>, SetupAux=()>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    pub fn compute_nova_cross_term_error(&self) -> Projective<G1> {
        let (_, com_t) = R1CS_commit_T(&self.shape,
                                       &self.commitment_pp,
                                       &self.running_accumulator.0,
                                       &self.running_accumulator.1,
                                       &self.current_accumulator.0,
                                       &self.current_accumulator.1,
        ).unwrap();

        com_t
    }

    pub fn compute_beta(&self) -> (F, Transcript<F>) {
        // turn the cross term error for nova into affine
        let affine: Affine<G1> = CurveGroup::into_affine(self.compute_nova_cross_term_error());
        let (com_T_x, com_T_y) = get_affine_coords(&affine);

        // make a new transcript and add with the following order: running accumulator instance + current accumulator instance + cross term error
        let mut transcript = Transcript::new(b"new transcript");
        transcript.append_scalars(b"label", self.running_accumulator.0.to_sponge_field_elements().as_slice());
        transcript.append_scalars(b"label", self.current_accumulator.0.to_sponge_field_elements().as_slice());
        transcript.append_scalars_non_native(b"label", &[com_T_x, com_T_y]);

        // derive the challenge
        let beta = transcript.challenge_scalar(b"challenge");

        // this transcript is then used to compute cycle fold challenges
        (beta, transcript)
    }

    pub fn compute_final_accumulator(&self, beta: &F) -> (RelaxedR1CSInstance<G1, C1>, RelaxedR1CSWitness<G1>, Projective<G1>) {
        let (nova_cross_term_error, nova_cross_term_error_commitment) = R1CS_commit_T(
            &self.shape,
            &self.commitment_pp,
            &self.running_accumulator.0,
            &self.running_accumulator.1,
            &self.current_accumulator.0,
            &self.current_accumulator.1,
        ).unwrap();

        // folding two instances
        let folded_instance = self.running_accumulator.0.fold(
            &self.current_accumulator.0,
            &nova_cross_term_error_commitment,
            &beta,
        ).unwrap();

        // folding two witnesses
        let folded_witness = self.running_accumulator.1.fold(
            &self.current_accumulator.1,
            &nova_cross_term_error,
            &beta,
        ).unwrap();

        (folded_instance, folded_witness, nova_cross_term_error_commitment)
    }

    pub fn compute_ova_auxiliary_input_E(&self, beta: &F) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let g1 = self.compute_nova_cross_term_error();
        let g2 = affine_to_projective(self.running_accumulator.0.commitment_E.into());

        let g_out = (g1 * beta) + g2;

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1,
            g2,
            g_out,
            r: cast_field::<G1::ScalarField, G1::BaseField>(*beta),
            flag: true,
        }, &self.ova_commitment_pp[0..self.ova_shape.num_vars].to_vec()).unwrap()
    }

    pub fn compute_ova_auxiliary_input_W(&self, beta: &F) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let g1 = affine_to_projective(self.current_accumulator.0.commitment_W.into());
        let g2 = affine_to_projective(self.running_accumulator.0.commitment_W.into());
        let g_out = (g1 * beta) + g2;

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1,
            g2,
            g_out,
            r: cast_field::<G1::ScalarField, G1::BaseField>(*beta),
            flag: true,
        }, &self.ova_commitment_pp[0..self.ova_shape.num_vars].to_vec()).unwrap()
    }

    pub fn compute_ova_final_instance(&self) -> ((RelaxedOvaInstance<G2, C2>, RelaxedOvaWitness<G2>), (C2::Commitment, C2::Commitment), (F, F)) {
        // get the random challenge beta
        let (beta, mut transcript) = self.compute_beta();

        // get the ova instance/witness for computing the new witness
        let (ova_instance_w, ova_witness_w) = self.compute_ova_auxiliary_input_W(&beta);

        // get the cross term proof
        let (cross_term_error_w, cross_term_error_commitment_w) = Ova_commit_T(
            &self.ova_shape,
            &self.ova_commitment_pp[self.ova_shape.num_vars..].to_vec(),
            &self.ova_running_instance,
            &self.ova_running_witness,
            &ova_instance_w,
            &ova_witness_w,
        ).unwrap();

        // add the new cross term error to the transcript
        let coordinates = get_affine_coords::<G2::BaseField, G2>(&cross_term_error_commitment_w.into_affine());
        transcript.append_scalars(b"add scalars", &[
            coordinates.0,
            coordinates.1,
        ]);

        // derive beta_1
        let beta_1 = transcript.challenge_scalar(b"challenge");

        // currently we use the same beta as randomness, this can later change
        let beta_1_non_native = cast_field::<G1::ScalarField, G1::BaseField>(beta_1);

        // compute the folded instance and folded witness
        let folded_instance = self.ova_running_instance.fold(
            &ova_instance_w,
            &cross_term_error_commitment_w,
            &beta_1_non_native,
        ).expect("folding instance error");

        let folded_witness = self.ova_running_witness.fold(
            &ova_witness_w,
            &cross_term_error_w,
            &beta_1_non_native,
        ).expect("folding witness error");

        // compute ova instance / witness for the error term
        let (ova_instance_e, ova_witness_e) = self.compute_ova_auxiliary_input_E(&beta);

        // compute the cross term error
        let (cross_term_error_e, cross_term_error_commitment_e) = Ova_commit_T(
            &self.ova_shape,
            &self.ova_commitment_pp[self.ova_shape.num_vars..].to_vec(),
            &folded_instance,
            &folded_witness,
            &ova_instance_e,
            &ova_witness_e,
        ).unwrap();

        // add the new cross term error to the transcript
        let coordinates = get_affine_coords::<G2::BaseField, G2>(&cross_term_error_commitment_e.into_affine());
        transcript.append_scalars(b"add scalars", &[
            coordinates.0,
            coordinates.1,
        ]);

        // derive beta_2
        let beta_2 = transcript.challenge_scalar(b"challenge");

        // currently we use the same beta as randomness, this can later change
        let beta_2_non_native = cast_field::<G1::ScalarField, G1::BaseField>(beta_2);

        // compute the next folded instance / witness
        let final_folded_instance = folded_instance.fold(
            &ova_instance_e,
            &cross_term_error_commitment_e,
            &beta_2_non_native,
        ).expect("folding instance error");

        let final_folded_witness = folded_witness.fold(
            &ova_witness_e,
            &cross_term_error_e,
            &beta_2_non_native,
        ).expect("folding witness error");

        (
            (final_folded_instance, final_folded_witness),
            (cross_term_error_commitment_w, cross_term_error_commitment_e),
            (beta_1, beta_2)
        )
    }

    pub fn rand(structure: (usize, usize, usize)) -> NovaProver<F, G1, G2, C1, C2> {
        // the shape of the secondary curve R1CS instance
        let ova_shape = setup_shape::<G1, G2>().unwrap();
        // the main shape
        let (num_constraints, num_io, num_vars) = structure;

        // get commitment_pp
        let ova_commitment_pp: Vec<Affine<G2>> = C2::setup(ova_shape.num_vars + ova_shape.num_constraints, b"test", &());
        let commitment_pp: Vec<Affine<G1>> = C1::setup(num_vars, b"test", &());

        let (shape, instance, witness) = get_random_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &commitment_pp);

        // assert it's satisfied
        shape.is_satisfied(&instance, &witness, &commitment_pp).expect("unsatisfied r1cs");

        // generate a relaxed instance/witness this time
        let (relaxed_shape, relaxed_instance, relaxed_witness) = get_random_relaxed_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &commitment_pp);

        // assert the shape is equal to the previous shape
        assert_eq!(shape, relaxed_shape);

        // make sure the instance is satisfied
        shape.is_relaxed_satisfied(&relaxed_instance, &relaxed_witness, &commitment_pp).expect("unsatisfied r1cs");

        NovaProver {
            shape,
            commitment_pp,
            current_accumulator: (instance, witness),
            running_accumulator: (relaxed_instance, relaxed_witness),
            ova_commitment_pp,
            ova_running_instance: RelaxedOvaInstance::new(&ova_shape),
            ova_running_witness: RelaxedOvaWitness::zero(&ova_shape),
            ova_shape,
        }
    }

    pub fn rand_from_constraint_system(cs: ConstraintSystemRef<F>) -> NovaProver<F, G1, G2, C1, C2> {
        // the shape of the secondary curve R1CS instance
        let ova_shape = setup_shape::<G1, G2>().unwrap();

        // get commitment_pp
        let ova_commitment_pp: Vec<Affine<G2>> = C2::setup(ova_shape.num_vars + ova_shape.num_constraints, b"test", &());
        let commitment_pp: Vec<Affine<G1>> = C1::setup(cs.num_witness_variables(), b"test", &());

        let (shape, instance, witness) = convert_constraint_system_into_instance_witness::<F, C1, G1>(
            cs.clone(),
            &commitment_pp,
        );

        // assert it's satisfied
        shape.is_satisfied(&instance, &witness, &commitment_pp).expect("unsatisfied r1cs");

        // generate a relaxed instance/witness this time
        let (relaxed_instance, relaxed_witness) = shape.random_relaxed_r1cs(
            &commitment_pp,
            &mut thread_rng(),
        );

        // make sure the instance is satisfied
        shape.is_relaxed_satisfied(&relaxed_instance, &relaxed_witness, &commitment_pp).expect("unsatisfied r1cs");

        NovaProver {
            shape,
            commitment_pp,
            current_accumulator: (instance, witness),
            running_accumulator: (relaxed_instance, relaxed_witness),
            ova_commitment_pp,
            ova_running_instance: RelaxedOvaInstance::new(&ova_shape),
            ova_running_witness: RelaxedOvaWitness::zero(&ova_shape),
            ova_shape,
        }
    }
}


#[cfg(test)]
mod test {
    use ark_ec::CurveGroup;
    use crate::constant_for_curves::{ScalarField, C1, C2, E, G1, G2};
    use crate::hash::poseidon::PoseidonHashVar;
    use crate::nova::nova::verifier_circuit::NovaVerifierCircuit;
    use crate::nova::nova::verifier_circuit_var::NovaVerifierCircuitVar;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisMode};
    use ark_serialize::CanonicalSerialize;
    use ark_std::{end_timer, start_timer, UniformRand};
    use rand::thread_rng;
    use crate::gadgets::r1cs::r1cs::commit_T;
    use crate::nova::nova::prover::NovaProver;
    use crate::commitment::{CommitmentScheme};
    use crate::transcript::transcript::Transcript;

    type F = ScalarField;

    #[test]
    fn test_nova_compute_final_accumulator() {
        let prover: NovaProver<F, G1, G2, C1, C2> = NovaProver::rand((10, 3, 17));

        let beta = F::rand(&mut thread_rng());
        let folded_accumulator = prover.compute_final_accumulator(&beta);

        prover.shape.is_relaxed_satisfied(&folded_accumulator.0, &folded_accumulator.1, &prover.commitment_pp).unwrap();

        let circuit_e = prover.compute_ova_auxiliary_input_E(&beta);
        let circuit_w = prover.compute_ova_auxiliary_input_W(&beta);

        // assert both instances are satisfied
        prover.ova_shape.is_ova_satisfied(&circuit_e.0, &circuit_e.1, &prover.ova_commitment_pp).unwrap();
        prover.ova_shape.is_ova_satisfied(&circuit_w.0, &circuit_w.1, &prover.ova_commitment_pp).unwrap();

        // make sure they're consistent with the folded_accumulator
        let secondary_circuit_W = circuit_w.0.parse_secondary_io::<G1>().unwrap();
        assert_eq!(secondary_circuit_W.g_out, folded_accumulator.0.commitment_W);
        let secondary_circuit_E = circuit_e.0.parse_secondary_io::<G1>().unwrap();
        assert_eq!(secondary_circuit_E.g_out, folded_accumulator.0.commitment_E);

        let ((final_ova_instance, final_ova_witness), (_, _), (_, _)) = prover.compute_ova_final_instance();
        prover.ova_shape.is_relaxed_ova_satisfied(&final_ova_instance, &final_ova_witness, &prover.ova_commitment_pp).unwrap();

        println!("ova instance len: {} bytes", final_ova_instance.compressed_size());
        println!("ova witness len: {} bytes", final_ova_witness.compressed_size());
    }

    #[test]
    fn bench_nova() {
        const POSEIDON_NUM: usize = 0;

        // generate some random prover which for a small R1CS which is supposed to generate some sample NovaVerifierCircuit
        // this prover is sort of faking the previous step of IVC, but since we don't have it I just use a prover for a small R1CS
        let prover: NovaProver<F, G1, G2, C1, C2> = NovaProver::rand((10, 3, 17));

        // generate a constraint system
        let cs = ConstraintSystem::<F>::new_ref();

        // generate the Nova verifier circuit
        let nova_verifier_var: NovaVerifierCircuitVar<F, G1, G2, C1, C2> = {
            // the non-circuit version
            let nova_verifier: NovaVerifierCircuit<F, G1, G2, C1, C2> = NovaVerifierCircuit::initialise(prover);
            nova_verifier.verify();

            // the circuit version and output it
            let nova_verifier_var: NovaVerifierCircuitVar<F, G1, G2, C1, C2> = NovaVerifierCircuitVar::new_variable(
                cs.clone(),
                || Ok(nova_verifier.clone()),
                AllocationMode::Input,
            ).unwrap();

            nova_verifier_var
        };

        // generate constraints
        nova_verifier_var.generate_constraints(cs.clone()).expect("error");

        println!("Num of constraints for Nova verifier: {}, cs_satisfied: {}", cs.num_constraints(), cs.is_satisfied().unwrap());

        // pad it with some random poseidon hash
        let mut hash = PoseidonHashVar::new(cs.clone());
        for _ in 0..POSEIDON_NUM {
            // get a random element
            let r = FpVar::new_variable(cs.clone(), || Ok(F::rand(&mut thread_rng())), AllocationMode::Witness).unwrap();
            // update sponge with this random element
            hash.update_sponge(vec![r]);
            // output the hash
            let _ = hash.output();
        }

        println!("Number of constraints of the augmented circuit: {}", cs.num_constraints());

        // the cs has to be finalised before the matrices A, B and C are generated
        cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
        cs.finalize();

        // now generate a prover for the augmented circuit
        let prover: NovaProver<F, G1, G2, C1, C2> = NovaProver::rand_from_constraint_system(cs.clone());

        // now prover include the running instance + the final instance, so we can benchmark each of these

        // ********************************** Benchmarking The Prover **********************************
        // prover cost: commit to current witness + compute and commit to cross term error
        let prover_timer = start_timer!(|| "prover");

        let commit_w_timer = start_timer!(|| "commit witness");
        let _ = C1::commit(&prover.commitment_pp, prover.current_accumulator.1.W.as_slice());
        end_timer!(commit_w_timer);

        let commit_t_timer = start_timer!(|| "commit crossterm");
        let (_, com_t) = commit_T(
            &prover.shape,
            &prover.commitment_pp,
            &prover.running_accumulator.0,
            &prover.running_accumulator.1,
            &prover.current_accumulator.0,
            &prover.current_accumulator.1
        ).unwrap();
        end_timer!(commit_t_timer);

        let mut transcript = Transcript::new(b"label");
        let randomness_timer = start_timer!(|| "randomness");
        transcript.append_point::<E>(b"com_t", &CurveGroup::into_affine(com_t));
        transcript.append_scalars(b"com_t", prover.running_accumulator.0.to_sponge_field_elements().as_slice());
        transcript.append_scalars(b"com_t", prover.current_accumulator.0.to_sponge_field_elements().as_slice());
        let _ = transcript.challenge_scalar(b"challenge");
        end_timer!(randomness_timer);

        end_timer!(prover_timer);

        println!("accumulator size: {}", prover.current_accumulator.compressed_size());

        // ********************************** Benchmarking The Decider **********************************
        // decider cost: decide both the running accumulator and the current accumulator
        let verifier_timer = start_timer!(|| "verifier");
        prover.shape.is_relaxed_satisfied(
            &prover.running_accumulator.0,
            &prover.running_accumulator.1,
            &prover.commitment_pp,
        ).unwrap();
        prover.shape.is_satisfied(
            &prover.current_accumulator.0,
            &prover.current_accumulator.1,
            &prover.commitment_pp,
        ).unwrap();
        end_timer!(verifier_timer);
    }
}

#![allow(non_snake_case)]
use ark_ec::pairing::Pairing;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_relations::r1cs::{ConstraintSystem, SynthesisMode};
use ark_serialize::CanonicalSerialize;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use sqrtn_pcs::accumulation_circuit::prover::AccumulatorVerifierCircuitProver;
use sqrtn_pcs::accumulation_circuit::verifier_circuit::AccumulatorVerifierVar;
use sqrtn_pcs::augmented_circuit::augmented_circuit::AugmentedCircuitVar;
use sqrtn_pcs::constant_for_curves::{ScalarField as F, C2, E, G1, G2};
use sqrtn_pcs::kzh::kzh2::{KZH2, KZH2SRS};
use sqrtn_pcs::kzh::KZH;
use sqrtn_pcs::kzh_fold::kzh2_fold::Accumulator2 as Accumulator;
use sqrtn_pcs::nexus_spartan::commitment_traits::ToAffine;
use sqrtn_pcs::nexus_spartan::committed_relaxed_snark::CRSNARKKey;
use sqrtn_pcs::nexus_spartan::crr1cs::{is_sat, produce_synthetic_crr1cs, CRR1CSInstance, CRR1CSShape, CRR1CSWitness};
use sqrtn_pcs::nexus_spartan::crr1csproof::CRR1CSProof;
use sqrtn_pcs::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationAccVerifier, MatrixEvaluationAccVerifierVar};
use sqrtn_pcs::nexus_spartan::partial_verifier::partial_verifier::SpartanPartialVerifier;
use sqrtn_pcs::nexus_spartan::partial_verifier::partial_verifier_var::SpartanPartialVerifierVar;
use sqrtn_pcs::nova::cycle_fold::coprocessor::setup_shape;
use sqrtn_pcs::transcript::transcript::Transcript;
use sqrtn_pcs::transcript::transcript_var::TranscriptVar;

fn bench_augmented_circuit(c: &mut Criterion) {
    let poseidon_iterations_vec = [
        0, 100, 1000, 2000
    ];



    for poseidon_iterations in poseidon_iterations_vec {
        let (pcs_srs, spartan_shape, spartan_instance, spartan_proof, rx, ry) = {
            let num_vars = 131072;
            let num_cons = num_vars;
            let num_inputs = 10;

            // this generates a new instance/witness for spartan as well as PCS parameters
            let (spartan_shape, spartan_instance, spartan_witness, spartan_key) = produce_synthetic_crr1cs::<E, KZH2<E>>(num_cons, num_vars, num_inputs);

            assert!(is_sat(&spartan_shape, &spartan_instance, &spartan_witness, &spartan_key.gens_r1cs_sat).unwrap());

            let pcs_srs = spartan_key.gens_r1cs_sat.clone();

            let mut prover_transcript = Transcript::new(b"example");

            // Get `proof_i` and random evaluation point (r_x, r_y)
            let (spartan_proof, rx, ry) = CRR1CSProof::prove(
                &spartan_shape,
                &spartan_instance,
                spartan_witness,
                &spartan_key.gens_r1cs_sat,
                &mut prover_transcript,
            );

            (pcs_srs, spartan_shape, spartan_instance, spartan_proof, rx, ry)
        };

        // fresh transcripts to be used by the prover and verifier
        let mut prover_transcript = Transcript::new(b"example");
        let verifier_transcript_clone = prover_transcript.clone();
        let cs = ConstraintSystem::<F>::new_ref();

        let partial_verifier_var = {
            let mut verifier_transcript = prover_transcript.clone();
            // Get A(r_x, r_y), B(r_x, r_y), C(r_x, r_y)
            let current_A_B_C_evaluations = spartan_shape.inst.inst.evaluate(&rx, &ry);

            let partial_verifier = SpartanPartialVerifier::initialise(
                &spartan_proof,
                spartan_shape.get_num_vars(),
                spartan_shape.get_num_cons(),
                (spartan_instance.input.assignment, {
                    let com_w: <E as Pairing>::G1Affine = spartan_instance.comm_W.clone().to_affine();
                    com_w
                }),
                &current_A_B_C_evaluations,
                &mut prover_transcript,
            );

            partial_verifier.verify(&mut verifier_transcript);

            let partial_verifier_var = SpartanPartialVerifierVar::new_variable(
                cs.clone(),
                || Ok(partial_verifier.clone()),
                AllocationMode::Input,
            ).unwrap();

            partial_verifier_var
        };

        let acc_verifier_var = {
            let acc_srs = Accumulator::setup(pcs_srs.clone(), &mut thread_rng());

            // Get the KZH opening proof from the Spartan proof
            let opening_proof = spartan_proof.proof_eval_vars_at_ry.clone();

            // Commitment to witness polynomial
            let commitment_w = spartan_instance.comm_W.clone();

            // Sanity check: verify the opening proof
            KZH2::verify(
                &pcs_srs,
                &ry[1..],
                &spartan_proof.eval_vars_at_ry,
                &commitment_w,
                &opening_proof,
            );

            let (x, y) = {
                let split_input = KZH2::split_input(&pcs_srs, &ry[1..]);
                let x = split_input[0].clone();
                let y = split_input[1].clone();

                (x, y)
            };

            // Get accumulator from the opening proof
            let acc_instance = Accumulator::new_accumulator_instance_from_fresh_kzh_instance(
                &acc_srs,
                &commitment_w.C,
                x.as_slice(),
                y.as_slice(),
                &spartan_proof.eval_vars_at_ry,
            );

            let acc_witness = Accumulator::new_accumulator_witness_from_fresh_kzh_witness(
                &acc_srs,
                opening_proof,
                x.as_slice(),
                y.as_slice(),
            );

            let current_acc = Accumulator::new(&acc_instance, &acc_witness);

            // println!("proof size: {}", proof.compressed_size());
            println!("acc size: {}", current_acc.compressed_size());

            // Check that the accumulator is valid
            assert!(
                Accumulator::decide(
                    &acc_srs,
                    &current_acc,
                )
            );

            // use a random accumulator as the running one
            let running_acc = Accumulator::rand(&acc_srs, &mut thread_rng());

            // the shape of the R1CS instance
            let ova_shape = setup_shape::<G1, G2>().unwrap();

            // get trivial running instance
            let (ova_running_instance, ova_running_witness) = AccumulatorVerifierCircuitProver::<G1, G2, C2, E, F>::get_trivial_cycle_fold_running_instance_witness(&ova_shape);

            // get commitment_pp
            let ova_commitment_pp = AccumulatorVerifierCircuitProver::<G1, G2, C2, E, F>::get_commitment_pp(&ova_shape);

            let kzh_acc_verifier_prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, F> = AccumulatorVerifierCircuitProver::new(
                &acc_srs,
                ova_commitment_pp,
                running_acc,
                current_acc.clone(),
                ova_running_instance,
                ova_running_witness,
                prover_transcript,
            );

            // assert it's formated correctly
            kzh_acc_verifier_prover.is_satisfied();

            let acc_verifier_var = AccumulatorVerifierVar::<G1, G2, C2>::new::<E>(cs.clone(), kzh_acc_verifier_prover);

            acc_verifier_var
        };

        let matrix_evaluation_verifier_var = {
            let matrix_eval_acc_verifier = MatrixEvaluationAccVerifier::random_from_eval_point(
                &spartan_shape,
                rx.clone(),
                ry.clone(),
                &mut thread_rng(),
            );

            let matrix_evaluation_verifier_var = MatrixEvaluationAccVerifierVar::new_variable(
                cs.clone(),
                || Ok(matrix_eval_acc_verifier.clone()),
                AllocationMode::Input,
            ).unwrap();

            matrix_evaluation_verifier_var
        };

        // construct the augmented circuit
        let augmented_circuit = AugmentedCircuitVar {
            spartan_partial_verifier: partial_verifier_var,
            kzh_acc_verifier: acc_verifier_var,
            matrix_evaluation_verifier: matrix_evaluation_verifier_var,
        };

        let mut transcript_var = TranscriptVar::from_transcript(cs.clone(), verifier_transcript_clone);

        // run the verification function on augmented circuit
        let _ = augmented_circuit.verify::<E>(cs.clone(), &mut transcript_var, poseidon_iterations);

        assert!(cs.is_satisfied().unwrap());
        println!("augmented circuit constraints: {}", cs.num_constraints());

        // Set the mode to Prove before we convert it for spartan
        cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
        cs.finalize();

        ////////// Now run the spartan prover on the augmented circuit /////////////////

        // convert to the corresponding Spartan types
        let shape = CRR1CSShape::<F>::convert::<G1>(cs.clone());

        // get the number the minimum size we need for committing to the constraint system
        let min_num_vars = CRSNARKKey::<E, KZH2<E>>::get_min_num_vars(shape.get_num_cons(), shape.get_num_vars(), shape.get_num_inputs());
        let SRS: KZH2SRS<E> = KZH2::setup(min_num_vars + 1, &mut thread_rng());


        let bench_name = format!("spartan+commitment time: number of poseidon calls {}", poseidon_iterations);
        c.bench_function(&bench_name, |b| {
            let witness = CRR1CSWitness::<F>::convert(cs.clone());
            let mut new_prover_transcript = Transcript::new(b"example");

            b.iter(|| {
                // committing to the witness
                let instance: CRR1CSInstance<E, KZH2<E>> = CRR1CSInstance::convert(cs.clone(), &SRS);
                // spartan prover: sumchecks and stuff
                let _ = CRR1CSProof::prove(
                    &shape,
                    &instance,
                    witness.clone(),
                    &SRS,
                    &mut new_prover_transcript,
                );
            })
        });

        let bench_name = format!("ABC evals: number of poseidon calls {}", poseidon_iterations);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = shape.inst.inst.evaluate(&rx, &ry);
            })
        });

        let bench_name = format!("verify: number of poseidon calls {}", poseidon_iterations);
        c.bench_function(&bench_name, |b| {
            let instance: CRR1CSInstance<E, KZH2<E>> = CRR1CSInstance::convert(cs.clone(), &SRS);
            let witness = CRR1CSWitness::<F>::convert(cs.clone());

            let mut new_prover_transcript = Transcript::new(b"example");
            let (proof, rx, ry) = CRR1CSProof::prove(
                &shape,
                &instance,
                witness,
                &SRS,
                &mut new_prover_transcript,
            );

            // evaluate matrices A B C
            let inst_evals = shape.inst.inst.evaluate(&rx, &ry);

            b.iter(|| {
                let mut new_verifier_transcript = Transcript::new(b"example");

                assert!(proof
                    .verify(
                        shape.get_num_vars(),
                        shape.get_num_cons(),
                        &instance,
                        &inst_evals,
                        &mut new_verifier_transcript,
                    )
                    .is_ok());
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = augmented_circuit_benches;
    config = custom_criterion_config();
    targets = bench_augmented_circuit
}

criterion_main!(augmented_circuit_benches);





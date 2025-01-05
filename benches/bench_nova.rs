use ark_ec::CurveGroup;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisMode};
use ark_serialize::CanonicalSerialize;
use ark_std::UniformRand;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use sqrtn_pcs::constant_for_curves::{C1, C2, E, G1, G2, ScalarField as F};
use sqrtn_pcs::gadgets::r1cs::r1cs::commit_T;
use sqrtn_pcs::hash::poseidon::PoseidonHashVar;
use sqrtn_pcs::nova::nova::prover::NovaProver;
use sqrtn_pcs::nova::nova::verifier_circuit::NovaVerifierCircuit;
use sqrtn_pcs::nova::nova::verifier_circuit_var::NovaVerifierCircuitVar;
use sqrtn_pcs::transcript::transcript::Transcript;
use sqrtn_pcs::commitment::{CommitmentScheme};

fn bench_nova_prover(c: &mut Criterion) {
    let poseidon_values = [
        0, 10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
        1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000
    ];

    for poseidon_num in poseidon_values {
        // generate some random prover which for a small R1CS which is supposed to generate some sample NovaVerifierCircuit
        // this prover is sort of faking the previous step of IVC, but since we don't have it I just use a prover for a small R1CS
        let prover: NovaProver<F, G1, G2, C1, C2> = NovaProver::rand((10, 3, 17));

        // generate a constraint system
        let cs = ConstraintSystem::<F>::new_ref();

        // generate the Nova verifier circuit
        let nova_verifier_var: NovaVerifierCircuitVar<F, G1, G2, C1, C2> = {
            // the non-circuit version
            let nova_verifier: NovaVerifierCircuit<F, G1, G2, C1, C2> = NovaVerifierCircuit::initialise(prover.clone());
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
        for _ in 0..poseidon_num {
            // get a random element
            let r = FpVar::new_variable(cs.clone(), || Ok(F::rand(&mut thread_rng())), AllocationMode::Witness).unwrap();
            // update sponge with this random element
            hash.update_sponge(vec![r]);
            // output the hash
            let _ = hash.output();
        }

        println!("Number of constraints of the augmented circuit: {}", cs.num_constraints());

        // sum of these two numbers is the accumulator size for Nova
        println!("current accumulator size: {}", prover.current_accumulator.compressed_size());
        println!("running accumulator size: {}", prover.running_accumulator.compressed_size());

        // the cs has to be finalised before the matrices A, B and C are generated
        cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
        cs.finalize();

        // now generate a prover for the augmented circuit
        let prover: NovaProver<F, G1, G2, C1, C2> = NovaProver::rand_from_constraint_system(cs.clone());

        let bench_name = format!("prover time: number of poseidon calls {}", poseidon_num);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                // commit to the current accumulator
                let _ = C1::commit(&prover.commitment_pp, prover.current_accumulator.1.W.as_slice());

                // compute and commit to the cross term error
                let (_, com_t) = commit_T(
                    &prover.shape,
                    &prover.commitment_pp,
                    &prover.running_accumulator.0,
                    &prover.running_accumulator.1,
                    &prover.current_accumulator.0,
                    &prover.current_accumulator.1
                ).unwrap();

                // compute random combination that is to be used for linear combination
                let mut transcript = Transcript::new(b"label");
                transcript.append_point::<E>(b"com_t", &CurveGroup::into_affine(com_t));
                transcript.append_scalars(b"com_t", prover.running_accumulator.0.to_sponge_field_elements().as_slice());
                transcript.append_scalars(b"com_t", prover.current_accumulator.0.to_sponge_field_elements().as_slice());
                let _ = transcript.challenge_scalar(b"challenge");
            })
        });


        let bench_name = format!("verifier time: number of poseidon calls {}", poseidon_num);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
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
            })
        });

    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(20)
}

// Benchmark group setup
criterion_group! {
    name = nova_benches;
    config = custom_criterion_config();
    targets = bench_nova_prover
}

criterion_main!(nova_benches);

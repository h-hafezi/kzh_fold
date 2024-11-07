#![allow(non_snake_case)]
#![allow(unused_imports)]

extern crate criterion;

use criterion::{black_box, Criterion, criterion_group, criterion_main};
use ark_bn254::G1Projective;
use ark_ec::short_weierstrass::{Affine, Projective};
use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef, SynthesisMode};
use ark_std::UniformRand;
use rand::thread_rng;
use sqrtn_pcs::accumulation::accumulator::Accumulator;
use sqrtn_pcs::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use sqrtn_pcs::accumulation_circuit::prover::{get_random_prover, AccumulatorVerifierCircuitProver};
use sqrtn_pcs::accumulation_circuit::verifier_circuit::AccumulatorVerifierVar;
use sqrtn_pcs::commitment::CommitmentScheme;
use sqrtn_pcs::constant_for_curves::{BaseField, E, G1, G2, ScalarField};
use sqrtn_pcs::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use sqrtn_pcs::gadgets::non_native::util::cast_field;
use sqrtn_pcs::hash::pederson::PedersenCommitment;
use sqrtn_pcs::pcs::multilinear_pcs::{PCSEngine, PolynomialCommitmentSRS};
use sqrtn_pcs::transcript::transcript_var::TranscriptVar;

type C2 = PedersenCommitment<Projective<G2>>;

fn setup_benchmark() -> Vec<ScalarField> {
    // a constraint system
    let cs = ConstraintSystem::<ScalarField>::new_ref();

    // initialise the accumulate verifier circuit
    let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E, ScalarField> = get_random_prover();
    let verifier = AccumulatorVerifierVar::<G1, G2, C2>::new::<E>(
        cs.clone(),
        prover.clone()
    );

    println!("number of constraint for initialisation: {}", cs.num_constraints());

    let mut transcript_var = TranscriptVar::from_transcript(
        cs.clone(),
        prover.initial_transcript.clone()
    );

    // run the accumulation
    let _ = verifier.accumulate(&mut transcript_var);

    println!("number of constraint after accumulation: {}", cs.num_constraints());

    // assert the constraint system is satisfied
    assert!(cs.is_satisfied().unwrap());

    // these are required to called CRR1CSShape::convert
    cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
    cs.finalize();

    // write the witness into a file
    let cs_borrow = cs.borrow().unwrap();
    let witness = cs_borrow.witness_assignment.clone();

    witness
}

fn bench_committing_witness(c: &mut Criterion) {
    let witness = setup_benchmark();
    // Setup witness and Pedersen commitment params
    let pp: Vec<Affine<G1>> = PedersenCommitment::<Projective<G1>>::setup(witness.len(), b"test", &());

    // Criterion benchmark for PedersenCommitment::commit
    c.bench_function("PedersenCommitment::commit", |b| {
        b.iter(|| {
            let _c: Projective<G1> = PedersenCommitment::commit(black_box(&pp), black_box(&witness));
        })
    });
}

criterion_group!(benches, bench_committing_witness);
criterion_main!(benches);

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
use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
use rand::thread_rng;
use sqrtn_pcs::accumulation::accumulator::Accumulator;
use sqrtn_pcs::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use sqrtn_pcs::accumulation_circuit::prover::AccumulatorVerifierCircuitProver;
use sqrtn_pcs::accumulation_circuit::verifier_circuit::AccumulatorVerifierVar;
use sqrtn_pcs::commitment::CommitmentScheme;
use sqrtn_pcs::constant_for_curves::{BaseField, E, G1, G2, ScalarField};
use sqrtn_pcs::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use sqrtn_pcs::gadgets::non_native::util::convert_field_one_to_field_two;
use sqrtn_pcs::hash::pederson::PedersenCommitment;
use sqrtn_pcs::nova::cycle_fold::coprocessor_constraints::{R1CSInstanceVar, RelaxedR1CSInstanceVar};
use sqrtn_pcs::pcs::multilinear_pcs::{PolyCommit, PolyCommitTrait, SRS};

type C2 = PedersenCommitment<Projective<G2>>;

pub fn randomness_different_formats(cs: ConstraintSystemRef<ScalarField>, beta: ScalarField) -> (
    BaseField,
    FpVar<ScalarField>,
    NonNativeFieldVar<BaseField, ScalarField>
) {
    let beta_base = convert_field_one_to_field_two::<ScalarField, BaseField>(beta);
    let beta_var = FpVar::new_variable(
        ns!(cs, "beta var"),
        || Ok(beta.clone()),
        AllocationMode::Witness,
    ).unwrap();
    let beta_var_non_native = NonNativeFieldVar::new_variable(
        ns!(cs, "beta var non-native"),
        || Ok(beta_base.clone()),
        AllocationMode::Witness,
    ).unwrap();
    (beta_base, beta_var, beta_var_non_native)
}

fn setup_benchmark() -> Vec<ScalarField> {
    // specifying degrees of polynomials
    let (n, m) = (4, 4);

    // get a random srs
    let srs = {
        let srs_pcs: SRS<E> = PolyCommit::<E>::setup(n, m, &mut thread_rng());
        Accumulator::setup(srs_pcs.clone(), &mut thread_rng())
    };

    // a constraint system
    let cs = ConstraintSystem::<ScalarField>::new_ref();

    let verifier = AccumulatorVerifierVar::rand(&srs, cs.clone());

    println!("number of constraint for initialisation: {}", cs.num_constraints());
    verifier.accumulate();
    println!("number of constraint for initialisation: {}", cs.num_constraints());
    assert!(cs.is_satisfied().unwrap());

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

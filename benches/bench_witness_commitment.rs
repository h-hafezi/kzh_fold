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
use sqrtn_pcs::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use sqrtn_pcs::accumulation_circuit::prover::tests::get_random_prover;
use sqrtn_pcs::accumulation_circuit::verifier_circuit::AccumulatorVerifierVar;
use sqrtn_pcs::commitment::CommitmentScheme;
use sqrtn_pcs::constant_for_curves::{BaseField, G1, ScalarField};
use sqrtn_pcs::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use sqrtn_pcs::gadgets::non_native::util::convert_field_one_to_field_two;
use sqrtn_pcs::hash::pederson::PedersenCommitment;
use sqrtn_pcs::nova::cycle_fold::coprocessor_constraints::{R1CSInstanceVar, RelaxedR1CSInstanceVar};

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
    // a constraint system
    let cs = ConstraintSystem::<ScalarField>::new_ref();

    // get a random initialized prover
    let prover = get_random_prover();

    // the randomness in different formats
    let beta_scalar = prover.beta.clone();
    let (_, beta_var, beta_var_non_native) = randomness_different_formats(cs.clone(), beta_scalar);


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

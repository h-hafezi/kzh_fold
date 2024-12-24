#![allow(non_snake_case)]
#![allow(unused_imports)]

use criterion::{Criterion, criterion_group, criterion_main};
use rand::thread_rng;

use sqrtn_pcs::kzh_fold::kzh2_fold::{Acc2SRS, Accumulator2};
use sqrtn_pcs::constant_for_curves::E;
use sqrtn_pcs::kzh::KZH;
use sqrtn_pcs::kzh::kzh2::{KZH2, KZH2SRS};
use sqrtn_pcs::math::Math;
use sqrtn_pcs::transcript::transcript::Transcript;

fn get_srs(degree: usize) -> Acc2SRS<E> {
    let srs_pcs: KZH2SRS<E> = KZH2::setup(degree, &mut thread_rng());

    // return the result
    Accumulator2::setup(srs_pcs.clone(), &mut thread_rng())
}


fn bench_setup(c: &mut Criterion) {
    // Witness size: 2^2, 2^4 ..., 2^20
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for degree in num_variables {
        let bench_name = format!("setup for degrees n={}", degree);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = get_srs(degree);
            });
        });
    }
}

fn bench_prove(c: &mut Criterion) {
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for degree in num_variables {
        let srs = get_srs(degree);
        let acc_1 = Accumulator2::rand(&srs, &mut thread_rng());
        let acc_2 = Accumulator2::rand(&srs, &mut thread_rng());
        let bench_name = format!("prove for degree n={}", degree);
        let mut transcript = Transcript::new(b"some label");
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator2::prove(&srs, &acc_1, &acc_2, &mut transcript);
            })
        });
    }
}

fn bench_verify(c: &mut Criterion) {
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for degree in num_variables {
        let srs = get_srs(degree);
        let acc_1 = Accumulator2::rand(&srs, &mut thread_rng());
        let acc_2 = Accumulator2::rand(&srs, &mut thread_rng());
        let mut prover_transcript = Transcript::new(b"some label");
        let mut verifier_transcript = prover_transcript.clone();
        let (_, _, Q) = Accumulator2::prove(&srs, &acc_1, &acc_2, &mut prover_transcript);
        let bench_name = format!("verify for degrees n={}", degree);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator2::verify(&srs, &acc_1.instance, &acc_2.instance, Q, &mut verifier_transcript);
            })
        });
    }
}

fn bench_decide(c: &mut Criterion) {
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for degree in num_variables {
        let srs = get_srs(degree);
        let acc = Accumulator2::rand(&srs, &mut thread_rng());
        let bench_name = format!("decide for degrees n={}", degree);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator2::decide(&srs, &acc);
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = kzh2_fold_benches;
    config = custom_criterion_config();
    targets =  bench_verify, bench_prove, bench_decide, bench_setup
}

criterion_main!(kzh2_fold_benches);

#![allow(non_snake_case)]
#![allow(unused_imports)]

use criterion::{Criterion, criterion_group, criterion_main};
use rand::thread_rng;

use sqrtn_pcs::accumulation::accumulator::{AccSRS, Accumulator};
use sqrtn_pcs::constant_for_curves::E;
use sqrtn_pcs::kzh::kzh2::{KZH2, KZH2SRS};
use sqrtn_pcs::transcript::transcript::Transcript;

fn get_srs(degree_x: usize, degree_y: usize) -> AccSRS<E> {
    let srs_pcs: KZH2SRS<E> = KZH2::setup_1(degree_x, degree_y, &mut thread_rng());

    // return the result
    Accumulator::setup(srs_pcs.clone(), &mut thread_rng())
}


fn bench_setup(c: &mut Criterion) {
    // Witness size: 2^2, 2^4 ..., 2^20
    let degrees = vec![(2, 2), (4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let bench_name = format!("setup for degrees n={} * m={} (witness size: {})", degree_x, degree_y, degree_x*degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = get_srs(degree_x, degree_y);
            });
        });
    }
}

fn bench_prove(c: &mut Criterion) {
    let degrees = vec![(2, 2), (4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let srs = get_srs(degree_x, degree_y);
        let acc_1 = Accumulator::rand(&srs, &mut thread_rng());
        let acc_2 = Accumulator::rand(&srs, &mut thread_rng());
        let bench_name = format!("prove for degrees n={} * m={} (witness size: {})", degree_x, degree_y, degree_x*degree_y);
        let mut transcript = Transcript::new(b"some label");
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator::prove(&srs, &acc_1, &acc_2, &mut transcript);
            })
        });
    }
}

fn bench_verify(c: &mut Criterion) {
    let degrees = vec![(2, 2), (4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let srs = get_srs(degree_x, degree_y);
        let acc_1 = Accumulator::rand(&srs, &mut thread_rng());
        let acc_2 = Accumulator::rand(&srs, &mut thread_rng());
        let mut prover_transcript = Transcript::new(b"some label");
        let mut verifier_transcript = prover_transcript.clone();
        let (_, _, Q) = Accumulator::prove(&srs, &acc_1, &acc_2, &mut prover_transcript);
        let bench_name = format!("verify for degrees n={} * m={} (witness size: {})", degree_x, degree_y, degree_x*degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator::verify(&srs, &acc_1.instance, &acc_2.instance, Q, &mut verifier_transcript);
            })
        });
    }
}

fn bench_decide(c: &mut Criterion) {
    let degrees = vec![(2, 2), (4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let srs = get_srs(degree_x, degree_y);
        let acc = Accumulator::rand(&srs, &mut thread_rng());
        let bench_name = format!("decide for degrees n={} * m={} (witness size: {})", degree_x, degree_y, degree_x*degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator::decide(&srs, &acc);
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = acc_benches;
    config = custom_criterion_config();
    targets =  bench_verify, bench_prove, bench_decide, bench_setup
}

criterion_main!(acc_benches);

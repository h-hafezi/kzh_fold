use ark_serialize::CanonicalSerialize;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use sqrtn_pcs::constant_for_curves::E;
use sqrtn_pcs::kzh::KZH;
use sqrtn_pcs::kzh::kzh3::{KZH3, KZH3SRS};
use sqrtn_pcs::kzh_fold::kzh3_fold::{Acc3SRS, Accumulator3};
use sqrtn_pcs::transcript::transcript::Transcript;

fn get_srs(degree: usize) -> Acc3SRS<E> {
    let srs_pcs: KZH3SRS<E> = KZH3::setup(degree, &mut thread_rng());

    // return the result
    Accumulator3::setup(srs_pcs.clone(), &mut thread_rng())
}

fn bench_prove(c: &mut Criterion) {
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for degree in num_variables {
        let srs = get_srs(degree);
        let acc_1 = Accumulator3::rand(&srs);
        let acc_2 = Accumulator3::rand(&srs);
        let bench_name = format!("prove for degree n={}", degree);
        let mut transcript = Transcript::new(b"some label");
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator3::prove(&srs, &acc_1, &acc_2, &mut transcript);
            })
        });
    }
}


fn bench_verify(c: &mut Criterion) {
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for degree in num_variables {
        let srs = get_srs(degree);
        let acc_1 = Accumulator3::rand(&srs);
        let acc_2 = Accumulator3::rand(&srs);
        let mut prover_transcript = Transcript::new(b"some label");
        let mut verifier_transcript = prover_transcript.clone();
        let (_, _, proof) = Accumulator3::prove(&srs, &acc_1, &acc_2, &mut prover_transcript);
        let bench_name = format!("verify for degrees n={}", degree);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator3::verify(&srs, &acc_1.instance, &acc_2.instance, &proof, &mut verifier_transcript);
            })
        });
    }
}

fn bench_decide(c: &mut Criterion) {
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for degree in num_variables {
        let srs = get_srs(degree);
        let acc = Accumulator3::rand(&srs);
        let bench_name = format!("decide for degrees n={}", degree);

        println!("kzh3-fold accumulator length in bytes: {} for degree {degree}", acc.compressed_size());

        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = Accumulator3::decide(&srs, &acc);
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = kzh3_fold_benches;
    config = custom_criterion_config();
    targets =  bench_verify, bench_prove, bench_decide
}

criterion_main!(kzh3_fold_benches);

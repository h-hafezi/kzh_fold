#![allow(non_snake_case)]
#![allow(unused_imports)]

use ark_poly::EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use ark_std::UniformRand;
use criterion::{Criterion, criterion_group, criterion_main};
use rand::thread_rng;

use sqrtn_pcs::constant_for_curves::{E, ScalarField as F};
use sqrtn_pcs::kzh::KZH;
use sqrtn_pcs::kzh::kzh3::{KZH3, KZH3SRS};
use sqrtn_pcs::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

fn bench(c: &mut Criterion) {
    let num_variables = vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    for n in num_variables {
        // get srs
        let srs: KZH3SRS<E> = KZH3::setup(n, &mut thread_rng());

        // random bivariate polynomial
        let polynomial = MultilinearPolynomial::rand(n, &mut thread_rng());

        let bench_name = format!("kzh3 commit for num_variables n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                KZH3::commit(&srs, &polynomial)
            })
        });

        // commit to the polynomial
        let com = KZH3::commit(&srs, &polynomial);

        // open the commitment
        let input: Vec<_> = std::iter::repeat_with(|| F::rand(&mut thread_rng()))
            .take(n)
            .collect();

        let bench_name = format!("kzh3 opening for num_variables n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                KZH3::open(&srs, input.as_slice(), &com, &polynomial);
            })
        });

        let open = KZH3::open(&srs, input.as_slice(), &com, &polynomial);
        println!("witness length in bytes: {} for degree {n}", open.compressed_size());

        let z = polynomial.evaluate(&input);

        let bench_name = format!("kzh3 verifying for num_variables n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                KZH3::verify(&srs, input.as_slice(), &z, &com, &open);
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = kzh3_benches;
    config = custom_criterion_config();
    targets =  bench
}

criterion_main!(kzh3_benches);

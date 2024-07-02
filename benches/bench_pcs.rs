#![allow(non_snake_case)]
#![allow(unused_imports)]
use criterion::{criterion_group, criterion_main, Criterion};
use criterion::BenchmarkId;
use ark_bn254::{Bn254};

use sqrtn_pcs::pcs::{Commitment, OpeningProof, CoeffFormPCS, SRS};
use sqrtn_pcs::{bivariate_poly::BivariatePolynomial};

fn bench_commit(c: &mut Criterion) {
    let rng = &mut rand::thread_rng();

    for logsize in 1..=10 {
        let n = 1<< logsize;
        let srs = CoeffFormPCS::<Bn254>::generate_srs_for_testing(n, rng);
        let bivariate_poly = BivariatePolynomial::random(rng, n);

        let mut group = c.benchmark_group("PCS commit binary");

        group.bench_with_input(BenchmarkId::from_parameter(n), &logsize, |b, _| {
            b.iter(|| CoeffFormPCS::<Bn254>::commit(&bivariate_poly, &srs))
        });
    }
}

fn bench_open(_c: &mut Criterion) {
    unimplemented!();
}

fn bench_verify(_c: &mut Criterion) {
    unimplemented!();
}


criterion_group! {
    name=pcs_benches;
    config=Criterion::default();
    targets=bench_commit,bench_open,bench_verify
}
criterion_main!(pcs_benches);

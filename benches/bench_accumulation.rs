#![allow(non_snake_case)]
#![allow(unused_imports)]
use criterion::{criterion_group, criterion_main, Criterion};
use criterion::BenchmarkId;
use ark_bn254::{Bn254};

use sqrtn_pcs::pcs::{Commitment, OpeningProof, CoeffFormPCS, SRS};
use sqrtn_pcs::{bivariate_poly::BivariatePolynomial};

fn bench_prover(_c: &mut Criterion) {
    unimplemented!();
}

fn bench_verifier(_c: &mut Criterion) {
    unimplemented!();
}


criterion_group! {
    name=accumulation_benches;
    config=Criterion::default();
    targets=bench_prover,bench_verifier
}
criterion_main!(accumulation_benches);

#![allow(non_snake_case)]
#![allow(unused_imports)]

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::UniformRand;
use ark_bn254::{Bn254, Fr};
use criterion::{Criterion, criterion_group, criterion_main};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand::thread_rng;

use sqrtn_pcs::halo_infinite::private_aggregation::{prove, verify};
use sqrtn_pcs::halo_infinite::private_aggregation::tests::{prepare_polynomials_and_srs};
use sqrtn_pcs::constant_for_curves::{E, ScalarField};
use sqrtn_pcs::transcript::transcript::Transcript;

type F = ScalarField;

fn bench_acc_prover(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0u64);
    let N = 2;
    // Gets too slow after 2^15...
    let powers_of_two = (1..=20).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &d in &powers_of_two {
        let (vec_f, _vec_f_commitments, vec_omega, domain, ck, _vk) = prepare_polynomials_and_srs(N, d, &mut rng);

        let bench_name = format!("prove for degree {}", d);
        c.bench_function(&bench_name, |b| {
            let mut transcript_prover = Transcript::new(b"some label");
            b.iter(|| {
                prove::<E, F>(&vec_f, &vec_omega, &domain, &mut transcript_prover, &ck)
            })
        });
    }
}

fn bench_acc_verifier(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0u64);
    let N = 2;
    let powers_of_two = (1..=20).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &d in &powers_of_two {
        let (vec_f, vec_f_commitments, vec_omega, domain, ck, vk) = prepare_polynomials_and_srs(N, d, &mut rng);

        let mut transcript_prover = Transcript::new(b"some label");
        let proof = prove::<E, F>(&vec_f, &vec_omega, &domain, &mut transcript_prover, &ck);

        let bench_name = format!("verify for degree {}", d);
        c.bench_function(&bench_name, |b| {
            let mut transcript_verifier = Transcript::new(b"some label");
            b.iter(|| {
                verify::<E, F>(&proof, &vec_f_commitments, &vec_omega, &domain, &mut transcript_verifier, &vk)
            })
        });
    }
}

criterion_group!(halo_infinite_benches,
                 bench_acc_verifier, bench_acc_prover);
criterion_main!(halo_infinite_benches);

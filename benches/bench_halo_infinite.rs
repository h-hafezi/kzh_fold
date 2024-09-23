#![allow(non_snake_case)]
#![allow(unused_imports)]

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::UniformRand;
use ark_bn254::{Bn254, Fr};
use criterion::{Criterion, criterion_group, criterion_main};
use rand::rngs::StdRng;
use transcript::IOPTranscript;
use rand::SeedableRng;
use rand::thread_rng;

use sqrtn_pcs::halo_infinite::private_aggregation::{prove, verify};
use sqrtn_pcs::halo_infinite::private_aggregation::tests::{prepare_polynomials_and_srs};
use sqrtn_pcs::constant_for_curves::{E, ScalarField};

fn bench_acc_prover(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0u64);
    let N = 2;
    let powers_of_two = (1..=20).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &d in &powers_of_two {
        let (vec_f, _vec_f_commitments, vec_omega, domain, ck, _vk) = prepare_polynomials_and_srs(N, d, &mut rng);

        let bench_name = format!("prove for degree {}", d);
        c.bench_function(&bench_name, |b| {
            let mut transcript_prover = IOPTranscript::<Fr>::new(b"pa");
            b.iter(|| {
                prove::<Bn254>(&vec_f, &vec_omega, &domain, &mut transcript_prover, &ck)
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

        let mut transcript_prover = IOPTranscript::<Fr>::new(b"pa");
        let proof = prove::<Bn254>(&vec_f, &vec_omega, &domain, &mut transcript_prover, &ck);

        let bench_name = format!("verify for degree {}", d);
        c.bench_function(&bench_name, |b| {
            let mut transcript_verifier = IOPTranscript::<Fr>::new(b"pa");
            b.iter(|| {
                verify::<Bn254>(&proof, &vec_f_commitments, &vec_omega, &domain, &mut transcript_verifier, &vk)
            })
        });
    }
}

criterion_group!(halo_infinite_benches, bench_acc_prover, bench_acc_verifier);
criterion_main!(halo_infinite_benches);

#![allow(non_snake_case)]
#![allow(unused_imports)]

use ark_bn254::{Bn254, Fr};
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ec::AffineRepr;
use ark_std::UniformRand;
use criterion::{Criterion, criterion_group, criterion_main, black_box};
use rand::{Rng, thread_rng};
use sqrtn_pcs::constant_for_curves::{E, ScalarField, G1, G1Affine, G1Projective};

fn bench_msm(c: &mut Criterion) {
    let mut rng = thread_rng();
    let sizes = (18..=24).map(|i| 1 << i).collect::<Vec<_>>();

    for &size in &sizes {
        // Generate random bases
        let bases = (0..size)
            .map(|_| G1Affine::rand(&mut rng))
            .collect::<Vec<_>>();

        // Generate scalars that are 4-bit unsigned integers (in [0..15])
        let scalars = (0..size)
            .map(|_| {
                ScalarField::from(4 as u64)
            })
            .collect::<Vec<_>>();

        let bench_name = format!("MSM for size {}", size);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                // Perform the MSM
                let _ =
                    G1Projective::msm_unchecked(bases.as_slice(), scalars.as_slice());
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = msm_benches;
    config = custom_criterion_config();
    targets = bench_msm
}

criterion_main!(msm_benches);

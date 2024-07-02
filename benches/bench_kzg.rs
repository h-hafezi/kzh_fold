#![allow(non_snake_case)]
#![allow(unused_imports)]

use ark_bn254::{Bn254, Fr};
use ark_ec::pairing::Pairing;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_std::UniformRand;
use criterion::{Criterion, criterion_group, criterion_main};
use rand::{Rng, thread_rng};
use sqrtn_pcs::constant_for_curves::{E, ScalarField};
use sqrtn_pcs::kzg::{KZG10, trim, KZGPowers, KZGUniversalParams, KZGVerifierKey};
use sqrtn_pcs::math::Math;

type Poly = DensePolynomial<<E as Pairing>::ScalarField>;


fn bench(c: &mut Criterion) {
    let mut rng = thread_rng();
    let degrees: Vec<usize> = (10..21) // Range from 10 to 20
        .map(|x| 1 << x as usize)
        .collect();

    for degree in degrees {
        let params = KZG10::<E, Poly>::setup(degree, false, &mut rng).expect("Setup failed");

        // trim the public parameters in to ck and vk
        let (ck, vk) = trim(&params, degree);

        // a random polynomial
        let polynomial = Poly::rand(degree, &mut rng);

        let bench_name = format!("commit for number of variables {}", degree.log_2());
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                KZG10::<E, Poly>::commit(&ck, &polynomial, None, Some(&mut rng)).expect("Commitment failed")
            })
        });

        // commit to the polynomial
        let (comm, r) = KZG10::<E, Poly>::commit(&ck, &polynomial, None, Some(&mut rng)).expect("Commitment failed");

        let bench_name = format!("open for number of variables {}", degree.log_2());
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let point = ScalarField::rand(&mut rng);
                KZG10::<E, Poly>::open(&ck, &polynomial, point, &r).expect("Proof generation failed")
            })
        });

        // open the polynomial at a random point
        let point = ScalarField::rand(&mut rng);
        let proof = KZG10::<E, Poly>::open(&ck, &polynomial, point, &r).expect("Proof generation failed");
        let value = polynomial.evaluate(&point);

        let bench_name = format!("verify for number of variables {}", degree.log_2());
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                KZG10::<E, Poly>::check(&vk, &comm, point, value, &proof).expect("Verification failed")
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = kzg_benches;
    config = custom_criterion_config();
    targets = bench
}

criterion_main!(kzg_benches);

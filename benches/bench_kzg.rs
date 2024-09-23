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
use sqrtn_pcs::kzg::{KZG10, KZGPowers, KZGUniversalParams, KZGVerifierKey};

type Poly = DensePolynomial<<E as Pairing>::ScalarField>;

pub(crate) fn trim(pp: &KZGUniversalParams<E>, mut supported_degree: usize) -> (KZGPowers<E>, KZGVerifierKey<E>) {
    if supported_degree == 1 {
        supported_degree += 1;
    }
    let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
    let powers_of_gamma_g = (0..=supported_degree)
        .map(|i| pp.powers_of_gamma_g[&i])
        .collect();

    let powers = KZGPowers {
        powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
        powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
    };
    let vk = KZGVerifierKey {
        g: pp.powers_of_g[0],
        gamma_g: pp.powers_of_gamma_g[&0],
        h: pp.h,
        beta_h: pp.beta_h,
        prepared_h: pp.prepared_h.clone(),
        prepared_beta_h: pp.prepared_beta_h.clone(),
    };
    (powers, vk)
}

fn rand<R: Rng>(d: usize, rng: &mut R) -> Poly {
    let mut random_coeffs = Vec::new();
    for _ in 0..=d {
        random_coeffs.push(ScalarField::rand(rng));
    }
    Poly::from_coefficients_vec(random_coeffs)
}

fn bench_setup(c: &mut Criterion) {
    let mut rng = thread_rng();
    let sqrt_degrees = vec![4, 8, 16, 32, 64, 128, 256, 512, 1024];

    for &sqrt_degree in &sqrt_degrees {
        let degree = sqrt_degree * sqrt_degree;
        let bench_name = format!("setup for degree {}", degree);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                KZG10::<E, Poly>::setup(degree, false, &mut rng).expect("Setup failed")
            })
        });
    }
}

fn bench_commit(c: &mut Criterion) {
    let mut rng = thread_rng();
    let sqrt_degrees = vec![4, 8, 16, 32, 64, 128, 256, 512, 1024];

    for &sqrt_degree in &sqrt_degrees {
        let degree = sqrt_degree * sqrt_degree;
        let params = KZG10::<E, Poly>::setup(degree, false, &mut rng).expect("Setup failed");
        let (ck, _vk) = trim(&params, degree);

        let polynomial = Poly::rand(degree, &mut rng);
        let bench_name = format!("commit for degree {}", degree);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                KZG10::<E, Poly>::commit(&ck, &polynomial, None, Some(&mut rng)).expect("Commitment failed")
            })
        });
    }
}

fn bench_open(c: &mut Criterion) {
    let mut rng = thread_rng();
    let sqrt_degrees = vec![4, 8, 16, 32, 64, 128, 256, 512, 1024];

    for &sqrt_degree in &sqrt_degrees {
        let degree = sqrt_degree * sqrt_degree;
        let params = KZG10::<E, Poly>::setup(degree, false, &mut rng).expect("Setup failed");
        let (ck, _vk) = trim(&params, degree);

        let polynomial = Poly::rand(degree, &mut rng);
        let (_, r) = KZG10::<E, Poly>::commit(&ck, &polynomial, None, Some(&mut rng)).expect("Commitment failed");

        let bench_name = format!("open for degree {}", degree);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let point = ScalarField::rand(&mut rng);
                KZG10::<E, Poly>::open(&ck, &polynomial, point, &r).expect("Proof generation failed")
            })
        });
    }
}

fn bench_verify(c: &mut Criterion) {
    let mut rng = thread_rng();
    let sqrt_degrees = vec![4, 8, 16, 32, 64, 128, 256, 512, 1024];

    for &sqrt_degree in &sqrt_degrees {
        let degree = sqrt_degree * sqrt_degree;
        let params = KZG10::<E, Poly>::setup(degree, false, &mut rng).expect("Setup failed");
        let (ck, vk) = trim(&params, degree);

        let polynomial = Poly::rand(degree, &mut rng);
        let (comm, r) = KZG10::<E, Poly>::commit(&ck, &polynomial, None, Some(&mut rng)).expect("Commitment failed");

        let point = ScalarField::rand(&mut rng);
        let proof = KZG10::<E, Poly>::open(&ck, &polynomial, point, &r).expect("Proof generation failed");
        let value = polynomial.evaluate(&point);

        let bench_name = format!("verify for degree {}", degree);
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
    targets = bench_setup, bench_commit, bench_open, bench_verify
}

criterion_main!(kzg_benches);

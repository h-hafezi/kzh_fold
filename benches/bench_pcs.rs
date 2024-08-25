#![allow(non_snake_case)]
#![allow(unused_imports)]

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::UniformRand;
use criterion::{Criterion, criterion_group, criterion_main};
use rand::thread_rng;
use sqrtn_pcs::constant_for_curves::{E, ScalarField};
use sqrtn_pcs::polynomial::bivariate_poly::{BivariatePolynomial};
use sqrtn_pcs::polynomial::lagrange_basis::LagrangeBasis;
use sqrtn_pcs::polynomial_commitment::bivariate_pcs::{PolyCommit, PolyCommitTrait, SRS};


fn bench_setup(c: &mut Criterion) {
    let mut rng = thread_rng();
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];

    for (degree_x, degree_y) in degrees {
        let bench_name = format!("setup for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter_custom(|iters| {
                let mut total_time = std::time::Duration::new(0, 0);
                let mut size = 0;
                for _ in 0..iters {
                    let start = std::time::Instant::now();
                    let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut rng);
                    total_time += start.elapsed();
                    // Measure the heap-allocated size of the srs object
                    size = srs.size_of();
                }
                // Print the size of the SRS object once at the end of the last iteration
                println!("Size of SRS, degree n={} * m={}: {} bytes", degree_x, degree_y, size);
                total_time
            });
        });
    }
}

fn bench_commit(c: &mut Criterion) {
    let mut rng = thread_rng();
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut rng);
        let poly_commit = PolyCommit { srs };
        let polynomial = BivariatePolynomial::random(&mut rng, domain_x.clone(), domain_y.clone(), degree_x, degree_y);
        let bench_name = format!("commit for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                poly_commit.commit(&polynomial)
            })
        });
    }
}

fn bench_open(c: &mut Criterion) {
    let mut rng = thread_rng();
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];

    for (degree_x, degree_y) in degrees {
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut rng);
        let poly_commit = PolyCommit { srs };
        let polynomial = BivariatePolynomial::random(&mut rng, domain_x.clone(), domain_y.clone(), degree_x, degree_y);
        let com = poly_commit.commit(&polynomial);

        let bench_name = format!("open for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |x| {
            let b = ScalarField::rand(&mut rng);
            x.iter(|| {
                poly_commit.open(&polynomial, com.clone(), &b)
            })
        });
    }
}

fn bench_verify(c: &mut Criterion) {
    let mut rng = thread_rng();
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];

    for (degree_x, degree_y) in degrees {
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut rng);
        let poly_commit = PolyCommit { srs };
        let polynomial = BivariatePolynomial::random(&mut rng, domain_x.clone(), domain_y.clone(), degree_x, degree_y);
        let com = poly_commit.commit(&polynomial);

        let bench_name = format!("verify for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |x| {
            let b = ScalarField::rand(&mut rng);
            let c = ScalarField::rand(&mut rng);
            let y = polynomial.evaluate(&b, &c);
            let open = poly_commit.open(&polynomial, com.clone(), &b);
            x.iter(|| {
                poly_commit.verify(LagrangeBasis { domain: domain_x.clone() }, &com, &open, &b, &c, &y)
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = pcs_benches;
    config = custom_criterion_config();
    targets =  bench_setup, bench_commit, bench_open, bench_verify
}

criterion_main!(pcs_benches);

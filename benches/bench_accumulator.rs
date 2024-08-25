#![allow(non_snake_case)]
#![allow(unused_imports)]

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::UniformRand;
use criterion::{Criterion, criterion_group, criterion_main};
use rand::thread_rng;

use sqrtn_pcs::accumulation::accumulator::{AccSRS, Accumulator};
use sqrtn_pcs::constant_for_curves::{E, ScalarField};
use sqrtn_pcs::polynomial::bivariate_poly::{BivariatePolynomial};
use sqrtn_pcs::polynomial::lagrange_basis::LagrangeBasis;
use sqrtn_pcs::polynomial_commitment::bivariate_pcs::{PolyCommit, PolyCommitTrait, SRS};

fn get_srs(degree_x: usize, degree_y: usize) -> AccSRS<E> {
    let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
    let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();

    // define the srs
    let pc_srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut thread_rng());

    // set accumulator srs
    let lagrange_x = LagrangeBasis { domain: domain_x.clone() };
    let lagrange_y = LagrangeBasis { domain: domain_y.clone() };
    Accumulator::setup(
        degree_x,
        degree_y,
        lagrange_x,
        lagrange_y,
        pc_srs,
        &mut thread_rng(),
    )
}

fn get_satisfying_accumulator(srs: &AccSRS<E>) -> Accumulator<E> {

    // random bivariate polynomials
    let polynomial_1 = BivariatePolynomial::random(
        &mut thread_rng(),
        srs.lagrange_basis_x.domain.clone(),
        srs.lagrange_basis_y.domain.clone(),
        srs.degree_x,
        srs.degree_y,
    );

    let polynomial_2 = BivariatePolynomial::random(
        &mut thread_rng(),
        srs.lagrange_basis_x.domain.clone(),
        srs.lagrange_basis_y.domain.clone(),
        srs.degree_x,
        srs.degree_y,
    );

    // random points and evaluation
    let b_1 = ScalarField::rand(&mut thread_rng());
    let c_1 = ScalarField::rand(&mut thread_rng());
    let y_1 = polynomial_1.evaluate(&b_1, &c_1);
    let b_2 = ScalarField::rand(&mut thread_rng());
    let c_2 = ScalarField::rand(&mut thread_rng());
    let y_2 = polynomial_2.evaluate(&b_2, &c_2);

    // define the polynomial commitment scheme
    let poly_commit = PolyCommit { srs: srs.pc_srs.clone() };

    // commit to the polynomials
    let com_1 = poly_commit.commit(&polynomial_1);
    let com_2 = poly_commit.commit(&polynomial_2);

    // open the commitment
    let open_1 = poly_commit.open(&polynomial_1, com_1.clone(), &b_1);
    let open_2 = poly_commit.open(&polynomial_2, com_2.clone(), &b_2);


    // get accumulator instance/proof from polynomial instance/opening
    let instance_1 = Accumulator::new_accumulator_instance_from_proof(&srs, &com_1.C, &b_1, &c_1, &y_1);
    let witness_1 = Accumulator::new_accumulator_witness_from_proof(&srs, open_1.clone(), &b_1, &c_1);
    let instance_2 = Accumulator::new_accumulator_instance_from_proof(&srs, &com_2.C, &b_2, &c_2, &y_2);
    let witness_2 = Accumulator::new_accumulator_witness_from_proof(&srs, open_2.clone(), &b_2, &c_2);

    // define accumulators
    let acc_1 = Accumulator::new_accumulator(&instance_1, &witness_1);
    let acc_2 = Accumulator::new_accumulator(&instance_2, &witness_2);

    // asserting decide without accumulation
    assert!(Accumulator::decide(&srs, &acc_1));
    assert!(Accumulator::decide(&srs, &acc_2));

    let beta = ScalarField::rand(&mut thread_rng());

    let (acc_instance, acc_witness, _Q) = Accumulator::prove(&srs, &beta, &acc_1, &acc_2);

    return Accumulator {
        witness: acc_witness,
        instance: acc_instance,
    };
}

fn bench_setup(c: &mut Criterion) {
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];

    for (degree_x, degree_y) in degrees {
        let bench_name = format!("setup for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                get_srs(degree_x, degree_y);
            });
        });
    }
}

fn bench_prove(c: &mut Criterion) {
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let srs = get_srs(degree_x, degree_y);
        let acc_1 = get_satisfying_accumulator(&srs);
        let acc_2 = get_satisfying_accumulator(&srs);
        let bench_name = format!("prove for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let beta = ScalarField::rand(&mut thread_rng());
                Accumulator::prove(&srs, &beta, &acc_1, &acc_2);
            })
        });
    }
}

fn bench_verify(c: &mut Criterion) {
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let srs = get_srs(degree_x, degree_y);
        let acc_1 = get_satisfying_accumulator(&srs);
        let acc_2 = get_satisfying_accumulator(&srs);
        let beta = ScalarField::rand(&mut thread_rng());
        let (_, _, Q) = Accumulator::prove(&srs, &beta, &acc_1, &acc_2);
        let bench_name = format!("verify for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                Accumulator::verify(&acc_1.instance, &acc_2.instance, Q, &beta);
            })
        });
    }
}

fn bench_decide(c: &mut Criterion) {
    let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)];
    for (degree_x, degree_y) in degrees {
        let srs = get_srs(degree_x, degree_y);
        let acc = get_satisfying_accumulator(&srs);
        let bench_name = format!("decide for degree n={} * m={}", degree_x, degree_y);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                Accumulator::decide(&srs, &acc);
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = acc_benches;
    config = custom_criterion_config();
    targets =  bench_setup, bench_verify, bench_prove, bench_decide
}

criterion_main!(acc_benches);

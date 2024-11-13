use ark_std::UniformRand;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use sqrtn_pcs::constant_for_curves::{ScalarField, E};
use sqrtn_pcs::math::Math;
use sqrtn_pcs::pcs::kzh3::kzh3::{commit, open, verify, KZH3SRS};
use sqrtn_pcs::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

fn bench_commit(c: &mut Criterion) {
    let degrees = vec![(4, 4, 4), (8, 8, 8)];
    for (degree_x, degree_y, degree_z) in degrees {
        // build the srs
        let srs: KZH3SRS<E> = KZH3SRS::setup(degree_x, degree_y, degree_z, &mut thread_rng());
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        // random  polynomial
        let polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        let bench_name = format!("commit for degrees x={} * y={} * z={} (witness size: {})", degree_x, degree_y, degree_z, degree_x * degree_y * degree_z);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                // commit to the polynomial
                let _ = commit(&srs, &polynomial);
            })
        });
    }
}

fn bench_open(c: &mut Criterion) {
    let degrees = vec![(4, 4, 4), (8, 8, 8)];
    for (degree_x, degree_y, degree_z) in degrees {
        // build the srs
        let srs: KZH3SRS<E> = KZH3SRS::setup(degree_x, degree_y, degree_z, &mut thread_rng());
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        // random  polynomial
        let polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());
        let com = commit(&srs, &polynomial);

        let bench_name = format!("commit for degrees x={} * y={} * z={} (witness size: {})", degree_x, degree_y, degree_z, degree_x * degree_y * degree_z);
        c.bench_function(&bench_name, |b| {
            let x :Vec<ScalarField> = {
                let mut res = Vec::new();
                for _ in 0..degree_x.log_2() {
                    res.push(ScalarField::rand(&mut thread_rng()));
                }
                res
            };

            let y :Vec<ScalarField> = {
                let mut res = Vec::new();
                for _ in 0..degree_y.log_2() {
                    res.push(ScalarField::rand(&mut thread_rng()));
                }
                res
            };

            let z :Vec<ScalarField> = {
                let mut res = Vec::new();
                for _ in 0..degree_z.log_2() {
                    res.push(ScalarField::rand(&mut thread_rng()));
                }
                res
            };

            b.iter(|| {
                open(&srs, &polynomial, x.as_slice(), y.as_slice())
            })
        });
    }
}

fn bench_verify(c: &mut Criterion) {
    let degrees = vec![(4, 4, 4), (8, 8, 8)];
    for (degree_x, degree_y, degree_z) in degrees {
        // build the srs
        let srs: KZH3SRS<E> = KZH3SRS::setup(degree_x, degree_y, degree_z, &mut thread_rng());
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        // random  polynomial
        let polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());
        let com = commit(&srs, &polynomial);

        let bench_name = format!("commit for degrees x={} * y={} * z={} (witness size: {})", degree_x, degree_y, degree_z, degree_x * degree_y * degree_z);
        c.bench_function(&bench_name, |b| {
            let x :Vec<ScalarField> = {
                let mut res = Vec::new();
                for _ in 0..degree_x.log_2() {
                    res.push(ScalarField::rand(&mut thread_rng()));
                }
                res
            };

            let y :Vec<ScalarField> = {
                let mut res = Vec::new();
                for _ in 0..degree_y.log_2() {
                    res.push(ScalarField::rand(&mut thread_rng()));
                }
                res
            };

            let z :Vec<ScalarField> = {
                let mut res = Vec::new();
                for _ in 0..degree_z.log_2() {
                    res.push(ScalarField::rand(&mut thread_rng()));
                }
                res
            };

            let concat = {
                let mut res = vec![];
                res.extend(x.clone());
                res.extend(y.clone());
                res.extend(z.clone());
                res
            };

            let eval = polynomial.evaluate(&concat);

            let open = open(&srs, &polynomial, x.as_slice(), y.as_slice());

            b.iter(|| {
                verify(&srs, &com, x.as_slice(), y.as_slice(), z.as_slice(), &open, eval);
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
    targets =  bench_commit,
}

criterion_main!(kzh3_benches);

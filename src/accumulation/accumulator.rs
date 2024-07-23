use std::ops::{Add, Mul, Sub};
use ark_bn254::G1Affine;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ff::{FftField, Field, Zero};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::UniformRand;
use rand::RngCore;
use crate::lagrange_basis::{LagrangeBasis, LagrangeTraits};
use crate::pcs::{Commitment, OpeningProof, PolyCommit, PolyCommitTrait, SRS};
use crate::univariate_poly::UnivariatePolynomial;
use crate::utils::{inner_product, is_power_of_two, power};


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccSRS<E: Pairing> {
    degree_x: usize,
    degree_y: usize,
    lagrange_basis_x: LagrangeBasis<E::ScalarField>,
    lagrange_basis_y: LagrangeBasis<E::ScalarField>,

    // vector of size degree_x
    k_vec_b: Vec<E::G1Affine>,

    // vector of size degree_y
    k_vec_c: Vec<E::G1Affine>,

    k_prime: E::G1Affine,
    pc_srs: SRS<E>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccInstance<E: Pairing> {
    C: E::G1Affine,
    T: E::G1Affine,
    b: E::ScalarField,
    c: E::ScalarField,
    y: E::ScalarField,
    z_b: E::ScalarField,
    z_c: E::ScalarField,
    E: E::G1Affine,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccWitness<E: Pairing> {
    // size of degree_x
    pub vec_D: Vec<E::G1Affine>,
    pub f_star_poly: UnivariatePolynomial<E::ScalarField>,
    // size of degree_x
    pub vec_b: Vec<E::ScalarField>,
    // size of degree_y
    pub vec_c: Vec<E::ScalarField>,
}


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Accumulator<E: Pairing> {
    pub srs: AccSRS<E>,
    pub witness: AccWitness<E>,
    pub instance: AccInstance<E>,
}

pub trait AccumulatorTrait<E: Pairing> {
    fn setup<T: RngCore>(degree_x: usize,
                         degree_y: usize,
                         lagrange_basis_x: LagrangeBasis<E::ScalarField>,
                         lagrange_basis_y: LagrangeBasis<E::ScalarField>,
                         pc_srs: SRS<E>,
                         rng: &mut T,
    ) -> AccSRS<E>;

    fn new_accumulator(srs: &AccSRS<E>, instance: &AccInstance<E>, witness: &AccWitness<E>) -> Accumulator<E>;

    fn new_accumulator_instance_from_proof(srs: &AccSRS<E>, C: &E::G1Affine, b: &E::ScalarField, c: &E::ScalarField, y: &E::ScalarField) -> AccInstance<E>;

    fn new_accumulator_witness_from_proof(srs: &AccSRS<E>, proof: OpeningProof<E>, b: &E::ScalarField, c: &E::ScalarField) -> AccWitness<E>;

    fn prove(acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> (AccInstance<E>, AccWitness<E>, E::G1Affine);

    fn verify(instance_1: &AccInstance<E>, instance_2: &AccInstance<E>, Q: E::G1Affine) -> AccInstance<E>;

    fn decide(acc: &Accumulator<E>) -> bool;

    /// Compute e_i <- b_i * (b - w_i) - z_b : ∀i ∈ [0, n]
    /// Check that e'_j <- c_j * (c - w_j) - z_c : ∀j ∈ [0, m]
    /// Compute e'' <- ⟨f*, c⟩ - y
    /// Compute E_G <- ⟨f*, (H_1, ..., H_m)⟩ - ⟨b, D⟩
    /// Outputs V(b, c, y, z_b, z_c, f*, b, c, D) = ⟨e' || e'' || K⟩ + E_G
    fn helper_function_V(acc: &Accumulator<E>) -> E::G1Affine;

    fn helper_function_Q(acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> E::G1Affine;
}

impl<E: Pairing> AccumulatorTrait<E> for Accumulator<E> {
    fn setup<T: RngCore>(degree_x: usize,
                         degree_y: usize,
                         lagrange_basis_x: LagrangeBasis<E::ScalarField>,
                         lagrange_basis_y: LagrangeBasis<E::ScalarField>,
                         pc_srs: SRS<E>,
                         rng: &mut T,
    ) -> AccSRS<E> {
        // asserting the degree_x and degree_y are a power of two
        assert!(is_power_of_two(degree_x), "degree_x (upper bound) must be a power of two");
        assert!(is_power_of_two(degree_y), "degree_y (upper bound) must be a power of two");
        // return the result
        return AccSRS {
            degree_x,
            degree_y,
            lagrange_basis_x,
            lagrange_basis_y,
            pc_srs,
            k_vec_b: {
                let mut elements = Vec::new();
                for _ in 0..degree_x {
                    elements.push(E::G1Affine::rand(rng));
                }
                elements
            },
            k_vec_c: {
                let mut elements = Vec::new();
                for _ in 0..degree_y {
                    elements.push(E::G1Affine::rand(rng));
                }
                elements
            },
            // TODO: to change
            k_prime: E::G1Affine::zero(),
        };
    }

    fn new_accumulator(srs: &AccSRS<E>, instance: &AccInstance<E>, witness: &AccWitness<E>) -> Accumulator<E> {
        Accumulator {
            srs: srs.clone(),
            witness: witness.clone(),
            instance: instance.clone(),
        }
    }

    fn new_accumulator_instance_from_proof(srs: &AccSRS<E>, C: &E::G1Affine, b: &E::ScalarField, c: &E::ScalarField, y: &E::ScalarField) -> AccInstance<E> {
        let vec_b = srs.lagrange_basis_x.evaluate(b);
        let vec_c = srs.lagrange_basis_y.evaluate(c);

        // asserting the sizes are correct
        assert_eq!(vec_b.len(), srs.degree_x, "invalid size");
        assert_eq!(vec_c.len(), srs.degree_y, "invalid size");

        // using function power to do optimal field multiplications
        let z_b = power(*b, srs.degree_x) - E::ScalarField::ONE;
        let z_c = power(*c, srs.degree_y) - E::ScalarField::ONE;

        let T = E::G1::msm_unchecked(srs.k_vec_b.as_slice(), vec_b.as_slice()).add(E::G1::msm_unchecked(srs.k_vec_c.as_slice(), vec_c.as_slice()));

        return AccInstance {
            C: *C,
            T: T.into(),
            b: *b,
            c: *c,
            y: *y,
            z_b,
            z_c,
            E: E::G1Affine::zero(),
        };
    }

    fn new_accumulator_witness_from_proof(srs: &AccSRS<E>, proof: OpeningProof<E>, b: &E::ScalarField, c: &E::ScalarField) -> AccWitness<E> {
        let vec_b = srs.lagrange_basis_x.evaluate(b);
        let vec_c = srs.lagrange_basis_y.evaluate(c);

        // asserting the sizes are correct
        assert_eq!(vec_b.len(), srs.degree_x, "invalid size");
        assert_eq!(vec_c.len(), srs.degree_y, "invalid size");
        assert_eq!(proof.vec_D.len(), srs.degree_x, "invalid size");

        return AccWitness {
            vec_D: proof.vec_D,
            f_star_poly: proof.f_star_poly,
            vec_b,
            vec_c,
        };
    }

    fn prove(acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> (AccInstance<E>, AccWitness<E>, E::G1Affine) {
        // making sure the lagrange basis for two f_star_poly are the same
        assert_eq!(acc_1.witness.f_star_poly.lagrange_basis, acc_2.witness.f_star_poly.lagrange_basis, "lagrange basis need to be equal");

        // assert that two accumulators have the same srs
        assert_eq!(acc_1.srs, acc_2.srs, "the accumulators must have the same srs");

        // unwrap the instances and witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        // compute the quotient variable Q
        let Q: E::G1Affine = Self::helper_function_Q(acc_1, acc_2);

        // generate random values beta
        let beta: E::ScalarField = E::ScalarField::from(2u128);
        let beta_minus_one: E::ScalarField = E::ScalarField::ONE - beta;

        // get the accumulated new_instance
        let new_instance = Self::verify(instance_1, instance_2, Q);

        // get the accumulated witness
        let new_witness = AccWitness {
            vec_D: {
                witness_1.vec_D.iter()
                    .zip(witness_2.vec_D.iter())
                    .map(
                        |(&d_1, &d_2)|
                        d_1.mul(beta_minus_one).add(d_2.mul(beta)).into_affine()
                    )
                    .collect()
            },
            f_star_poly: UnivariatePolynomial {
                evaluations: {
                    witness_1.f_star_poly.evaluations.iter()
                        .zip(witness_2.f_star_poly.evaluations.iter())
                        .map(
                            |(&a, &b)|
                            a * (beta_minus_one) + (b * beta)
                        )
                        .collect()
                },
                lagrange_basis: witness_1.f_star_poly.lagrange_basis.clone(),
            },
            vec_b: {
                witness_1.vec_b.iter()
                    .zip(witness_2.vec_b.iter())
                    .map(
                        |(&b_1, &b_2)|
                        b_1 * (beta_minus_one) + (b_2 * beta)
                    )
                    .collect()
            },
            vec_c: {
                witness_1.vec_c.iter()
                    .zip(witness_2.vec_c.iter())
                    .map(
                        |(&c_1, &c_2)|
                        c_1 * (beta_minus_one) + (c_2 * beta)
                    )
                    .collect()
            },
        };
        return (new_instance, new_witness, Q);
    }

    fn verify(instance_1: &AccInstance<E>, instance_2: &AccInstance<E>, Q: E::G1Affine) -> AccInstance<E> {
        let beta = E::ScalarField::from(2u128);
        let beta_minus_one: E::ScalarField = E::ScalarField::ONE - beta;

        let new_error_term: E::G1Affine = {
            let mut res = instance_1.E.mul(beta_minus_one);
            res = res.add(instance_2.E.mul(beta));
            res.add(Q.mul(beta_minus_one * beta)).into()
        };

        AccInstance {
            C: {
                let res = instance_1.C.mul(beta_minus_one);
                res.add(instance_2.C.mul(beta)).into()
            },
            T: {
                let res = instance_1.T.mul(beta_minus_one);
                res.add(instance_2.T.mul(beta)).into()
            },
            b: {
                let res = instance_1.b * beta_minus_one;
                res + (instance_2.b * beta)
            },
            c: {
                let res = instance_1.c * beta_minus_one;
                res + (instance_2.c * beta)
            },
            y: {
                let res = instance_1.y * beta_minus_one;
                res + (instance_2.y * beta)
            },
            z_b: {
                let res = instance_1.z_b * beta_minus_one;
                res + (instance_2.z_b * beta)
            },
            z_c: {
                let res = instance_1.z_c * beta_minus_one;
                res + (instance_2.z_c * beta)
            },
            E: new_error_term,
        }
    }

    fn decide(acc: &Accumulator<E>) -> bool {
        let srs = &acc.srs;
        let instance = &acc.instance;
        let witness = &acc.witness;

        // first condition
        let pairing_lhs = E::multi_pairing(&witness.vec_D, &srs.pc_srs.vec_V);
        let pairing_rhs = E::pairing(instance.C, srs.pc_srs.V_prime);

        // second condition
        let ip_rhs = instance.T;
        let ip_lhs = {
            let res = E::G1::msm_unchecked(srs.k_vec_b.as_slice(), witness.vec_b.as_slice());
            res.add(E::G1::msm_unchecked(srs.k_vec_c.as_slice(), witness.vec_c.as_slice()))
        };

        // third condition
        let verify_lhs = Self::helper_function_V(acc);
        let verify_rhs = instance.E;

        println!("here: {} {} {}", verify_rhs == verify_lhs.into(), ip_lhs == ip_rhs.into(), ip_lhs == ip_rhs.into());
        return (verify_rhs == verify_lhs.into()) && (ip_lhs == ip_rhs.into()) && (pairing_lhs == pairing_rhs);
    }

    fn helper_function_V(acc: &Accumulator<E>) -> E::G1Affine {
        let srs = &acc.srs;
        let instance = &acc.instance;
        let witness = &acc.witness;

        let vec_e_b: Vec<E::ScalarField> = {
            let mut vec_e = Vec::with_capacity(witness.vec_b.len());
            let vec_b = &witness.vec_b;
            let vec_w = srs.lagrange_basis_x.domain.elements();
            for (&b_i, omega_i) in vec_b.iter().zip(vec_w.into_iter()) {
                let constant_i = srs.lagrange_basis_x.domain.size_inv() * omega_i;
                let e_i = ((b_i * (instance.b - omega_i)) / (constant_i)) - instance.z_b;
                vec_e.push(e_i);
            }
            vec_e
        };

        let vec_e_c: Vec<E::ScalarField> = {
            let mut vec_e = Vec::with_capacity(witness.vec_c.len());
            let vec_c = &witness.vec_c;
            let vec_w = srs.lagrange_basis_y.domain.elements();
            for (&c_i, omega_i) in vec_c.iter().zip(vec_w.into_iter()) {
                let constant_i = srs.lagrange_basis_y.domain.size_inv() * omega_i;
                let e_i = ((c_i * (instance.c - omega_i)) / (constant_i)) - instance.z_c;
                vec_e.push(e_i);
            }
            vec_e
        };

        let e_prime: E::ScalarField = inner_product(witness.f_star_poly.evaluations.as_slice(), witness.vec_c.as_slice()) - instance.y;

        let E_G: E::G1 = {
            let lhs = E::G1::msm_unchecked(srs.pc_srs.vec_H.as_slice(), witness.f_star_poly.evaluations.as_slice());
            let rhs = {
                E::G1::msm_unchecked(witness.vec_D.as_slice(), witness.vec_b.as_slice())
            };
            lhs.sub(rhs)
        };

        let mut res: E::G1 = E_G.clone();
        res = res.add(E::G1::msm_unchecked(srs.k_vec_b.as_slice(), vec_e_b.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_vec_c.as_slice(), vec_e_c.as_slice()));
        res.add(srs.k_prime.mul(e_prime)).into()
    }

    fn helper_function_Q(acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> E::G1Affine {
        // unwrap the instances/witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        // value of 2 in the field
        let two = E::ScalarField::from(2u8);

        // build the accumulator from linear combination to run helper_function_V on it
        let temp_acc = Accumulator {
            srs: acc_1.srs.clone(),
            witness: AccWitness {
                vec_D: {
                    witness_2.vec_D.iter().zip(witness_1.vec_D.iter())
                        .map(|(&d2, &d1)| d2.mul(two).sub(d1).into_affine())
                        .collect()
                },
                f_star_poly: UnivariatePolynomial {
                    evaluations: {
                        witness_2.f_star_poly.evaluations.iter().zip(witness_1.f_star_poly.evaluations.iter())
                            .map(|(&e2, &e1)| two * e2 - e1)
                            .collect()
                    },
                    lagrange_basis: witness_1.f_star_poly.lagrange_basis.clone(),
                },
                vec_b: {
                    witness_2.vec_b.iter().zip(witness_1.vec_b.iter())
                        .map(|(&b2, &b1)| two * b2 - b1)
                        .collect()
                },
                vec_c: {
                    witness_2.vec_c.iter().zip(witness_1.vec_c.iter())
                        .map(|(&c2, &c1)| two * c2 - c1)
                        .collect()
                },
            },
            instance: AccInstance {
                // not used by helper function
                C: E::G1Affine::zero(),
                // not used by the helper function
                T: E::G1Affine::zero(),
                b: two * instance_2.b - instance_1.b,
                c: two * instance_2.c - instance_1.c,
                y: two * instance_2.y - instance_1.y,
                z_b: two * instance_2.z_b - instance_1.z_b,
                z_c: two * instance_2.z_c - instance_1.z_c,
                // not used by the helper function
                E: E::G1Affine::zero(),
            },
        };

        // 1/2 in the scalar field
        let one_over_two = two.inverse().unwrap();

        let mut res = Self::helper_function_V(&temp_acc);
        res = res.mul(one_over_two).into();
        res = res.add(instance_1.E).into();
        res = res.sub(instance_2.E).into();
        res.sub(instance_2.E).into()
    }
}

#[cfg(test)]
mod tests {
    use ark_bn254::{Bn254, Fr, Fq, G1Affine, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::accumulation::accumulator::{Accumulator, AccumulatorTrait};
    use crate::bivariate_poly::{BivariatePolynomial, BivariatePolynomialTrait};
    use crate::lagrange_basis::LagrangeBasis;
    use crate::pcs::{Commitment, OpeningProof, PolyCommit, PolyCommitTrait, SRS};

    type E = Bn254;
    type G1 = G1Affine;
    type G2 = G2Affine;
    type ScalarField = Fr;

    fn get_poly_parameters(degree_x: usize, degree_y: usize) -> (SRS<E>, ScalarField, ScalarField, ScalarField, Commitment<E>, OpeningProof<E>, LagrangeBasis<ScalarField>, LagrangeBasis<ScalarField>) {
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut thread_rng());
        // define the polynomial commitment
        let poly_commit = PolyCommit { srs: srs.clone() };
        // random bivariate polynomial
        let polynomial = BivariatePolynomial::random(&mut thread_rng(), domain_x, domain_y, degree_x, degree_y);
        // random points and evaluation
        let b = ScalarField::rand(&mut thread_rng());
        let c = ScalarField::rand(&mut thread_rng());
        let y = polynomial.evaluate(&b, &c);
        // commit to the polynomial
        let com = poly_commit.commit(&polynomial);
        // open the commitment
        let open = poly_commit.open(&polynomial, com.clone(), &b);
        // assert correctness of pcs
        assert!(poly_commit.verify(LagrangeBasis { domain: domain_x }, &com, &open, &b, &c, &y));
        // return
        return (srs, b, c, y, com, open, LagrangeBasis { domain: domain_x }, LagrangeBasis { domain: domain_y });
    }
    #[test]
    fn raw_decide_test() {
        let degree_x = 4usize;
        let degree_y = 16usize;
        // get a polynomial and its related fields
        let (pc_srs, b, c, y, com, open, lagrange_x, lagrange_y) = get_poly_parameters(degree_x, degree_y);
        // setup the srs
        let srs = Accumulator::setup(degree_x, degree_y, lagrange_x, lagrange_y, pc_srs, &mut thread_rng());
        // get accumulator instance/proof from polynomial instance/opening
        let instance = Accumulator::new_accumulator_instance_from_proof(&srs, &com.C, &b, &c, &y);
        let witness = Accumulator::new_accumulator_witness_from_proof(&srs, open.clone(), &b, &c);
        // assign a new accumulator instance
        let acc = Accumulator::new_accumulator(&srs, &instance, &witness);
        // decide it
        assert!(Accumulator::decide(&acc));
    }

    #[test]
    fn raw_decide_test_1() {
        let degree_x = 4usize;
        let degree_y = 16usize;
        // get a polynomial and its related fields
        let (pc_srs, b, c, y, com, open, lagrange_x, lagrange_y) = get_poly_parameters(degree_x, degree_y);
        // setup the srs
        let srs = Accumulator::setup(degree_x, degree_y, lagrange_x, lagrange_y, pc_srs, &mut thread_rng());
        // get accumulator instance/proof from polynomial instance/opening
        let instance_1 = Accumulator::new_accumulator_instance_from_proof(&srs, &com.C, &b, &c, &y);
        let witness_1 = Accumulator::new_accumulator_witness_from_proof(&srs, open.clone(), &b, &c);
        // assign a new accumulator instance
        let acc_1 = Accumulator::new_accumulator(&srs, &instance_1, &witness_1);
        let acc_2 = Accumulator::new_accumulator(&srs, &instance_1, &witness_1);
        let (acc_instance, acc_witness, Q) = Accumulator::prove(&acc_1, &acc_2);
        let acc_instance_prime = Accumulator::verify(&instance_1, &instance_1, Q);
        // asserting instances are equal
        assert_eq!(acc_instance, acc_instance_prime);
        // deciding the accumulator
        let acc = Accumulator::new_accumulator(&srs, &acc_instance, &acc_witness);
        assert!(Accumulator::decide(&acc));
    }

    #[test]
    fn accumulation_test_2() {
        // set polynomial degree
        let degree_x = 4usize;
        let degree_y = 16usize;
        // define unity roots for degree_x and degree_y
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        // define the srs
        let pc_srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut thread_rng());
        // define the polynomial commitment scheme
        let poly_commit = PolyCommit { srs: pc_srs.clone() };
        // random bivariate polynomials
        let polynomial_1 = BivariatePolynomial::random(&mut thread_rng(), domain_x, domain_y, degree_x, degree_y);
        let polynomial_2 = BivariatePolynomial::random(&mut thread_rng(), domain_x, domain_y, degree_x, degree_y);
        // random points and evaluation
        let b_1 = ScalarField::rand(&mut thread_rng());
        let c_1 = ScalarField::rand(&mut thread_rng());
        let y_1 = polynomial_1.evaluate(&b_1, &c_1);
        let b_2 = ScalarField::rand(&mut thread_rng());
        let c_2 = ScalarField::rand(&mut thread_rng());
        let y_2 = polynomial_2.evaluate(&b_2, &c_2);
        // commit to the polynomials
        let com_1 = poly_commit.commit(&polynomial_1);
        let com_2 = poly_commit.commit(&polynomial_2);
        // open the commitment
        let open_1 = poly_commit.open(&polynomial_1, com_1.clone(), &b_1);
        let open_2 = poly_commit.open(&polynomial_2, com_2.clone(), &b_2);
        // set accumulator srs
        let lagrange_x = LagrangeBasis { domain: domain_x };
        let lagrange_y = LagrangeBasis { domain: domain_y };
        let srs = Accumulator::setup(degree_x, degree_y, lagrange_x, lagrange_y, pc_srs, &mut thread_rng());
        // get accumulator instance/proof from polynomial instance/opening
        let instance_1 = Accumulator::new_accumulator_instance_from_proof(&srs, &com_1.C, &b_1, &c_1, &y_1);
        let witness_1 = Accumulator::new_accumulator_witness_from_proof(&srs, open_1.clone(), &b_1, &c_1);
        let instance_2 = Accumulator::new_accumulator_instance_from_proof(&srs, &com_2.C, &b_2, &c_2, &y_2);
        let witness_2 = Accumulator::new_accumulator_witness_from_proof(&srs, open_2.clone(), &b_2, &c_2);
        // define accumulators
        let acc_1 = Accumulator::new_accumulator(&srs, &instance_1, &witness_1);
        let acc_2 = Accumulator::new_accumulator(&srs, &instance_2, &witness_2);
        // asserting decide without accumulation
        //assert!(Accumulator::decide(&acc_1));
        //assert!(Accumulator::decide(&acc_2));
        // accumulate proof
        let (acc_instance, acc_witness, Q) = Accumulator::prove(&acc_1, &acc_2);
        let acc_instance_prime = Accumulator::verify(&instance_1, &instance_2, Q);
        // asserting instances are equal
        assert_eq!(acc_instance, acc_instance_prime);
        // deciding the accumulator
        let acc = Accumulator::new_accumulator(&srs, &acc_instance, &acc_witness);
        assert!(Accumulator::decide(&acc));
    }
}

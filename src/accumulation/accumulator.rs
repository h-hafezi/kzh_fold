use std::ops::{Add, Mul, Neg, Sub};

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ff::{AdditiveGroup, FftField, Field, PrimeField, Zero};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::UniformRand;
use rand::{Rng, RngCore};

use crate::gadgets::non_native::util::convert_field_one_to_field_two;
use crate::hash::poseidon::{PoseidonHash, PoseidonHashTrait};
use crate::polynomial::lagrange_basis::{LagrangeBasis};
use crate::polynomial::univariate_poly::UnivariatePolynomial;
use crate::polynomial_commitment::pcs::{OpeningProof, PolyCommit, PolyCommitTrait, SRS};
use crate::utils::{inner_product, is_power_of_two, power};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccSRS<E: Pairing> {
    pub degree_x: usize,
    pub degree_y: usize,
    pub lagrange_basis_x: LagrangeBasis<E::ScalarField>,
    pub lagrange_basis_y: LagrangeBasis<E::ScalarField>,

    // vector of size degree_x
    pub k_vec_b: Vec<E::G1Affine>,

    // vector of size degree_y
    pub k_vec_c: Vec<E::G1Affine>,

    pub k_prime: E::G1Affine,
    pub pc_srs: SRS<E>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccInstance<E: Pairing> {
    pub C: E::G1Affine,
    pub T: E::G1Affine,
    pub E: E::G1Affine,
    pub b: E::ScalarField,
    pub c: E::ScalarField,
    pub y: E::ScalarField,
    pub z_b: E::ScalarField,
    pub z_c: E::ScalarField,
}

impl<E: Pairing> AccInstance<E>
where
    E::ScalarField: PrimeField,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    pub fn to_sponge_field_elements(&self) -> Vec<E::ScalarField> {
        let mut dest = Vec::new();
        // Define a closure to handle the conversion of affine points to scalars
        let convert_affine_to_scalars = |point: E::G1Affine| {
            if point.is_zero() {
                (E::ScalarField::ONE, E::ScalarField::ZERO)
            } else {
                // Extract x and y coordinates and convert them
                let x = convert_field_one_to_field_two::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(point.x().unwrap());
                let y = convert_field_one_to_field_two::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(point.y().unwrap());
                (x, y)
            }
        };

        // Use the closure for C, T, and E
        let (c_x, c_y) = convert_affine_to_scalars(self.C);
        let (t_x, t_y) = convert_affine_to_scalars(self.T);
        let (e_x, e_y) = convert_affine_to_scalars(self.E);

        // Extend the destination vector with the computed values
        dest.extend(vec![c_x, c_y, t_x, t_y, e_x, e_y]);

        // Extend with other scalar fields
        dest.extend(vec![self.b, self.c, self.y, self.z_b, self.z_c]);

        dest
    }
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
    pub witness: AccWitness<E>,
    pub instance: AccInstance<E>,
}

impl<E: Pairing> Accumulator<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
{
    pub fn setup<T: RngCore>(degree_x: usize,
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
            k_prime: E::G1Affine::rand(rng),
        };
    }

    pub fn new_accumulator(instance: &AccInstance<E>, witness: &AccWitness<E>) -> Accumulator<E> {
        Accumulator {
            witness: witness.clone(),
            instance: instance.clone(),
        }
    }

    pub fn compute_randomness(instance_1: &AccInstance<E>, instance_2: &AccInstance<E>, Q: E::G1Affine) -> E::ScalarField {
        let mut sponge = Vec::new();
        sponge.extend(instance_1.to_sponge_field_elements());
        sponge.extend(instance_2.to_sponge_field_elements());

        // Define a closure to handle the conversion of affine points to scalars
        let convert_affine_to_scalars = |point: E::G1Affine| {
            if point.is_zero() {
                vec![E::ScalarField::ONE, E::ScalarField::ZERO]
            } else {
                // Extract x and y coordinates and convert them
                let x = convert_field_one_to_field_two::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(point.x().unwrap());
                let y = convert_field_one_to_field_two::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(point.y().unwrap());
                vec![x, y]
            }
        };
        sponge.extend(convert_affine_to_scalars(Q));

        let mut hash_object = PoseidonHash::new();
        hash_object.update_sponge(sponge);
        hash_object.output()
    }

    pub fn new_accumulator_instance_from_proof(srs: &AccSRS<E>, C: &E::G1Affine, b: &E::ScalarField, c: &E::ScalarField, y: &E::ScalarField) -> AccInstance<E> {
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

    pub fn new_accumulator_witness_from_proof(srs: &AccSRS<E>, proof: OpeningProof<E>, b: &E::ScalarField, c: &E::ScalarField) -> AccWitness<E> {
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

    pub fn prove(srs: &AccSRS<E>, acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> (Accumulator<E>, E::G1Affine)
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
    {
        // making sure the lagrange basis for two f_star_poly are the same
        assert_eq!(acc_1.witness.f_star_poly.lagrange_basis, acc_2.witness.f_star_poly.lagrange_basis, "lagrange basis need to be equal");

        // unwrap the instances and witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        // compute the quotient variable Q
        let Q: E::G1Affine = Self::helper_function_Q(srs, acc_1, acc_2);

        let beta = Accumulator::compute_randomness(instance_1, instance_2, Q);

        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        // get the accumulated new_instance
        let new_instance = Self::verify(instance_1, instance_2, Q);

        // get the accumulated witness
        let new_witness = AccWitness {
            vec_D: {
                witness_1.vec_D.iter()
                    .zip(witness_2.vec_D.iter())
                    .map(
                        |(&d_1, &d_2)|
                        d_1.mul(one_minus_beta).add(d_2.mul(beta)).into_affine()
                    )
                    .collect()
            },
            f_star_poly: UnivariatePolynomial {
                evaluations: {
                    witness_1.f_star_poly.evaluations.iter()
                        .zip(witness_2.f_star_poly.evaluations.iter())
                        .map(
                            |(&a, &b)|
                            a * (one_minus_beta) + (b * beta)
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
                        b_1 * (one_minus_beta) + (b_2 * beta)
                    )
                    .collect()
            },
            vec_c: {
                witness_1.vec_c.iter()
                    .zip(witness_2.vec_c.iter())
                    .map(
                        |(&c_1, &c_2)|
                        c_1 * (one_minus_beta) + (c_2 * beta)
                    )
                    .collect()
            },
        };
        let acc_prime = Accumulator {
            witness: new_witness,
            instance: new_instance,
        };
        (acc_prime, Q)
    }

    pub fn verify(instance_1: &AccInstance<E>, instance_2: &AccInstance<E>, Q: E::G1Affine) -> AccInstance<E> {
        let beta = Accumulator::compute_randomness(instance_1, instance_2, Q);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        let new_error_term: E::G1Affine = {
            let mut res = instance_1.E.mul(one_minus_beta);
            res = res.add(instance_2.E.mul(beta));
            res.add(Q.mul(one_minus_beta * beta)).into()
        };

        AccInstance {
            C: {
                let res = instance_1.C.mul(one_minus_beta);
                res.add(instance_2.C.mul(beta)).into()
            },
            T: {
                let res = instance_1.T.mul(one_minus_beta);
                res.add(instance_2.T.mul(beta)).into()
            },
            b: {
                let res = instance_1.b * one_minus_beta;
                res + (instance_2.b * beta)
            },
            c: {
                let res = instance_1.c * one_minus_beta;
                res + (instance_2.c * beta)
            },
            y: {
                let res = instance_1.y * one_minus_beta;
                res + (instance_2.y * beta)
            },
            z_b: {
                let res = instance_1.z_b * one_minus_beta;
                res + (instance_2.z_b * beta)
            },
            z_c: {
                let res = instance_1.z_c * one_minus_beta;
                res + (instance_2.z_c * beta)
            },
            E: new_error_term,
        }
    }

    pub fn decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> bool {
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
        let verify_lhs = Self::helper_function_decide(srs, acc);
        let verify_rhs = instance.E;

        return (verify_rhs == verify_lhs.into()) && (ip_lhs == ip_rhs.into()) && (pairing_lhs == pairing_rhs);
    }

    pub fn helper_function_decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> E::G1Affine {
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

    pub fn helper_function_Q(srs: &AccSRS<E>, acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> E::G1Affine {
        // unwrap the instances/witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        // value of 2 in the field
        let two = E::ScalarField::from(2u128);

        // build the accumulator from linear combination to run helper_function_V on it
        let temp_acc = Accumulator {
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

        // -1/2 in the scalar field
        let one_over_two: E::ScalarField = two.neg().inverse().unwrap();


        let mut res = Self::helper_function_decide(&srs, &temp_acc);
        res = res.add(instance_1.E).into();
        res = res.sub(instance_2.E).into();
        res = res.sub(instance_2.E).into();
        res.mul(one_over_two).into()
    }
}

pub fn get_srs<E, T: Rng>(degree_x: usize, degree_y: usize, rng: &mut T) -> AccSRS<E>
where
    E: Pairing,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    <E as Pairing>::ScalarField: Absorb,
{
    let domain_x = GeneralEvaluationDomain::<E::ScalarField>::new(degree_x).unwrap();
    let domain_y = GeneralEvaluationDomain::<E::ScalarField>::new(degree_y).unwrap();

    // define the srs
    let pc_srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, rng);

    // set accumulator srs
    let lagrange_x = LagrangeBasis { domain: domain_x.clone() };
    let lagrange_y = LagrangeBasis { domain: domain_y.clone() };
    Accumulator::setup(
        degree_x,
        degree_y,
        lagrange_x,
        lagrange_y,
        pc_srs,
        rng,
    )
}

#[cfg(test)]
pub mod tests {
    use ark_ec::AffineRepr;
    use ark_ec::pairing::Pairing;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::UniformRand;
    use rand::thread_rng;

    use super::*;
    use crate::constant_for_curves::{E, ScalarField};
    use crate::polynomial::bivariate_poly::{BivariatePolynomial};
    use crate::polynomial::lagrange_basis::LagrangeBasis;
    use crate::polynomial_commitment::pcs::{Commitment, OpeningProof, PolyCommit, PolyCommitTrait, SRS};


    pub fn get_satisfying_accumulator(srs: &AccSRS<E>) -> Accumulator<E>
    {

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

        let (accumulator, _Q) = Accumulator::prove(&srs, &acc_1, &acc_2);

        accumulator
    }


    fn get_poly_parameters(
        degree_x: usize,
        degree_y: usize,
    ) -> (SRS<E>,
          ScalarField,
          ScalarField,
          ScalarField,
          Commitment<E>,
          OpeningProof<E>,
          LagrangeBasis<ScalarField>,
          LagrangeBasis<ScalarField>
    ) {
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut thread_rng());

        // define the polynomial commitment
        let poly_commit = PolyCommit { srs: srs.clone() };

        // random bivariate polynomial
        let polynomial = BivariatePolynomial::random(&mut thread_rng(), domain_x.clone(), domain_y.clone(), degree_x, degree_y);

        // random points and evaluation
        let b = ScalarField::rand(&mut thread_rng());
        let c = ScalarField::rand(&mut thread_rng());
        let y = polynomial.evaluate(&b, &c);

        // commit to the polynomial
        let com = poly_commit.commit(&polynomial);

        // open the commitment
        let open = poly_commit.open(&polynomial, com.clone(), &b);

        // assert correctness of pcs
        assert!(poly_commit.verify(LagrangeBasis { domain: domain_x.clone() }, &com, &open, &b, &c, &y));

        // return
        return (srs, b, c, y, com, open, LagrangeBasis { domain: domain_x.clone() }, LagrangeBasis { domain: domain_y.clone() });
    }

    #[test]
    fn raw_decide_test() {
        let degree_x = 4usize;
        let degree_y = 16usize;

        // get a polynomial and its related fields
        let (pc_srs,
            b,
            c,
            y,
            com,
            open,
            lagrange_x,
            lagrange_y
        ) = get_poly_parameters(degree_x, degree_y);

        // setup the srs
        let srs = Accumulator::setup(degree_x, degree_y, lagrange_x, lagrange_y, pc_srs, &mut thread_rng());

        // get accumulator instance/proof from polynomial instance/opening
        let instance = Accumulator::new_accumulator_instance_from_proof(&srs, &com.C, &b, &c, &y);
        let witness = Accumulator::new_accumulator_witness_from_proof(&srs, open.clone(), &b, &c);

        // assign a new accumulator instance
        let acc = Accumulator::new_accumulator(&instance, &witness);

        // decide it
        assert!(Accumulator::decide(&srs, &acc));
    }


    #[test]
    fn general_accumulation_test() {
        let degree_x = 16;
        let degree_y = 16;
        let srs = get_srs(degree_x, degree_y, &mut thread_rng());
        let acc_1 = get_satisfying_accumulator(&srs);
        let acc_2 = get_satisfying_accumulator(&srs);

        // asserting decide without accumulation
        assert!(Accumulator::decide(&srs, &acc_1));
        assert!(Accumulator::decide(&srs, &acc_2));

        let beta = ScalarField::rand(&mut thread_rng());

        // accumulate proof
        let (accumulator, Q) = Accumulator::prove(&srs, &acc_1, &acc_2);

        // accumulate verifier
        let instance_prime = Accumulator::verify(&acc_1.instance, &acc_2.instance, Q);

        // deciding the accumulator
        assert!(Accumulator::decide(&srs, &accumulator));

        assert_eq!(accumulator.instance, instance_prime);
    }
}

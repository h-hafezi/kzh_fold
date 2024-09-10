use std::ops::{Add, Mul, Neg, Sub};

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ff::{AdditiveGroup, Field, PrimeField, Zero};
use ark_poly::EvaluationDomain;
use ark_std::UniformRand;
use rand::RngCore;
use crate::polynomial::multilinear_polynomial::{compute_dot_product};
use crate::accumulation::{convert_affine_to_scalars, generate_random_elements};
use crate::accumulation::eq_tree::EqTree;
use crate::hash::poseidon::{PoseidonHash, PoseidonHashTrait};
use crate::pcs::multilinear_pcs::{OpeningProof, SRS};
use crate::polynomial::multilinear_polynomial::multilinear_poly::MultilinearPolynomial;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccSRS<E: Pairing> {
    // vector of size 2 * degree_x - 1
    pub k_x: Vec<E::G1Affine>,

    // vector of size 2 * degree_y - 1
    pub k_y: Vec<E::G1Affine>,

    pub k_prime: E::G1Affine,
    pub pc_srs: SRS<E>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccInstance<E: Pairing> {
    pub C: E::G1Affine,
    pub T: E::G1Affine,
    pub E: E::G1Affine,

    // vector of length log2(degree_x)
    pub x: Vec<E::ScalarField>,

    // vector of length log2(degree_y)
    pub y: Vec<E::ScalarField>,

    pub z: E::ScalarField,
}

impl<E: Pairing> AccInstance<E>
where
    E::ScalarField: PrimeField,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    pub fn to_sponge_field_elements(&self) -> Vec<E::ScalarField> {
        let mut dest = Vec::new();
        // Define a closure to handle the conversion of affine points to scalars

        // Use the closure for C, T, and E
        let (c_x, c_y) = convert_affine_to_scalars::<E>(self.C);
        let (t_x, t_y) = convert_affine_to_scalars::<E>(self.T);
        let (e_x, e_y) = convert_affine_to_scalars::<E>(self.E);

        // Extend the destination vector with the computed values
        dest.extend(vec![c_x, c_y, t_x, t_y, e_x, e_y]);

        // Extend with other scalar fields
        dest.extend(vec![self.z]);
        dest.extend(self.x.clone());
        dest.extend(self.y.clone());

        dest
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccWitness<E: Pairing> {
    // size of degree_x
    pub vec_D: Vec<E::G1Affine>,

    // MLP of degree_y
    pub f_star_poly: MultilinearPolynomial<E::ScalarField, E>,

    pub tree_x: EqTree<E::ScalarField>,

    pub tree_y: EqTree<E::ScalarField>,
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
    pub fn setup<T: RngCore>(pc_srs: SRS<E>, rng: &mut T) -> AccSRS<E> {
        // return the result
        AccSRS {
            pc_srs: pc_srs.clone(),
            k_x: generate_random_elements::<E, T>(2 * pc_srs.degree_x - 1, rng),
            k_y: generate_random_elements::<E, T>(2 * pc_srs.degree_y - 1, rng),
            k_prime: E::G1Affine::rand(rng),
        }
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
        let (p1, p2) = convert_affine_to_scalars::<E>(Q);
        sponge.extend(vec![p1, p2]);

        let mut hash_object = PoseidonHash::new();
        hash_object.update_sponge(sponge);
        hash_object.output()
    }

    pub fn new_accumulator_instance_from_proof(
        srs: &AccSRS<E>,
        C: &E::G1Affine,
        x: &Vec<E::ScalarField>,
        y: &Vec<E::ScalarField>,
        z: &E::ScalarField,
    ) -> AccInstance<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.pc_srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.pc_srs.degree_y, "invalid size of vector y");

        let tree_x = EqTree::new(x.as_slice());
        let tree_y = EqTree::new(y.as_slice());

        let mut T: E::G1 = E::G1::ZERO;
        T = T.add(E::G1::msm_unchecked(srs.k_x.as_slice(), tree_x.nodes.as_slice()));
        T = T.add(E::G1::msm_unchecked(srs.k_y.as_slice(), tree_y.nodes.as_slice()));

        assert_eq!(srs.k_x.len(), tree_x.nodes.len(), "invalid size of vector x");
        assert_eq!(srs.k_y.len(), tree_y.nodes.len(), "invalid size of vector y");


        AccInstance {
            C: *C,
            T: T.into(),
            E: E::G1Affine::zero(),
            x: x.clone(),
            y: y.clone(),
            z: z.clone(),
        }
    }

    pub fn new_accumulator_witness_from_proof(srs: &AccSRS<E>, proof: OpeningProof<E>, x: &Vec<E::ScalarField>, y: &Vec<E::ScalarField>) -> AccWitness<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.pc_srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.pc_srs.degree_y, "invalid size of vector y");
        assert_eq!(proof.vec_D.len(), srs.pc_srs.degree_x, "invalid proof size");

        AccWitness {
            vec_D: proof.vec_D,
            f_star_poly: proof.f_star_poly,
            tree_x: EqTree::new(x.as_slice()),
            tree_y: EqTree::new(y.as_slice()),
        }
    }

    pub fn prove(srs: &AccSRS<E>, acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> (AccInstance<E>, AccWitness<E>, E::G1Affine)
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
    {
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
            vec_D: witness_1.vec_D.iter().zip(witness_2.vec_D.iter())
                .map(|(&d_1, &d_2)| d_1.mul(one_minus_beta).add(d_2.mul(beta)).into_affine())
                .collect(),
            f_star_poly: MultilinearPolynomial {
                num_variables: witness_1.f_star_poly.num_variables,
                evaluation_over_boolean_hypercube: witness_1.f_star_poly.evaluations_over_boolean_domain().iter()
                    .zip(witness_2.f_star_poly.evaluations_over_boolean_domain().iter())
                    .map(
                        |(&a, &b)|
                        a * (one_minus_beta) + (b * beta)
                    )
                    .collect(),
                phantom: Default::default(),
            },
            tree_x: EqTree {
                nodes: witness_1.tree_x.nodes.iter()
                    .zip(witness_2.tree_x.nodes.iter())
                    .map(
                        |(&a, &b)|
                        a * (one_minus_beta) + (b * beta)
                    )
                    .collect(),
                depth: witness_1.tree_x.depth,
            },
            tree_y: EqTree {
                nodes: witness_1.tree_y.nodes.iter()
                    .zip(witness_2.tree_y.nodes.iter())
                    .map(
                        |(&a, &b)|
                        a * (one_minus_beta) + (b * beta)
                    )
                    .collect(),
                depth: witness_1.tree_y.depth,
            },
        };
        return (new_instance, new_witness, Q);
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
            x: instance_1.x.iter()
                .zip(instance_2.x.iter())
                .map(|(&e1, &e2)| e1 * one_minus_beta + e1 * beta)
                .collect(),
            y: instance_1.y.iter()
                .zip(instance_2.y.iter())
                .map(|(&e1, &e2)| e1 * one_minus_beta + e1 * beta)
                .collect(),
            z: instance_1.z * one_minus_beta + instance_2.z * beta,
            E: new_error_term,
        }
    }

    pub fn helper_function_Q(srs: &AccSRS<E>, acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> E::G1Affine {
        // unwrap the instances/witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        assert_eq!(witness_1.f_star_poly.num_variables, witness_2.f_star_poly.num_variables);
        assert_eq!(witness_1.tree_x.depth, witness_2.tree_x.depth);
        assert_eq!(witness_1.tree_y.depth, witness_2.tree_y.depth);
        assert_eq!(witness_1.vec_D.len(), witness_2.vec_D.len());

        // value of 2 in the field
        let two = E::ScalarField::from(2u128);

        // build the accumulator from linear combination to run helper_function_V on it
        let temp_acc = Accumulator {
            witness: AccWitness {
                vec_D: witness_2.vec_D.iter()
                    .zip(witness_1.vec_D.iter())
                    .map(|(&d2, &d1)| d2.mul(two).sub(d1).into_affine())
                    .collect(),
                f_star_poly: MultilinearPolynomial {
                    num_variables: witness_1.f_star_poly.num_variables,
                    evaluation_over_boolean_hypercube: witness_2.f_star_poly.evaluation_over_boolean_hypercube.iter()
                        .zip(witness_1.f_star_poly.evaluation_over_boolean_hypercube.iter())
                        .map(|(&e2, &e1)| e2 * two - e1)
                        .collect(),
                    phantom: Default::default(),
                },
                tree_x: EqTree {
                    nodes: witness_2.tree_x.nodes.iter()
                        .zip(witness_1.tree_x.nodes.iter())
                        .map(|(&w2, &w1)| w2 * two - w1)
                        .collect(),
                    depth: witness_1.tree_x.depth,
                },
                tree_y: EqTree {
                    nodes: witness_2.tree_y.nodes.iter()
                        .zip(witness_1.tree_y.nodes.iter())
                        .map(|(&w2, &w1)| w2 * two - w1)
                        .collect(),
                    depth: witness_1.tree_y.depth,
                },
            },
            instance: AccInstance {
                // not used by helper function
                C: E::G1Affine::zero(),
                T: E::G1Affine::zero(),
                E: E::G1Affine::zero(),
                // used by the helper function
                x: instance_2.x.iter()
                    .zip(instance_1.x.iter())
                    .map(|(&e2, &e1)| e2 * two - e1)
                    .collect(),
                y: instance_2.y.iter()
                    .zip(instance_1.y.iter())
                    .map(|(&e2, &e1)| e2 * two - e1)
                    .collect(),
                z: two * instance_2.z - instance_1.z,
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

    pub fn decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> bool {
        let instance = &acc.instance;
        let witness = &acc.witness;

        // first condition
        let pairing_lhs = E::multi_pairing(&witness.vec_D, &srs.pc_srs.vec_V);
        let pairing_rhs = E::pairing(instance.C, srs.pc_srs.V_prime);

        // second condition
        let ip_rhs = instance.T;
        let ip_lhs = {
            let res = E::G1::msm_unchecked(srs.k_x.as_slice(), witness.tree_x.nodes.as_slice());
            res.add(E::G1::msm_unchecked(srs.k_y.as_slice(), witness.tree_y.nodes.as_slice()))
        };

        // third condition
        let verify_lhs = Self::helper_function_decide(srs, acc);
        let verify_rhs = instance.E;

        println!("{}", verify_rhs == verify_lhs.into());
        println!("{}", ip_lhs == ip_rhs.into());
        println!("{}", pairing_lhs == pairing_rhs);

        return (verify_rhs == verify_lhs.into()) && (ip_lhs == ip_rhs.into()) && (pairing_lhs == pairing_rhs);
    }

    pub fn helper_function_decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> E::G1Affine {
        let instance = &acc.instance;
        let witness = &acc.witness;

        let e_prime: E::ScalarField = compute_dot_product(
            &witness.f_star_poly.evaluation_over_boolean_hypercube,
            &acc.witness.tree_y.get_leaves()
        ) - instance.z;
        let E_G = {
            let lhs = E::G1::msm_unchecked(
                srs.pc_srs.vec_H.as_slice(),
                witness.f_star_poly.evaluation_over_boolean_hypercube.as_slice(),
            );
            let rhs = E::G1::msm_unchecked(
                witness.vec_D.as_slice(),
                acc.witness.tree_x.get_leaves(),
            );
            lhs.sub(rhs)
        };

        let error_tree_x = acc.witness.tree_x.difference(acc.instance.x.as_slice());
        let error_tree_y = acc.witness.tree_y.difference(acc.instance.y.as_slice());

        let mut res: E::G1 = E_G.clone();
        res = res.add(E::G1::msm_unchecked(srs.k_x.as_slice(), error_tree_x.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_y.as_slice(), error_tree_y.nodes.as_slice()));
        res.add(srs.k_prime.mul(e_prime)).into()
    }
}

pub mod test {
    use ark_ec::pairing::Pairing;
    use ark_ff::AdditiveGroup;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::accumulation::accumulator::{AccSRS, Accumulator};
    use crate::constant_for_curves::{E, ScalarField};
    use crate::pcs::multilinear_pcs::{PolyCommit, PolyCommitTrait, SRS};
    use crate::polynomial::multilinear_polynomial::bivariate_multilinear::BivariateMultiLinearPolynomial;
    use crate::polynomial::multilinear_polynomial::math::Math;
    use crate::polynomial::multilinear_polynomial::multilinear_poly::MultilinearPolynomial;

    pub fn get_satisfying_accumulator(srs: &AccSRS<E>) -> Accumulator<E> {
        // define the polynomial commitment
        let poly_commit: PolyCommit<E> = PolyCommit { srs: srs.pc_srs.clone() };

        // random bivariate polynomial
        let polynomial1 = BivariateMultiLinearPolynomial::from_multilinear_to_bivariate_multilinear(
            MultilinearPolynomial::rand(srs.pc_srs.degree_x.log_2() + srs.pc_srs.degree_y.log_2(), &mut thread_rng()),
            srs.pc_srs.degree_x,
        );
        let polynomial2 = BivariateMultiLinearPolynomial::from_multilinear_to_bivariate_multilinear(
            MultilinearPolynomial::rand(srs.pc_srs.degree_x.log_2() + srs.pc_srs.degree_y.log_2(), &mut thread_rng()),
            srs.pc_srs.degree_x,
        );

        // random points and evaluation
        let x1 = {
            let mut res = Vec::with_capacity(srs.pc_srs.degree_x.log_2());
            for _ in 0..srs.pc_srs.degree_x.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };
        let x2 = {
            let mut res = Vec::with_capacity(srs.pc_srs.degree_x.log_2());
            for _ in 0..srs.pc_srs.degree_x.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };
        let y1 = {
            let mut res = Vec::with_capacity(srs.pc_srs.degree_y.log_2());
            for _ in 0..srs.pc_srs.degree_y.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };
        let y2 = {
            let mut res = Vec::with_capacity(srs.pc_srs.degree_y.log_2());
            for _ in 0..srs.pc_srs.degree_y.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        let whole_input_1 = {
            let mut res = vec![];
            res.extend(x1.clone());
            res.extend(y1.clone());
            res
        };
        let whole_input_2 = {
            let mut res = vec![];
            res.extend(x2.clone());
            res.extend(y2.clone());
            res
        };

        let z1 = polynomial1.poly.evaluate(&whole_input_1);
        let z2 = polynomial2.poly.evaluate(&whole_input_2);

        // commit to the polynomial
        let com1 = poly_commit.commit(&polynomial1);
        let com2 = poly_commit.commit(&polynomial2);

        // open the commitment
        let open1 = poly_commit.open(&polynomial1, com1.clone(), &x1);
        let open2 = poly_commit.open(&polynomial2, com2.clone(), &x2);

        // verify the proof
        assert!(poly_commit.verify(&com1, &open1, &x1, &y1, &z1));
        assert!(poly_commit.verify(&com2, &open2, &x2, &y2, &z2));

        let instance1 = Accumulator::new_accumulator_instance_from_proof(&srs, &com1.C, &x1, &y1, &z1);
        let witness1 = Accumulator::new_accumulator_witness_from_proof(&srs, open1, &x1, &y1);
        let instance2 = Accumulator::new_accumulator_instance_from_proof(&srs, &com2.C, &x2, &y2, &z2);
        let witness2 = Accumulator::new_accumulator_witness_from_proof(&srs, open2, &x2, &y2);

        let acc1 = Accumulator::new_accumulator(&instance1, &witness1);
        let acc2 = Accumulator::new_accumulator(&instance2, &witness2);
        assert!(Accumulator::decide(&srs, &acc1));
        assert!(Accumulator::decide(&srs, &acc2));

        let (new_instance, new_witness, Q) = Accumulator::prove(&srs, &acc1, &acc2);

        let new_acc = Accumulator::new_accumulator(&new_instance, &new_witness);
        assert!(Accumulator::decide(&srs, &new_acc));

        new_acc
    }

    #[test]
    fn test_end_to_end() {
        let degree_x = 4usize;
        let degree_y = 16usize;
        let srs_pcs: SRS<E> = PolyCommit::<E>::setup(degree_x, degree_y, &mut thread_rng());
        let srs = Accumulator::setup(srs_pcs.clone(), &mut thread_rng());

        get_satisfying_accumulator(&srs);
    }
}
use std::ops::{Add, Mul, Sub};

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_poly::EvaluationDomain;
use ark_std::UniformRand;
use rand::RngCore;

use crate::accumulation::{convert_affine_to_scalars, generate_random_elements};
use crate::accumulation::eq_tree::EqTree;
use crate::hash::poseidon::{PoseidonHash, PoseidonHashTrait};
use crate::pcs::multilinear_pcs::{OpeningProof, SRS};
use crate::polynomial::multilinear_polynomial::multilinear_poly::MultilinearPolynomial;
use crate::polynomial::traits::OneDimensionalPolynomial;
use crate::utils::{inner_product, is_power_of_two};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccSRS<E: Pairing> {
    pub degree_x: usize,
    pub degree_y: usize,

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
    pub fn setup<T: RngCore>(degree_x: usize,
                             degree_y: usize,
                             pc_srs: SRS<E>,
                             rng: &mut T,
    ) -> AccSRS<E> {
        // asserting the degree_x and degree_y are a power of two
        assert!(is_power_of_two(degree_x), "degree_x (upper bound) must be a power of two");
        assert!(is_power_of_two(degree_y), "degree_y (upper bound) must be a power of two");

        // return the result
        AccSRS {
            degree_x,
            degree_y,
            pc_srs,
            k_x: generate_random_elements::<E, T>(2 * degree_x - 1, rng),
            k_y: generate_random_elements::<E, T>(2 * degree_y - 1, rng),
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

    pub fn new_accumulator_instance_from_proof(srs: &AccSRS<E>, C: &E::G1Affine, x: &Vec<E::ScalarField>, y: &Vec<E::ScalarField>, z: &E::ScalarField) -> AccInstance<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.degree_y, "invalid size of vector y");

        let T = E::G1::msm_unchecked(srs.k_x.as_slice(), x.as_slice()).add(E::G1::msm_unchecked(srs.k_y.as_slice(), y.as_slice()));

        AccInstance {
            C: *C,
            T: T.into(),
            E: E::G1Affine::zero(),
            x: x.clone(),
            y: y.clone(),
            z: z.clone(),
        }
    }

    pub fn new_accumulator_witness_from_proof<U: OneDimensionalPolynomial<E>>(srs: &AccSRS<E>, proof: OpeningProof<E, U>, x: &Vec<E::ScalarField>, y: &Vec<E::ScalarField>) -> AccWitness<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.degree_y, "invalid size of vector y");
        assert_eq!(proof.vec_D.len(), srs.degree_x, "invalid proof size");

        AccWitness {
            vec_D: proof.vec_D,
            f_star_poly:
            tree_x: EqTree::new(x.as_slice()),
            tree_y: EqTree::new(y.as_slice()),
        }
    }

    pub fn decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> bool {
        let instance = &acc.instance;
        let witness = &acc.witness;

        true
    }

    pub fn helper_function_decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> E::G1Affine {
        let instance = &acc.instance;
        let witness = &acc.witness;

        let y_0 = acc.witness.tree_y.get_leaves().to_vec();
        let e_prime: E::ScalarField = witness.f_star_poly.evaluate(&y_0)- instance.z;

        let E_G: E::G1 = {
            let lhs = E::G1::msm_unchecked(
                srs.pc_srs.vec_H.as_slice(),
                witness.f_star_poly.evaluation_over_boolean_hypercube.as_slice()
            );
            let rhs = {
                let x_0 = acc.witness.tree_x.get_leaves().to_vec();
                E::G1::msm_unchecked(witness.vec_D.as_slice(), x_0.as_slice())
            };
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
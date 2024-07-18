/*use std::marker::PhantomData;

use crate::{pcs::{Commitment, OpeningProof, SRS}, univariate_poly::UnivariatePolynomial, utils::compute_powers};
use crate::accumulation::verifier::AccInstance;
use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ec::VariableBaseMSM;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccWitness<E: Pairing> {
    pub vec_D_i: Vec<E::G1Affine>,
    pub f_star_poly: UnivariatePolynomial<E::ScalarField>,
    pub vec_b: Vec<E::ScalarField>,
    pub vec_c: Vec<E::ScalarField>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccProver<E: Pairing> {
    _phantom: PhantomData<E>,
}

impl<E: Pairing> AccProver<E> {
    fn new_acc_witness(proof: &OpeningProof<E>, b: &E::ScalarField, c: &E::ScalarField) -> AccWitness<E> {
        let n = proof.vec_D_i.len();
        let _vec_b_powers = compute_powers(b, n);
        let _vec_c_powers = compute_powers(c, n);

        unimplemented!();
    }

    /// Accumulation prover: Given (A_1.X, A_1.W, A_2.X, A_2.W) compute (A_3.X, A_3.W)
    fn prove(_instance_1: &AccInstance<E>, _witness_1: &AccWitness<E>, _instance_2: &AccInstance<E>, _witness_2: &AccWitness<E>) -> (AccInstance<E>, AccWitness<E>) {
        unimplemented!();
    }
}



 */
/// A sqrtn polynomial commitment scheme
use std::{iter, marker::PhantomData, ops::Mul};

use ark_ec::AffineRepr;
use ark_ff::UniformRand;
use rand::RngCore;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;

use crate::{bivariate_poly::BivariatePolynomial, univariate_poly::UnivariatePolynomial, utils::compute_powers};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SRS<E: Pairing> {
    // just a helper variable (not really needed)
    pub n: usize,
    pub vec_H_i_j: Vec<E::G1Affine>,
    pub vec_H_j: Vec<E::G1Affine>,
    pub vec_V_i: Vec<E::G2Affine>,
    pub V_prime: E::G2Affine,
    // This is the vector G from the accumulation (TODO: abstraction leak!)
    pub vec_G_accumulation: Vec<E::G1Affine>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<E: Pairing> {
    C: E::G1Affine
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpeningProof<E: Pairing> {
    pub vec_D_i: Vec<E::G1Affine>,
    pub f_star_poly: UnivariatePolynomial<E::ScalarField>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoeffFormPCS<E: Pairing> {
    _phantom: PhantomData<E>,
}

impl<E: Pairing> CoeffFormPCS<E> {
    /// Create an SRS of size n^2
    pub fn generate_srs_for_testing<T: RngCore>(n: usize, rng: &mut T) -> SRS<E> {
        // Generate SRS parameters
        // G
        let generator_G1 = E::G1Affine::rand(rng);
        // V
        let generator_G2 = E::G2Affine::rand(rng);

        let vec_G_accumulation: Vec<E::G1Affine> =
            iter::repeat_with(|| E::G1Affine::rand(rng))
            .take(n)
            .collect();

        // Generate trapdoors
        let tau = E::ScalarField::rand(rng);
        let mu = E::ScalarField::rand(rng);
        let alpha = E::ScalarField::rand(rng);

        let powers_of_tau = compute_powers(&tau, n);
        let powers_of_mu = compute_powers(&mu, n);

        // Allocate space for the SRS
        let mut vec_H_i_j = Vec::with_capacity(n*n);
        let mut vec_H_j = Vec::with_capacity(n);
        let mut vec_V_i = Vec::with_capacity(n);


        // TODO: For big values of n, we need to multithread. see how bi-kzg uses rayon with parallelize()
        for i in 0..n {
            vec_H_j.push(generator_G1.mul(alpha * powers_of_mu[i]).into());
            vec_V_i.push(generator_G2.mul(powers_of_tau[i]).into());
             for j in 0..n {
                 vec_H_i_j.push(generator_G1.mul(powers_of_tau[i] * powers_of_mu[j]).into());
             }
        }
        let V_prime = generator_G2.mul(alpha).into();

        SRS {
            n,
            vec_H_i_j,
            vec_H_j,
            vec_V_i,
            V_prime,
            vec_G_accumulation,
        }
    }

    pub fn commit(poly: &BivariatePolynomial<E::ScalarField>, srs: &SRS<E>) -> Commitment<E> {
        assert_eq!(poly.degree, srs.n);
        Commitment {
            C: E::G1::msm_unchecked(&srs.vec_H_i_j, &poly.coeffs).into()
        }
    }

    /// Create a proof that f(b,c) = y
    pub fn open(poly: &BivariatePolynomial<E::ScalarField>, b: &E::ScalarField, srs: &SRS<E>) -> OpeningProof<E> {
        // Create row commitments
        let mut vec_D_i: Vec<E::G1Affine> = Vec::with_capacity(srs.n);
        for i in 0..srs.n {
            let mut D_i = E::G1Affine::zero();
            for j in 0..srs.n {
                D_i = (D_i + srs.vec_H_j[j].mul(poly.coeffs[i*srs.n + j])).into();
            }
            vec_D_i.push(D_i);
        }

        // Get the partial evaluation
        let f_star_poly = poly.partially_evaluate_at_x(&b);

        OpeningProof {
            vec_D_i,
            f_star_poly
        }
    }

    pub fn verify(C: &Commitment<E>, proof: &OpeningProof<E>, b: &E::ScalarField, c: &E::ScalarField, y: &E::ScalarField, srs: &SRS<E>) -> bool {
        // Step 1) Perform the pairing check
        let pairing_rhs = E::multi_pairing(proof.vec_D_i.clone(), &srs.vec_V_i); // XXX avoid clone()?
        // XXX this can be optimized further (do one pairing less) by taking its inverse and moving it to the multi_pairing
        let pairing_lhs = E::pairing(&C.C, &srs.V_prime);

        // Step 2) Perform the partial poly commitment check
        let powers_of_b = compute_powers(b, srs.n);
        let msm_lhs = E::G1::msm_unchecked(&srs.vec_H_j, &proof.f_star_poly.poly.coeffs);
        // XXX this can be optimized fruther by taking its inverse and doing it in one MSM
        let msm_rhs = E::G1::msm_unchecked(&proof.vec_D_i, &powers_of_b);

        // Step 3) Check the final result
        let y_expected =  proof.f_star_poly.evaluate(c);

        pairing_rhs == pairing_lhs && msm_lhs == msm_rhs && y_expected == *y
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use ark_std::test_rng;
    use ark_std::UniformRand;

    use ark_bn254::{Bn254, Fr, G1Projective};

    #[test]
    fn test_end_to_end() {
        let rng = &mut test_rng();
        let N = 4;

        let srs = CoeffFormPCS::<Bn254>::generate_srs_for_testing(N, rng);

        let bivariate_poly = BivariatePolynomial::random(rng, N);

        let commitment = CoeffFormPCS::<Bn254>::commit(&bivariate_poly, &srs);

        let b = Fr::rand(rng);
        let c = Fr::rand(rng);
        let y = bivariate_poly.evaluate(&b, &c);

        let opening_proof = CoeffFormPCS::<Bn254>::open(&bivariate_poly, &b, &srs);

        let verify_result = CoeffFormPCS::<Bn254>::verify(&commitment, &opening_proof, &b, &c, &y, &srs);
        assert!(verify_result)
    }
}

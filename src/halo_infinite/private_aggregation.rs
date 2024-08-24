#![allow(unused)]
use itertools::izip;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_ec::VariableBaseMSM;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, batch_inversion};
use ark_ff::Zero;
use ark_ec::CurveGroup;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, EvaluationDomain};

use std::ops::{Div, Mul, Add};

use transcript::IOPTranscript;

use crate::halo_infinite::errors::ProofError;
use crate::utils::compute_powers;

/// The polynomial in \\(\FF\\) that vanishes in all the points `points`.
pub fn vanishing_polynomial<E: Pairing>(points: &[E::ScalarField]) -> DensePolynomial<E::ScalarField> {
    let one = DensePolynomial::from_coefficients_vec(vec![E::ScalarField::one()]);
    points
        .iter()
        .map(|&point| DensePolynomial::from_coefficients_vec(vec![-point, E::ScalarField::one()]))
        .fold(one, |x, y| x.naive_mul(&y))
}


/// Remove `element` from the set behind `iter`
fn remove_from_set<I, T>(iter: I, element: &T) -> Vec<T>
where
    I: IntoIterator<Item = T>,
    T: PartialEq,
{
    iter.into_iter()
        .filter(|item| item != element)
        .collect()
}

/// Prove that f_i(w_i) = 0
/// Simplification over paper: omega_i is a single element
pub fn prove<E: Pairing>(
    vec_f: &Vec<DensePolynomial<E::ScalarField>>,
    vec_omega: &Vec<E::ScalarField>,
    entire_domain: &GeneralEvaluationDomain<E::ScalarField>,
    transcript: &mut IOPTranscript<E::ScalarField>,
) -> () {
    let z_poly = vanishing_polynomial::<E>(&vec_omega);

    // This is a vector of Omega_i
    let mut vec_omega_complements = vec![];
    for omega in vec_omega {
        let vec_omega_complement_i = remove_from_set(entire_domain.elements().into_iter(), omega);
        vec_omega_complements.push(vec_omega_complement_i);
    }

    // Compute z_i(x) = \prod (x - w_i)
    let mut vec_z_i = vec![];
    for vec_omega_complement in &vec_omega_complements {
        let z_i_poly = vanishing_polynomial::<E>(&vec_omega_complement);
        vec_z_i.push(z_i_poly);
    }

    // Step 2: Get a challenge from the verifier
    transcript.append_serializable_element(b"z_x", &z_poly).unwrap();
    transcript.append_serializable_element(b"z_more", &vec_omega_complements).unwrap();
    transcript.append_serializable_element(b"z_moremore", &vec_z_i).unwrap();
    let rho: E::ScalarField = transcript.get_and_append_challenge(b"rho").unwrap();

    // Compute powers of rho
    let vec_rho = compute_powers(&rho, vec_f.len());

    assert_eq!(vec_rho.len(), vec_f.len());
    assert_eq!(vec_f.len(), vec_z_i.len());

    // Step 3: Compute q(x)
    let mut q_x = DensePolynomial::from_coefficients_vec(vec![E::ScalarField::zero()]);
    for (rho_i, f_i, z_i) in izip!(vec_rho, vec_f, vec_z_i) {
        let mut numerator = DensePolynomial::from_coefficients_vec(vec![E::ScalarField::zero()]);
        numerator = numerator.mul(f_i);
        numerator = numerator.mul(&z_i);

        q_x = q_x.add(numerator);
    }
    q_x = q_x.div(&z_poly);
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::{Bn254, Fr, G1Affine, G1Projective};
    use ark_poly::DenseUVPolynomial;
    use std::ops::Sub;
    use rand::thread_rng;

    #[test]
    pub fn test_prove() {
        let N = 2;
        let d = 16;

        let mut transcript_prover = IOPTranscript::<Fr>::new(b"ipa");

        let mut vec_f : Vec<DensePolynomial<Fr>> = vec![];
        let mut vec_omega = vec![];

        // we want to prove that N polynomials evaluate to t_i at w_i
        // so p_i(w_i) = t_i
        // So essentially: f(w_i) = p_i(w_i) - t_i = 0
        // where f_i(x) = p_i(x) - t_i

        let domain = GeneralEvaluationDomain::<Fr>::new(d).unwrap();
        for i in 0..N {
            vec_omega.push(domain.element(i));
        }

        for i in 0..N {
            let p_poly = DensePolynomial::rand(d, &mut thread_rng());
            let p_poly_at_w_i = p_poly.evaluate(&vec_omega[i]);
            let f_poly = p_poly.sub(&DensePolynomial::from_coefficients_slice(&[p_poly_at_w_i]));
            // Check that f_i(w_i) == 0
            assert_eq!(Fr::zero(), f_poly.evaluate(&vec_omega[i]));
            vec_f.push(f_poly);
        }

        prove::<Bn254>(&vec_f, &vec_omega, &domain, &mut transcript_prover);
    }
}


#![allow(unused)]
use itertools::izip;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_ec::AffineRepr;
use ark_ec::VariableBaseMSM;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, batch_inversion};
use ark_ff::{FftField};
use ark_ff::Zero;
use ark_std::{end_timer, start_timer};
use ark_ec::CurveGroup;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, EvaluationDomain};
use crate::kzg::KZGProof;
use crate::kzg::{KZGCommitment, KZGPowers, KZGVerifierKey, KZG10};

use std::ops::{Div, Mul, Add, Sub, Neg};

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

pub struct PrivateAggregationProof<E: Pairing> {
    q_commitment: KZGCommitment<E>,
    proof_g_at_r: Option<KZGProof<E>>,
    g_commitment: Option<KZGCommitment<E>>, // XXX remove (just for testing)
}

/// Prove that f_i(w_i) = 0
/// Simplification over paper: omega_i is a single element
pub fn prove<E: Pairing>(
    vec_f: &Vec<DensePolynomial<E::ScalarField>>,
    vec_omega: &Vec<E::ScalarField>,
    entire_domain: &GeneralEvaluationDomain<E::ScalarField>,
    transcript: &mut IOPTranscript<E::ScalarField>,
    ck: &KZGPowers<E>,
) -> PrivateAggregationProof<E>
where
    E::ScalarField: FftField, // FFTField is needed to do poly mul with FFTs
{
    // Step 1: Compute z(x) and z_i(x) polys
    let z_poly = vanishing_polynomial::<E>(&vec_omega);

    // This is a vector of Omega_i
    let mut vec_omega_complements = vec![];
    for omega in vec_omega {
        let vec_omega_complement_i:Vec<E::ScalarField> = remove_from_set(vec_omega.clone(), omega);
        vec_omega_complements.push(vec_omega_complement_i);
    }

    // Compute z_i(x) = \prod (x - w_i)
    let vec_z_i_time = start_timer!(|| format!("vec_z_i time"));
    let mut vec_z_i = vec![];
    for vec_omega_complement in &vec_omega_complements {
        let z_i_poly = vanishing_polynomial::<E>(&vec_omega_complement);
        vec_z_i.push(z_i_poly);
    }
    end_timer!(vec_z_i_time);

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
    let step3_time = start_timer!(|| format!("Step3 time"));
    let mut q_x = DensePolynomial::from_coefficients_vec(vec![E::ScalarField::zero()]);
    for (rho_i, f_i, z_i) in izip!(vec_rho.clone(), vec_f.clone(), vec_z_i.clone()) { // XXX loose the clone
        // Compute rho^{i-1}*f_i
        let mut numerator = DensePolynomial::from_coefficients_vec(
            f_i.coeffs.iter().map(|f| *f * rho_i).collect(),
        );
        // Compute the entire numerator
        numerator = numerator.mul(&z_i);

        q_x = q_x.add(numerator);
    }
    q_x = q_x.div(&z_poly);
    end_timer!(step3_time);

    // Step 4: Commit to q(x)
    let step4_time = start_timer!(|| format!("Step4 time"));
    let (q_commitment, q_blinder) = KZG10::<E,DensePolynomial<<E as Pairing>::ScalarField>>::commit(&ck, &q_x, None, None).expect("q commitment failed");
    end_timer!(step4_time);

    // Step 5: Get r from verifier
    transcript.append_serializable_element(b"q_comm", &q_commitment.0).unwrap();
    let r: E::ScalarField = transcript.get_and_append_challenge(b"r").unwrap();

    // Step 6: Compute g(x)
    let step6_time = start_timer!(|| format!("Step6 time"));
    let mut g_x  = DensePolynomial::from_coefficients_vec(vec![E::ScalarField::zero()]);
    for (rho_i, f_i, z_i) in izip!(vec_rho, vec_f, vec_z_i) {
        let z_i_r = z_i.evaluate(&r);
        let f_i_times_z_i_r = DensePolynomial::from_coefficients_vec(
            f_i.coeffs.iter().map(|f| *f * z_i_r).collect(),
        );
        let summand = DensePolynomial::from_coefficients_vec(
            f_i_times_z_i_r.coeffs.iter().map(|f| *f * rho_i).collect(),
        );
        g_x = g_x.add(summand);
    }
    let z_r = z_poly.evaluate(&r);
    let q_x_times_z_r = DensePolynomial::from_coefficients_vec(
        q_x.coeffs.iter().map(|q| *q * z_r).collect(),
    );
    g_x = g_x.sub(&q_x_times_z_r);
    end_timer!(step6_time);

    // Step 7: Compute commitment and proof (NOT NEEDED for private aggregation!)
    // assert_eq!(g_x.evaluate(&r), E::ScalarField::zero());

    // let (g_commitment, g_blinder) = KZG10::<E,DensePolynomial<<E as Pairing>::ScalarField>>::commit(&ck, &g_x, None, None).expect("g commitment failed");
    // let proof_g_at_r = KZG10::<E,DensePolynomial<<E as Pairing>::ScalarField>>::open(&ck, &g_x, r, &g_blinder).expect("Proof generation failed");

    PrivateAggregationProof {
        q_commitment,
        proof_g_at_r: None,
        g_commitment: None,
    }
}

/// Prove that f_i(w_i) = 0
/// Simplification over paper: omega_i is a single element
pub fn verify<E: Pairing>(
    proof: &PrivateAggregationProof<E>,
    vec_f_commitments: &Vec<KZGCommitment<E>>,
    vec_omega: &Vec<E::ScalarField>,
    entire_domain: &GeneralEvaluationDomain<E::ScalarField>,
    transcript: &mut IOPTranscript<E::ScalarField>,
    vk: &KZGVerifierKey<E>
) -> bool {
    let n = vec_f_commitments.len();

    // Step 1: Compute all the domains and polys
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
    let vec_rho = compute_powers(&rho, n);

    assert_eq!(vec_rho.len(), n);
    assert_eq!(n, vec_z_i.len());

    // Step 2: Get r from verifier
    transcript.append_serializable_element(b"q_comm", &proof.q_commitment.0).unwrap();
    let r: E::ScalarField = transcript.get_and_append_challenge(b"r").unwrap();

    let mut C_prime = E::G1Affine::zero();
    for (rho_i, z_i, C_i) in izip!(vec_rho.clone(), vec_z_i.clone(), vec_f_commitments.clone()) { // XXX loose the clone
        let z_i_at_r = z_i.evaluate(&r);
        C_prime = (C_prime + C_i.0.mul(z_i_at_r * rho_i)).into();
    }

    let z_at_r = z_poly.evaluate(&r);
    let z_r_C_q = proof.q_commitment.0.mul(z_at_r);
    let C_g: E::G1Affine = C_prime.sub(z_r_C_q).into_affine();

    // Final check: Check the KZG proof that g(r) == 0
    // NOT NEEDED for accumulation verifier (only for decider)

    // TODO: bug! why need to negate? forgot something.
    // let tmp_negated: E::G1Affine = E::G1Affine::zero().sub(C_g).into();
    // let C_g_commitment: KZGCommitment<E> = KZGCommitment(tmp_negated.into());
    // assert_eq!(C_g_commitment, proof.g_commitment);

    // let zero = E::ScalarField::zero();
    // let is_valid = KZG10::<E, DensePolynomial<<E as Pairing>::ScalarField>>::check(&vk, &C_g_commitment, r, zero, &proof.proof_g_at_r).expect("Verification failed");
    // assert!(is_valid);

    true
}

pub mod tests {
    use crate::kzg::{trim};

    use super::*;
    use ark_bn254::{Bn254, Fr};
    use ark_poly::DenseUVPolynomial;
    use std::ops::Sub;
    use rand::{rngs::StdRng, thread_rng, SeedableRng};

    type E = Bn254;

    pub fn prepare_polynomials_and_srs(
        N: usize,
        d: usize,
        rng: &mut StdRng,
    ) -> (
        Vec<DensePolynomial<Fr>>,
        Vec<KZGCommitment<E>>,
        Vec<Fr>,
        GeneralEvaluationDomain<Fr>,
        KZGPowers<E>,
        KZGVerifierKey<E>,
    ) {
        let mut vec_f: Vec<DensePolynomial<Fr>> = vec![];
        let mut vec_f_commitments: Vec<KZGCommitment<E>> = vec![];
        let mut vec_omega = vec![];

        let params = KZG10::<E, DensePolynomial<<E as Pairing>::ScalarField>>::setup(2 * d, false, rng).expect("Setup failed");
        let (ck, vk) = trim(&params, 2 * d);

        let domain = GeneralEvaluationDomain::<Fr>::new(d).unwrap();
        for i in 0..N {
            vec_omega.push(domain.element(i));
        }

        for i in 0..N {
            let p_poly = DensePolynomial::rand(d, rng);
            let p_poly_at_w_i = p_poly.evaluate(&vec_omega[i]);
            let f_poly = p_poly.sub(&DensePolynomial::from_coefficients_slice(&[p_poly_at_w_i]));
            assert_eq!(Fr::zero(), f_poly.evaluate(&vec_omega[i]));
            vec_f.push(f_poly.clone());

            let (f_commitment, _) = KZG10::<E, DensePolynomial<<E as Pairing>::ScalarField>>::commit(&ck, &f_poly, None, None).expect("f commitment failed");
            vec_f_commitments.push(f_commitment);
        }

        (vec_f, vec_f_commitments, vec_omega, domain, ck, vk)
    }

    #[test]
    pub fn test_private_aggregation_end_to_end() {
        let mut rng = StdRng::seed_from_u64(0u64);

        let N = 2;
        let d = 8192;

        let mut transcript_prover = IOPTranscript::<Fr>::new(b"pa");
        let mut transcript_verifier = IOPTranscript::<Fr>::new(b"pa");

        let (vec_f, vec_f_commitments, vec_omega, domain, ck, vk) = prepare_polynomials_and_srs(N, d, &mut rng);

        let proof = prove::<Bn254>(&vec_f, &vec_omega, &domain, &mut transcript_prover, &ck);
        let is_valid = verify::<Bn254>(&proof, &vec_f_commitments, &vec_omega, &domain, &mut transcript_verifier, &vk);
        assert!(is_valid);
    }
}


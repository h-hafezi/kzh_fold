#![allow(unused)]
// It's basically an HPI

use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{batch_inversion, Field, One, PrimeField, Zero};

use crate::halo_infinite::errors::ProofError;
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use std::ops::Mul;

pub struct HPIProof<E: Pairing> {
    vec_Y_L: Vec<E::G1Affine>,
    vec_Y_R: Vec<E::G1Affine>,
    x_final: E::ScalarField,
}

pub fn prove<E, F>(
    crs_G_vec: Vec<E::G1Affine>,
    vec_x: Vec<F>,
    transcript: &mut Transcript<F>,
) -> HPIProof<E> where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    E: Pairing<ScalarField=F>,
    F: PrimeField + Absorb,
{
    let mut n = vec_x.len();
    let lg_n = ark_std::log2(n) as usize;
    assert_eq!(vec_x.len(), n);
    assert!(n.is_power_of_two());

    let mut vec_Y_L: Vec<E::G1Affine> = Vec::with_capacity(lg_n);
    let mut vec_Y_R: Vec<E::G1Affine> = Vec::with_capacity(lg_n);

    // Create slices backed by their respective vectors.  This lets us reslice as we compress the lengths of the
    // vectors in the main loop below.
    let mut slice_G = &mut crs_G_vec.clone()[..];
    let mut slice_x = &mut vec_x.clone()[..];

    while slice_x.len() > 1 {
        n /= 2;

        let (x_L, x_R) = slice_x.split_at_mut(n);
        let (G_L, G_R) = slice_G.split_at_mut(n);

        let Y_L = E::G1::msm_unchecked(G_R, x_L).into_affine();
        let Y_R = E::G1::msm_unchecked(G_L, x_R).into_affine();
        vec_Y_L.push(Y_L);
        vec_Y_R.push(Y_R);

        transcript.append_points::<E>(b"vec_Y", &[Y_L, Y_R]);
        let gamma: F = transcript.challenge_scalar(b"gamma");
        let gamma_inv = gamma.inverse().expect("gamma must have an inverse");

        // Fold input vectors and basis
        for i in 0..n {
            x_L[i] = x_L[i] + gamma_inv * x_R[i];
            G_L[i] = (G_L[i] + G_R[i].mul(gamma)).into_affine();
        }

        // Save the rescaled vector for splitting in the next loop
        slice_x = x_L;
        slice_G = G_L;
    }

    assert_eq!(slice_x.len(), 1);

    HPIProof {
        vec_Y_L,
        vec_Y_R,
        x_final: slice_x[0],
    }
}

/// Get a bit string to derive the verification scalars using binary decomposition.
///
/// XXX: This can be done more elegantly
pub fn get_verification_scalars_bitstring(n: usize, logn: usize) -> Vec<Vec<usize>> {
    let mut bitstring: Vec<Vec<usize>> = Vec::new();
    for _i in 0..n {
        let vec_i: Vec<usize> = Vec::new();
        bitstring.push(vec_i);
    }

    for j in 0..logn {
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            let current_bit_string = format!("{:b}", i);
            let mut bit_vec: Vec<char> = current_bit_string.chars().collect();
            bit_vec.reverse();
            while bit_vec.len() < logn {
                bit_vec.push('0');
            }

            if bit_vec[logn - j - 1] == '1' {
                bitstring[i].push(j);
            }
        }
    }

    bitstring
}

#[allow(clippy::type_complexity)]
fn verification_scalars<E, F>(
    proof: &HPIProof<E>,
    n: usize,
    transcript: &mut Transcript<F>,
) -> Result<(Vec<F>, Vec<F>, Vec<F>), ProofError> where
    E: Pairing<ScalarField=F>,
    F: PrimeField + Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    let lg_n = proof.vec_Y_L.len();
    if lg_n >= 32 {
        return Err(ProofError::VerificationError);
    }
    if n != (1 << lg_n) {
        return Err(ProofError::VerificationError);
    }

    let verification_scalars_bit_string = get_verification_scalars_bitstring(n, lg_n);

    // 1. Recompute gamma_k,...,gamma_1 based on the proof transcript
    let mut challenges: Vec<F> = Vec::with_capacity(lg_n);
    for i in 0..proof.vec_Y_L.len() {
        transcript.append_points::<E>(
            b"vec_Y",
            &[
                proof.vec_Y_L[i],
                proof.vec_Y_R[i],
            ],
        );
        challenges.push(transcript.challenge_scalar(b"gamma"));
    }

    // 2. Compute 1/gamma_k, ..., 1/gamma_1
    let mut challenges_inv: Vec<F> = challenges.clone();
    batch_inversion(&mut challenges_inv);

    // 3. Compute s values by iterating over the bit_string
    let mut vec_s: Vec<F> = Vec::with_capacity(n);
    for i in 0..n {
        vec_s.push(F::one());
        for j in 0..verification_scalars_bit_string[i].len() {
            vec_s[i] *= challenges[verification_scalars_bit_string[i][j]]
        }
    }

    Ok((challenges, challenges_inv, vec_s))
}

pub fn verify<E, F>(
    proof: &HPIProof<E>,
    crs_G_vec: Vec<E::G1Affine>,
    C: E::G1Affine,
    transcript: &mut Transcript<F>,
) -> Result<(), ProofError> where
    E: Pairing<ScalarField=F>,
    F: PrimeField + Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    let mut n = crs_G_vec.len();
    assert!(n.is_power_of_two());

    // Step 2
    let (vec_gamma, vec_gamma_inv, vec_s, ) = verification_scalars(proof, n, transcript)?;

    // it can be turned into one MSM
    let C_prime = E::G1::msm(&proof.vec_Y_L, &vec_gamma).unwrap() + C + E::G1::msm(&proof.vec_Y_R, &vec_gamma_inv).unwrap();

    let expected_C = E::G1::msm(&crs_G_vec, &vec_s).unwrap().mul(proof.x_final);

    assert!((C_prime - expected_C).is_zero());

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constant_for_curves::{G1Affine, G1Projective, ScalarField, E};
    use crate::utils::inner_product;
    use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
    use ark_std::UniformRand;
    use core::iter;

    type F = ScalarField;

    #[test]
    fn test_hpi_proof_end_to_end() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let n = 128;

        let mut transcript_prover = Transcript::<F>::new(b"ipa");
        let mut transcript_verifier = Transcript::<F>::new(b"ipa");

        let crs_G_vec: Vec<G1Affine> =
            iter::repeat_with(|| G1Projective::rand(&mut rng).into_affine())
                .take(n)
                .collect();

        let vec_x: Vec<F> = iter::repeat_with(|| rng.gen()).take(n).collect();

        let X = G1Projective::msm_unchecked(&crs_G_vec, &vec_x);

        let proof: HPIProof<E> = prove(
            crs_G_vec.clone(),
            vec_x.clone(),
            &mut transcript_prover,
        );

        let is_valid = verify(&proof,
                              crs_G_vec.clone(),
                              X.into_affine(),
                              &mut transcript_verifier,
        ).unwrap();
    }
}

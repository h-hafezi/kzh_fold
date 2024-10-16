/*#![allow(unused)]
// It's basically an HPI

use ark_ec::VariableBaseMSM;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, batch_inversion, Zero};
use ark_ec::CurveGroup;

use std::ops::Mul;

use crate::halo_infinite::errors::ProofError;
use crate::transcript::transcript::Transcript;

/// Return the inner product of two field vectors
pub fn inner_product<Fr: Field>(a: &[Fr], b: &[Fr]) -> Fr {
    assert_eq!(a.len(), b.len());
    let mut c: Fr = Fr::zero();
    for i in 0..a.len() {
        c += a[i] * b[i];
    }
    c
}

pub struct HPIProof<E: Pairing> {
//    B_c: E::G1Projective,
//    B_d: E::G1Projective,

    vec_Y_L: Vec<E::G1Affine>,
    vec_Y_R: Vec<E::G1Affine>,

    x_final: E::ScalarField,
}

pub fn prove<E: Pairing>(
    crs_G_vec: Vec<E::G1Affine>,

    vec_x: Vec<E::ScalarField>,
    // _y: E::ScalarField, // XXX

    transcript: &mut Transcript<E::ScalarField>,
) -> HPIProof<E> {
    let mut n = vec_x.len();
    let lg_n = ark_std::log2(n) as usize;
    assert_eq!(vec_x.len(), n);
    assert!(n.is_power_of_two());

    let mut vec_Y_L: Vec<E::G1Affine> = Vec::with_capacity(lg_n);
    let mut vec_Y_R: Vec<E::G1Affine> = Vec::with_capacity(lg_n);

    // Create slices backed by their respective vectors.  This lets us reslice as we compress the lengths of the
    // vectors in the main loop below.
    let mut slice_G = &mut crs_G_vec.clone()[..]; // XXX clone
    let mut slice_x = &mut vec_x.clone()[..];

    while slice_x.len() > 1 {
        n /= 2;

        let (x_L, x_R) = slice_x.split_at_mut(n);
        let (G_L, G_R) = slice_G.split_at_mut(n);

        let Y_L = E::G1::msm_unchecked(G_R, x_L).into_affine();
        let Y_R = E::G1::msm_unchecked(G_L, x_R).into_affine();
        vec_Y_L.push(Y_L);
        vec_Y_R.push(Y_R);

        transcript.append_serializable_element(b"vec_Y", &[Y_L, Y_R]).unwrap(); // XXX: Need to add the statement to the FS
        let gamma: E::ScalarField = transcript.get_and_append_challenge(b"gamma").unwrap();
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

/// Get a bitstring to derive the verification scalars using binary decomposition.
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
            let current_bitstring = format!("{:b}", i);
            let mut bit_vec: Vec<char> = current_bitstring.chars().collect();
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


// DOCDOC
#[allow(clippy::type_complexity)]
fn verification_scalars<E: Pairing>(
    proof: &HPIProof<E>,
    n: usize,
    transcript: &mut Transcript<E::ScalarField>,
) -> Result<(Vec<E::ScalarField>, Vec<E::ScalarField>, Vec<E::ScalarField>), ProofError> {
    let lg_n = proof.vec_Y_L.len();
    if lg_n >= 32 {
        return Err(ProofError::VerificationError);
    }
    if n != (1 << lg_n) {
        return Err(ProofError::VerificationError);
    }

    let verification_scalars_bitstring = get_verification_scalars_bitstring(n, lg_n);

    // 1. Recompute gamma_k,...,gamma_1 based on the proof transcript
    let mut challenges: Vec<E::ScalarField> = Vec::with_capacity(lg_n);
    for i in 0..proof.vec_Y_L.len() {
        transcript.append_serializable_element(
            b"vec_Y",
            &[
                proof.vec_Y_L[i],
                proof.vec_Y_R[i],
            ]
        ).unwrap();
        challenges.push(transcript.get_and_append_challenge(b"gamma").unwrap());
    }

    // 2. Compute 1/gamma_k, ..., 1/gamma_1
    let mut challenges_inv: Vec<E::ScalarField> = challenges.clone();
    batch_inversion(&mut challenges_inv);

    // 3. Compute s values by iterating over the bitstring
    let mut vec_s: Vec<E::ScalarField> = Vec::with_capacity(n);
    for i in 0..n {
        vec_s.push(E::ScalarField::one());
        for j in 0..verification_scalars_bitstring[i].len() {
            vec_s[i] *= challenges[verification_scalars_bitstring[i][j]]
        }
    }

    Ok((challenges, challenges_inv, vec_s))
}

pub fn verify<E: Pairing>(
    proof: &HPIProof<E>,
    crs_G_vec: Vec<E::G1Affine>,
    C: E::G1Affine,

    transcript: &mut Transcript<E::ScalarField>
) -> Result<(), ProofError> {
    let mut n = crs_G_vec.len();
    assert!(n.is_power_of_two());

    // Step 2
    let (vec_gamma, vec_gamma_inv, vec_s,) = verification_scalars(proof, n, transcript)?;

    // TODO can be turned into one MSM
    let C_prime = E::G1::msm(&proof.vec_Y_L, &vec_gamma).unwrap() + C + E::G1::msm(&proof.vec_Y_R, &vec_gamma_inv).unwrap();

    let expected_C = E::G1::msm(&crs_G_vec, &vec_s).unwrap().mul(proof.x_final);

    assert!((C_prime - expected_C).is_zero());

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
    use ark_std::UniformRand;
    use core::iter;

    use ark_bn254::{Bn254, Fr, G1Affine, G1Projective};

    #[test]
    fn test_inner_product() {
        let a = vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
        ];
        let b = vec![
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(5u64),
        ];
        assert_eq!(Fr::from(40u64), inner_product(&a, &b));
    }

    #[test]
    fn test_hpi_proof_end_to_end() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let n = 128;

        let mut transcript_prover = Transcript::<Fr>::new(b"ipa");
        let mut transcript_verifier = Transcript::<Fr>::new(b"ipa");

        let crs_G_vec: Vec<G1Affine> =
            iter::repeat_with(|| G1Projective::rand(&mut rng).into_affine())
                .take(n)
                .collect();

        let vec_x: Vec<Fr> = iter::repeat_with(|| rng.gen()).take(n).collect();

        let X = G1Projective::msm_unchecked(&crs_G_vec, &vec_x);

        let proof: HPIProof<Bn254> = prove(
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

 */

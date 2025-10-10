#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]

use ark_crypto_primitives::sponge::Absorb;
use super::errors::ProofVerifyError;
use crate::nexus_spartan::unipoly::unipoly::{CompressedUniPoly, UniPoly};
use ark_ec::CurveGroup;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::*;

use itertools::izip;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};

#[derive(CanonicalSerialize, CanonicalDeserialize, Debug, Clone)]
pub struct SumcheckInstanceProof<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPoly<F>>,
}

impl<F: PrimeField + Absorb> SumcheckInstanceProof<F> {
    pub fn new(compressed_polys: Vec<CompressedUniPoly<F>>) -> SumcheckInstanceProof<F> {
        SumcheckInstanceProof { compressed_polys }
    }

    pub fn verify<E>(
        &self,
        claim: F,
        num_rounds: usize,
        degree_bound: usize,
        transcript: &mut Transcript<F>,
    ) -> Result<(F, Vec<F>), ProofVerifyError>
    where
        E: Pairing<ScalarField = F>,
    {
        let mut e = claim;
        let mut r: Vec<F> = Vec::new();

        // verify that there is a univariate polynomial for each round
        assert_eq!(self.compressed_polys.len(), num_rounds);
        for i in 0..self.compressed_polys.len() {
            let poly = self.compressed_polys[i].decompress(&e);

            // verify degree bound
            assert_eq!(poly.degree(), degree_bound);

            // check if G_k(0) + G_k(1) = e
            assert_eq!(poly.eval_at_zero() + poly.eval_at_one(), e);

            // append the prover's message to the transcript
            <UniPoly<F> as AppendToTranscript<F>>::append_to_transcript(&poly, b"poly", transcript);

            //derive the verifier's challenge for the next round
            let r_i = Transcript::challenge_scalar(transcript, b"challenge_nextround");

            r.push(r_i);

            // evaluate the claimed degree-ell polynomial at r_i
            e = poly.evaluate(&r_i);
        }

        Ok((e, r))
    }
}


// This is used for the signature aggregation protocol!!!
impl<F: PrimeField + Absorb> SumcheckInstanceProof<F> {
    pub fn prove_cubic_four_terms<Func, G>(
        claim: &F,
        num_rounds: usize,
        poly_A: &mut MultilinearPolynomial<F>,
        poly_B: &mut MultilinearPolynomial<F>,
        poly_C: &mut MultilinearPolynomial<F>,
        poly_D: &mut MultilinearPolynomial<F>,
        comb_func: Func,
        transcript: &mut Transcript<F>,
    ) -> (Self, Vec<F>, Vec<F>)
    where
        Func: Fn(&F, &F, &F, &F) -> F,
        G: CurveGroup<ScalarField = F>,
    {
        let mut e = *claim;
        let mut r: Vec<F> = Vec::new();
        let mut cubic_polys: Vec<CompressedUniPoly<F>> = Vec::new();
        for _j in 0..num_rounds {
            let mut eval_point_0 = F::zero();
            let mut eval_point_2 = F::zero();
            let mut eval_point_3 = F::zero();

            let len = poly_A.len() / 2;
            for i in 0..len {
                // eval 0: bound_func is A(low)
                eval_point_0 += comb_func(&poly_A[i], &poly_B[i], &poly_C[i], &poly_D[i]);

                // eval 2: bound_func is -A(low) + 2*A(high)
                let poly_A_bound_point = poly_A[len + i] + poly_A[len + i] - poly_A[i];
                let poly_B_bound_point = poly_B[len + i] + poly_B[len + i] - poly_B[i];
                let poly_C_bound_point = poly_C[len + i] + poly_C[len + i] - poly_C[i];
                let poly_D_bound_point = poly_D[len + i] + poly_D[len + i] - poly_D[i];

                eval_point_2 += comb_func(
                    &poly_A_bound_point,
                    &poly_B_bound_point,
                    &poly_C_bound_point,
                    &poly_D_bound_point,
                );

                // eval 3: bound_func is -2A(low) + 3A(high); computed incrementally with bound_func applied to eval(2)
                let poly_A_bound_point = poly_A_bound_point + poly_A[len + i] - poly_A[i];
                let poly_B_bound_point = poly_B_bound_point + poly_B[len + i] - poly_B[i];
                let poly_C_bound_point = poly_C_bound_point + poly_C[len + i] - poly_C[i];
                let poly_D_bound_point = poly_D_bound_point + poly_D[len + i] - poly_D[i];

                eval_point_3 += comb_func(
                    &poly_A_bound_point,
                    &poly_B_bound_point,
                    &poly_C_bound_point,
                    &poly_D_bound_point,
                );
            }

            let evaluations = vec![eval_point_0, e - eval_point_0, eval_point_2, eval_point_3];
            let poly = UniPoly::from_evals(&evaluations);

            // append the prover's message to the transcript
            transcript.append_scalars(b"poly", poly.coeffs.as_slice());

            //derive the verifier's challenge for the next round
            let r_j: F = transcript.challenge_scalar(b"challenge_nextround");

            r.push(r_j);
            // bound all tables to the verifier's challenge
            poly_A.bound_poly_var_top(&r_j);
            poly_B.bound_poly_var_top(&r_j);
            poly_C.bound_poly_var_top(&r_j);
            poly_D.bound_poly_var_top(&r_j);
            e = poly.evaluate(&r_j);
            cubic_polys.push(poly.compress());
        }
        (
            SumcheckInstanceProof::new(cubic_polys),
            r,
            vec![poly_A[0], poly_B[0], poly_C[0], poly_D[0]],
        )
    }
}

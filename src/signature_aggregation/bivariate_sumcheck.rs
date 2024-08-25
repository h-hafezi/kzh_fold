use ark_ec::pairing::Pairing;
use transcript::IOPTranscript;
use crate::polynomial::bivariate_polynomial::bivariate_poly::BivariatePolynomial;
use crate::polynomial::bivariate_polynomial::univariate_poly::UnivariatePolynomial;

pub struct SumcheckProof<E: Pairing> {
    r_poly: UnivariatePolynomial<E::ScalarField>, // first round polynomial
    s_poly: UnivariatePolynomial<E::ScalarField>, // second round polynomial
    // more more more
}

pub fn prove<E: Pairing>(f_poly: &BivariatePolynomial<E::ScalarField>, transcript: &mut IOPTranscript<E::ScalarField>) -> (SumcheckProof<E>, (E::ScalarField, E::ScalarField))  {
    let r_poly = f_poly.sum_partial_evaluations_in_domain();

    // Squeeze alpha challenge
    transcript.append_serializable_element(b"r_poly", &r_poly).unwrap();
    let alpha = transcript.get_and_append_challenge(b"alpha").unwrap();

    let s_poly = f_poly.partially_evaluate_at_x(&alpha);

    transcript.append_serializable_element(b"s_poly", &s_poly).unwrap();
    let beta = transcript.get_and_append_challenge(b"beta").unwrap();

    let proof = SumcheckProof {
        r_poly,
        s_poly
    };

    (proof, (alpha, beta))
}

pub fn verify<E: Pairing>(proof: SumcheckProof<E>, sumcheck_result: E::ScalarField, transcript: &mut IOPTranscript<E::ScalarField>) -> bool {
    // Squeeze alpha challenge
    transcript.append_serializable_element(b"r_poly", &proof.r_poly).unwrap();
    let _alpha = transcript.get_and_append_challenge(b"alpha").unwrap();

    // Squeeze beta challenge
    transcript.append_serializable_element(b"s_poly", &proof.s_poly).unwrap();
    let _beta = transcript.get_and_append_challenge(b"beta").unwrap();

    assert_eq!(sumcheck_result, proof.r_poly.sum_evaluations_in_domain());

    // TODO need to handle evaluation of f(alpha, beta).
    // we will need a different sumcheck API. reduce/verify like in gemini

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
    use ark_std::UniformRand;
    use core::iter;
    use ark_ff::Zero;

    use crate::constant_for_curves::{ScalarField, E};

    pub fn sum_bivariate_poly_over_domain(poly: &BivariatePolynomial<ScalarField>) -> ScalarField {
        poly.evaluations.iter().cloned().sum()
    }

    #[test]
    fn test_bivariate_sumcheck_end_to_end() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let degree_x = 16usize;
        let degree_y = 4usize;
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();

        let mut transcript_prover = IOPTranscript::<ScalarField>::new(b"bsumcheck");
        let mut transcript_verifier = IOPTranscript::<ScalarField>::new(b"bsumcheck");

        let f_poly = BivariatePolynomial::random(&mut rng, domain_x, domain_y, degree_x, degree_y);
        let sumcheck_result = sum_bivariate_poly_over_domain(&f_poly);

        let (proof, (_,_)) = prove::<E>(&f_poly, &mut transcript_prover);

        verify(proof, sumcheck_result, &mut transcript_verifier);
    }
}

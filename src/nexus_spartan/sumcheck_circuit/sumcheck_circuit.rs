use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::nexus_spartan::unipoly::unipoly::CompressedUniPoly;
use crate::nexus_spartan::unipoly::unipoly_var::{CompressedUniPolyVar, UniPolyVar};
use crate::transcript::transcript::{Transcript};
use crate::transcript::transcript_var::{AppendToTranscriptVar, TranscriptVar};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;

pub struct SumcheckCircuit<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPoly<F>>,
    pub claim: F,
    pub num_rounds: usize,
    pub degree_bound: usize,
    pub transcript: Transcript<F>,
}

impl<F: PrimeField + Absorb> SumcheckCircuit<F> {
    pub fn verify<E: Pairing<ScalarField=F>>(&mut self) -> (F, Vec<F>) {
        let proof = SumcheckInstanceProof::new(self.compressed_polys.clone());
        proof.verify::<E>(self.claim, self.num_rounds, self.degree_bound, &mut self.transcript).unwrap()
    }
}

pub struct SumCheckCircuitVar<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPolyVar<F>>,
    claim: FpVar<F>,
    num_rounds: usize,
    degree_bound: usize,
    transcript: TranscriptVar<F>,
}

impl<F: PrimeField + Absorb> SumCheckCircuitVar<F> {
    pub fn verify(&mut self) -> (FpVar<F>, Vec<FpVar<F>>) {
        let mut e = self.claim.clone();
        let mut r: Vec<FpVar<F>> = Vec::new();

        // verify that there is a univariate polynomial for each round
        assert_eq!(self.compressed_polys.len(), self.num_rounds);
        for i in 0..self.compressed_polys.len() {
            let poly = self.compressed_polys[i].decompress(&e);

            // verify degree bound
            assert_eq!(poly.degree(), self.degree_bound);

            // check if G_k(0) + G_k(1) = e
            (poly.eval_at_zero() + poly.eval_at_one()).enforce_equal(&e).expect("equality error");

            // append the prover's message to the transcript
            UniPolyVar::append_to_transcript(&poly, b"poly", &mut self.transcript);

            //derive the verifier's challenge for the next round
            let r_i = TranscriptVar::challenge_scalar(&mut self.transcript, b"challenge_nextround");

            r.push(r_i.clone());

            // evaluate the claimed degree-ell polynomial at r_i
            e = poly.evaluate(&r_i);
        }
        (e, r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use crate::transcript::transcript::Transcript;
    use crate::transcript::transcript_var::TranscriptVar;
    use ark_ff::UniformRand;
    use ark_r1cs_std::prelude::*;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
    use rand::{thread_rng, Rng};

    // Mock implementation of CompressedUniPoly for test purposes
    impl<F: PrimeField + Absorb> CompressedUniPoly<F> {
        pub fn random<R: Rng>(rng: &mut R, degree_bound: usize) -> Self {
            // Generate a random univariate polynomial with the given degree bound
            let coeffs_except_linear_term: Vec<F> = (0..degree_bound).map(|_| F::rand(rng)).collect();
            CompressedUniPoly { coeffs_except_linear_term }
        }
    }

    type F = ScalarField;

    #[test]
    fn test_sumcheck_circuit() {
        // Set up the random number generator
        let mut rng = thread_rng();
        let label = b"test label";

        // Parameters for the test
        let num_rounds = 5;
        let degree_bound = 3;

        // Create a claim as a random field element
        let claim = F::rand(&mut rng);

        // Create random compressed univariate polynomials for the Sumcheck proof
        let compressed_polys: Vec<CompressedUniPoly<F>> = (0..num_rounds)
            .map(|_| CompressedUniPoly::random(&mut rng, degree_bound))
            .collect();

        // Initialize a fresh transcript for both circuits
        let transcript = Transcript::<F>::new(label);

        // Initialize the SumcheckCircuit with random values
        let mut sumcheck_circuit = SumcheckCircuit {
            compressed_polys: compressed_polys.clone(),
            claim: claim.clone(),
            num_rounds,
            degree_bound,
            transcript,
        };

        // Run the verification on the normal SumcheckCircuit
        let (e1, r1) = sumcheck_circuit.verify::<E>();

        // Initialize the constraint system
        let cs: ConstraintSystemRef<F> = ConstraintSystem::new_ref();

        // Create random compressed univariate polynomial vars for the SumCheckCircuitVar
        let compressed_polys_var: Vec<CompressedUniPolyVar<F>> = compressed_polys
            .iter()
            .map(|poly| CompressedUniPolyVar::new(cs.clone(), poly.clone(), AllocationMode::Witness))
            .collect();

        // Initialize the claim variable as a witness
        let claim_var = FpVar::new_witness(cs.clone(), || Ok(claim)).unwrap();

        // Initialize a fresh variable transcript for the variable circuit
        let transcript_var = TranscriptVar::<F>::new(cs.clone(), label);

        // Initialize the SumCheckCircuitVar with the same values
        let mut sumcheck_circuit_var = SumCheckCircuitVar {
            compressed_polys: compressed_polys_var,
            claim: claim_var,
            num_rounds,
            degree_bound,
            transcript: transcript_var,
        };

        // Run the verification on the SumCheckCircuitVar
        let (e2, r2) = sumcheck_circuit_var.verify();

        // Compare the results between the normal and variable circuits
        assert_eq!(e1, e2.value().unwrap());
        assert_eq!(r1.len(), r2.len());

        // check equality element wise
        for i in 0..r1.len() {
            assert_eq!(r1[i], r2[i].value().unwrap());
        }

        // Ensure that the constraint system is satisfied
        assert!(cs.is_satisfied().unwrap());
    }
}

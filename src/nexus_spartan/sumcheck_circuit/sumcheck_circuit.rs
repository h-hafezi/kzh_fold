use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::nexus_spartan::unipoly::unipoly::CompressedUniPoly;
use crate::nexus_spartan::unipoly::unipoly_var::{CompressedUniPolyVar, UniPolyVar};
use crate::transcript::transcript::Transcript;
use crate::transcript::transcript_var::{AppendToTranscriptVar, TranscriptVar};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use std::borrow::Borrow;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SumcheckCircuit<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPoly<F>>,
    pub claim: F,
    pub num_rounds: usize,
    pub degree_bound: usize,
}

impl<F: PrimeField + Absorb> SumcheckCircuit<F> {
    pub fn verify<E: Pairing<ScalarField=F>>(&self, transcript: &mut Transcript<F>) -> (F, Vec<F>) {
        let proof = SumcheckInstanceProof::new(self.compressed_polys.clone());
        proof.verify::<E>(self.claim, self.num_rounds, self.degree_bound, transcript).unwrap()
    }
}

pub struct SumcheckCircuitVar<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPolyVar<F>>,
    claim: FpVar<F>,
    num_rounds: usize,
    degree_bound: usize,
}

// Implement the R1CSVar trait for SumcheckCircuitVar
impl<F: PrimeField + Absorb> R1CSVar<F> for SumcheckCircuitVar<F> {
    type Value = SumcheckCircuit<F>;

    fn cs(&self) -> ConstraintSystemRef<F> {
        // Initialize a None constraint system
        let mut result = ConstraintSystemRef::None;

        // Combine the constraint systems of all compressed polynomials
        for poly_var in &self.compressed_polys {
            result = poly_var.cs().or(result);
        }

        // Combine with the constraint system of the claim
        result = self.claim.cs().or(result);

        result
    }

    // Returns the underlying value of the sumcheck circuit variables
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        // Collect the values of the compressed polynomials, unwrapping each value
        let compressed_poly_values = self
            .compressed_polys
            .iter()
            .map(|poly_var| poly_var.value().unwrap())
            .collect::<Vec<_>>();

        // Get the value of the claim, unwrapping it
        let claim_value = self.claim.value().unwrap();

        // Use the num_rounds and degree_bound directly as they are constants
        let num_rounds = self.num_rounds;
        let degree_bound = self.degree_bound;

        // Return the corresponding SumcheckCircuit value
        Ok(SumcheckCircuit {
            compressed_polys: compressed_poly_values,
            claim: claim_value,
            num_rounds,
            degree_bound,
        })
    }
}

impl<F: PrimeField + Absorb> AllocVar<SumcheckCircuit<F>, F> for SumcheckCircuitVar<F> {
    fn new_variable<T: Borrow<SumcheckCircuit<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        // Retrieve the input function's value
        let binding = f()?;
        let sumcheck_circuit = binding.borrow();

        // Allocate the `compressed_polys` vector by individually allocating each element
        let compressed_polys_var = sumcheck_circuit
            .compressed_polys
            .iter()
            .map(|poly| CompressedUniPolyVar::new_variable(cs.clone(), || Ok(poly), mode).unwrap())
            .collect::<Vec<_>>();

        let claim_var = FpVar::new_variable(
            cs.clone(),
            || Ok(sumcheck_circuit.claim),
            mode
        ).unwrap();

        // Directly set the `num_rounds` and `degree_bound` without allocation
        let num_rounds = sumcheck_circuit.num_rounds;
        let degree_bound = sumcheck_circuit.degree_bound;

        // Return the newly created SumcheckCircuitVar
        Ok(SumcheckCircuitVar {
            compressed_polys: compressed_polys_var,
            claim: claim_var,
            num_rounds,
            degree_bound,
        })
    }
}


impl<F: PrimeField + Absorb> SumcheckCircuitVar<F> {
    pub fn verify(&mut self, transcript: &mut TranscriptVar<F>) -> (FpVar<F>, Vec<FpVar<F>>) {
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
            UniPolyVar::append_to_transcript(&poly, b"poly", transcript);

            //derive the verifier's challenge for the next round
            let r_i = TranscriptVar::challenge_scalar(transcript, b"challenge_nextround");

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
        let mut transcript = Transcript::<F>::new(label);

        // Initialize the SumcheckCircuit with random values
        let mut sumcheck_circuit = SumcheckCircuit {
            compressed_polys: compressed_polys.clone(),
            claim: claim.clone(),
            num_rounds,
            degree_bound,
        };

        // Run the verification on the normal SumcheckCircuit
        let (e1, r1) = sumcheck_circuit.verify::<E>(&mut transcript);

        // Initialize the constraint system
        let cs: ConstraintSystemRef<F> = ConstraintSystem::new_ref();

        // Initialize a fresh variable transcript for the variable circuit
        let mut transcript_var = TranscriptVar::<F>::new(cs.clone(), label);

        // Initialize the SumCheckCircuitVar with the same values
        let mut sumcheck_circuit_var = SumcheckCircuitVar::new_variable(
            cs.clone(),
            || Ok(sumcheck_circuit.clone()),
            AllocationMode::Witness,
        ).unwrap();

        assert_eq!(sumcheck_circuit, sumcheck_circuit_var.value().unwrap());

        // Run the verification on the SumCheckCircuitVar
        let (e2, r2) = sumcheck_circuit_var.verify(&mut transcript_var);

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

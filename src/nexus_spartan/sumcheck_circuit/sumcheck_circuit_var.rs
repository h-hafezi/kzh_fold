use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
use crate::nexus_spartan::unipoly::unipoly_var::{CompressedUniPolyVar, UniPolyVar};
use crate::transcript::transcript_var::{AppendToTranscriptVar, TranscriptVar};
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

pub struct SumcheckCircuitVar<F: PrimeField + Absorb> {
    pub compressed_polys: Vec<CompressedUniPolyVar<F>>,
    pub claim: FpVar<F>,
    pub num_rounds: usize,
    pub degree_bound: usize,
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
            mode,
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
    pub fn verify(&self, transcript: &mut TranscriptVar<F>) -> (FpVar<F>, Vec<FpVar<F>>) {
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

use std::borrow::Borrow;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_relations::r1cs::{Namespace, SynthesisError};
use crate::nexus_spartan::partial_verifier::partial_verifier::PartialVerifier;
use crate::nexus_spartan::sparse_mlpoly::{SparsePoly, SparsePolyVar};
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::{SumcheckCircuitVar, SumcheckCircuit};
use crate::transcript::transcript::Transcript;
use crate::transcript::transcript_var::TranscriptVar;

pub struct PartialVerifierVar<F: PrimeField + Absorb> {
    /// io input, equivalent with
    /// let CRR1CSInstance { input: _input, comm_W, } = instance;
    pub input: Vec<FpVar<F>>,
    /// Sumcheck proof for the polynomial g(x) = \sum eq(tau,x) * (~Az~(x) * ~Bz~(x) - u * ~Cz~(x) - ~E~(x))
    pub sc_proof_phase1: SumcheckCircuitVar<F>,
    /// Evaluation claims for ~Az~(rx), ~Bz~(rx), and ~Cz~(rx).
    pub claims_phase2: (FpVar<F>, FpVar<F>, FpVar<F>),
    /// Sumcheck proof for the polynomial F(x) = ~Z(x)~ * ~ABC~(x), where ABC(x) = \sum_t ~M~(t,x) eq(r,t)
    /// for M a random linear combination of A, B, and C.
    pub sc_proof_phase2: SumcheckCircuitVar<F>,
    /// The claimed evaluation ~Z~(ry)
    pub eval_vars_at_ry: FpVar<F>,
    /// the transcript
    pub transcript: TranscriptVar<F>,
    /// matrix evaluations
    pub evals: (FpVar<F>, FpVar<F>, FpVar<F>),
    /// shape
    pub num_vars: usize,
    pub num_cons: usize,
}

/*
impl<F: PrimeField + Absorb> AllocVar<PartialVerifier<F>, F> for PartialVerifierVar<F> {
    fn new_variable<T: Borrow<PartialVerifier<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        // Convert to Namespace<F>
        let ns = cs.into();
        // Get the constraint system reference
        let cs = ns.cs();

        // Fetch the instance of `PartialVerifier<F>`
        let partial_verifier = f()?.borrow();

        // Allocate the `input` vector of FpVars
        let input = Vec::<FpVar<F>>::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.input),
            mode
        )?;

        // Allocate the sumcheck proof phase 1
        let sc_proof_phase1 = SumCheckCircuitVar::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.sc_proof_phase1),
            mode
        )?;

        // Allocate the claims_phase2 tuple
        let claims_phase2 = (
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.claims_phase2.0), mode)?,
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.claims_phase2.1), mode)?,
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.claims_phase2.2), mode)?
        );

        // Allocate the sumcheck proof phase 2
        let sc_proof_phase2 = SumCheckCircuitVar::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.sc_proof_phase2),
            mode
        )?;

        // Allocate the evaluation of variables at ry
        let eval_vars_at_ry = FpVar::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.eval_vars_at_ry),
            mode
        )?;

        // Allocate the transcript
        let transcript = TranscriptVar::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.transcript),
            mode
        )?;

        // Allocate matrix evaluations
        let evals = (
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.evals.0), mode)?,
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.evals.1), mode)?,
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.evals.2), mode)?
        );

        // Create the final PartialVerifierVar instance
        Ok(PartialVerifierVar {
            input,
            sc_proof_phase1,
            claims_phase2,
            sc_proof_phase2,
            eval_vars_at_ry,
            transcript,
            evals,
            num_vars: partial_verifier.num_vars,
            num_cons: partial_verifier.num_cons,
        })
    }
}
 */

use crate::nexus_spartan::partial_verifier::partial_verifier::PartialVerifier;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuitVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

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
    /// matrix evaluations
    pub evals: (FpVar<F>, FpVar<F>, FpVar<F>),
    /// shape
    pub num_vars: usize,
    pub num_cons: usize,
}


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
        let binding = f()?;
        let partial_verifier = binding.borrow();

        // Allocate the `input` vector of FpVars
        let input = Vec::<FpVar<F>>::new_variable(
            cs.clone(),
            || Ok(partial_verifier.input.clone()),
            mode,
        )?;

        // Allocate the sumcheck proof phase 1
        let sc_proof_phase1 = SumcheckCircuitVar::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.sc_proof_phase1),
            mode,
        )?;

        // Allocate the claims_phase2 tuple
        let claims_phase2 = (
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.claims_phase2.0), mode)?,
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.claims_phase2.1), mode)?,
            FpVar::new_variable(cs.clone(), || Ok(&partial_verifier.claims_phase2.2), mode)?
        );

        // Allocate the sumcheck proof phase 2
        let sc_proof_phase2 = SumcheckCircuitVar::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.sc_proof_phase2),
            mode,
        )?;

        // Allocate the evaluation of variables at ry
        let eval_vars_at_ry = FpVar::new_variable(
            cs.clone(),
            || Ok(&partial_verifier.eval_vars_at_ry),
            mode,
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
            evals,
            num_vars: partial_verifier.num_vars,
            num_cons: partial_verifier.num_cons,
        })
    }
}

// Implement the R1CSVar trait for PartialVerifierVar
impl<F: PrimeField + Absorb> R1CSVar<F> for PartialVerifierVar<F> {
    type Value = PartialVerifier<F>;

    // Combine the constraint systems of all the components
    fn cs(&self) -> ConstraintSystemRef<F> {
        let mut cs_ref = ConstraintSystemRef::None;

        // Combine the constraint system from input vector
        for input_var in &self.input {
            cs_ref = input_var.cs().or(cs_ref);
        }

        // Combine the constraint systems of both sumcheck circuits
        cs_ref = self.sc_proof_phase1.cs().or(cs_ref);
        cs_ref = self.sc_proof_phase2.cs().or(cs_ref);

        // Combine the constraint systems of the evaluation claims (Az, Bz, Cz)
        cs_ref = self.claims_phase2.0.cs().or(cs_ref);
        cs_ref = self.claims_phase2.1.cs().or(cs_ref);
        cs_ref = self.claims_phase2.2.cs().or(cs_ref);

        // Combine the constraint system of eval_vars_at_ry
        cs_ref = self.eval_vars_at_ry.cs().or(cs_ref);

        // Combine the constraint systems of the evals tuple (evaluations)
        cs_ref = self.evals.0.cs().or(cs_ref);
        cs_ref = self.evals.1.cs().or(cs_ref);
        cs_ref = self.evals.2.cs().or(cs_ref);

        cs_ref
    }

    // Extract the value of all components
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let input_val = self.input.iter().map(|var| var.value().unwrap()).collect();
        let sc_proof_phase1_val = self.sc_proof_phase1.value()?;
        let claims_phase2_val = (
            self.claims_phase2.0.value()?,
            self.claims_phase2.1.value()?,
            self.claims_phase2.2.value()?,
        );
        let sc_proof_phase2_val = self.sc_proof_phase2.value()?;
        let eval_vars_at_ry_val = self.eval_vars_at_ry.value()?;
        let evals_val = (
            self.evals.0.value()?,
            self.evals.1.value()?,
            self.evals.2.value()?,
        );

        // Return the full PartialVerifierVar value
        Ok(PartialVerifier {
            input: input_val,
            sc_proof_phase1: sc_proof_phase1_val,
            claims_phase2: claims_phase2_val,
            sc_proof_phase2: sc_proof_phase2_val,
            eval_vars_at_ry: eval_vars_at_ry_val,
            evals: evals_val,
            num_vars: self.num_vars,
            num_cons: self.num_cons,
        })
    }
}


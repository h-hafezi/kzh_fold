use crate::accumulation_circuit::affine_to_projective;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::math::Math;
use crate::nexus_spartan::partial_verifier::partial_verifier::SpartanPartialVerifier;
use crate::nexus_spartan::sparse_polynomial::sparse_polynomial_var::SparsePolyVar;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit_var::SumcheckCircuitVar;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

pub struct PartialVerifierVar<F: PrimeField + Absorb, G1>
where
    G1: SWCurveConfig<ScalarField=F> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
{
    /// io input, equivalent with, CRR1CSInstance { input: _input, comm_W, } = instance;
    pub instance: (Vec<FpVar<F>>, NonNativeAffineVar<G1>),
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

impl<F: PrimeField + Absorb, G1> PartialVerifierVar<F, G1>
where
    G1: SWCurveConfig<ScalarField=F> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
{
    pub fn verify(&self, transcript: &mut TranscriptVar<F>) -> (Vec<FpVar<F>>, Vec<FpVar<F>>) {
        TranscriptVar::append_scalars(transcript, b"input", self.instance.0.as_slice());
        TranscriptVar::append_scalars_non_native(
            transcript,
            b"input",
            &[self.instance.1.x.clone(), self.instance.1.y.clone()],
        );

        let n = self.num_vars;

        let (num_rounds_x, num_rounds_y) = (self.num_cons.log_2(), (2 * self.num_vars).log_2());

        // check number of round to be consistent with sc_proof_phase1 and sc_proof_phase2
        assert_eq!(self.sc_proof_phase1.num_rounds, num_rounds_x);
        assert_eq!(self.sc_proof_phase2.num_rounds, num_rounds_y);

        // derive the verifier's challenge tau
        let tau = TranscriptVar::challenge_vector(
            transcript,
            b"challenge_tau",
            num_rounds_x,
        );

        // consistency check for sc_proof_phase1
        assert_eq!(self.sc_proof_phase1.degree_bound, 3);
        self.sc_proof_phase1.claim.enforce_equal(&FpVar::zero()).expect("equality error");

        let (claim_post_phase1, rx) = self.sc_proof_phase1.verify(transcript);

        // perform the intermediate sum-check test with claimed Az, Bz, Cz, and E
        let (Az_claim, Bz_claim, Cz_claim) = &self.claims_phase2;

        TranscriptVar::append_scalar(transcript, b"Az_claim", &Az_claim);
        TranscriptVar::append_scalar(transcript, b"Bz_claim", &Bz_claim);
        TranscriptVar::append_scalar(transcript, b"Cz_claim", &Cz_claim);

        let mut taus_bound_rx: FpVar<F> = FpVar::one();
        for i in 0..rx.len() {
            taus_bound_rx *= &rx[i] * &tau[i] + (FpVar::one() - &rx[i]) * (FpVar::one() - &tau[i]);
        }

        let expected_claim_post_phase1 = (Az_claim * Bz_claim - Cz_claim) * taus_bound_rx;
        expected_claim_post_phase1.enforce_equal(&claim_post_phase1).expect("equality error");

        // derive three public challenges and then derive a joint claim
        let r_A = TranscriptVar::challenge_scalar(transcript, b"challenege_Az");
        let r_B = TranscriptVar::challenge_scalar(transcript, b"challenege_Bz");
        let r_C = TranscriptVar::challenge_scalar(transcript, b"challenege_Cz");

        // r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;
        let claim_phase2 = &r_A * Az_claim + &r_B * Bz_claim + &r_C * Cz_claim;

        // consistency check for sc_proof_phase1
        assert_eq!(self.sc_proof_phase2.degree_bound, 2);
        self.sc_proof_phase2.claim.enforce_equal(&claim_phase2).expect("equality error");

        let (claim_post_phase2, ry) = self.sc_proof_phase2.verify(transcript);

        // Compute (1, io)(r_y) so that we can use it to compute Z(r_y)
        let poly_input_eval = {
            // constant term: one
            let mut input_as_sparse_poly_entries = vec![FpVar::one()];
            // remaining inputs:
            input_as_sparse_poly_entries.extend(
                (0..self.instance.0.len())
                    .map(|i| self.instance.0[i].clone())
                    .collect::<Vec<FpVar<F>>>(),
            );
            SparsePolyVar::new(n.log_2(), input_as_sparse_poly_entries).evaluate(&ry[1..])
        };

        // compute Z(r_y): eval_Z_at_ry = (F::one() - ry[0]) * self.eval_vars_at_ry + ry[0] * poly_input_eval
        let eval_Z_at_ry = (FpVar::one() - &ry[0]) * &self.eval_vars_at_ry + &ry[0] * poly_input_eval;

        // perform the final check in the second sum-check protocol
        let (eval_A_r, eval_B_r, eval_C_r) = &self.evals;
        let expected_claim_post_phase2 = eval_Z_at_ry * (&r_A * eval_A_r + &r_B * eval_B_r + &r_C * eval_C_r);

        expected_claim_post_phase2.enforce_equal(&claim_post_phase2).expect("equality error");

        (rx, ry)
    }
}

impl<F: PrimeField + Absorb, G1, E> AllocVar<SpartanPartialVerifier<F, E>, F> for PartialVerifierVar<F, G1>
where
    G1: SWCurveConfig<ScalarField=F> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
{
    fn new_variable<T: Borrow<SpartanPartialVerifier<F, E>>>(
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
            || Ok(partial_verifier.instance.0.clone()),
            mode,
        )?;

        let com_w = NonNativeAffineVar::new_variable(
            ns!(cs, "C"),
            || Ok(affine_to_projective(partial_verifier.instance.1.clone())),
            mode,
        ).unwrap();


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
            instance: (input, com_w),
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

#[cfg(test)]
mod tests {
    use ark_relations::r1cs::ConstraintSystem;

    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use crate::nexus_spartan::partial_verifier::partial_verifier::tests::partial_verifier_test_helper;

    #[test]
    pub fn test_partial_verifier_circuit() {
        let (partial_verifier, transcript) = partial_verifier_test_helper::<E, MultilinearPolynomial<ScalarField>, ScalarField>();
        let cs = ConstraintSystem::<ScalarField>::new_ref();
        let partial_verifier_var = PartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(partial_verifier.clone()),
            AllocationMode::Input,
        ).unwrap();

        let mut transcript_var = TranscriptVar::from_transcript(cs.clone(), transcript);

        let (_r_x, _r_y) = partial_verifier_var.verify(&mut transcript_var);
        println!("constraint count: {} {}", cs.num_instance_variables(), cs.num_witness_variables());
        assert!(cs.is_satisfied().unwrap());
    }
}


use crate::math::Math;
use crate::nexus_spartan::crr1cs::is_sat;
use crate::nexus_spartan::crr1csproof::CRR1CSProof;
use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
use crate::nexus_spartan::sparse_polynomial::sparse_polynomial::SparsePoly;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ff::PrimeField;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SpartanPartialVerifier<F, E>
where
    F: PrimeField + Absorb,
    E: Pairing<ScalarField=F>,
{
    /// io input, equivalent with
    /// let CRR1CSInstance { input: _input, comm_W, } = instance;
    pub instance: (Vec<F>, E::G1Affine),
    /// Sumcheck proof for the polynomial g(x) = \sum eq(tau,x) * (~Az~(x) * ~Bz~(x) - u * ~Cz~(x) - ~E~(x))
    pub sc_proof_phase1: SumcheckCircuit<F>,
    /// Evaluation claims for ~Az~(rx), ~Bz~(rx), and ~Cz~(rx).
    pub claims_phase2: (F, F, F),
    /// Sumcheck proof for the polynomial F(x) = ~Z(x)~ * ~ABC~(x), where ABC(x) = \sum_t ~M~(t,x) eq(r,t)
    /// for M a random linear combination of A, B, and C.
    pub sc_proof_phase2: SumcheckCircuit<F>,
    /// The claimed evaluation ~Z~(ry)
    pub eval_vars_at_ry: F,
    /// matrix evaluations
    pub evals: (F, F, F),
    /// shape
    pub num_vars: usize,
    pub num_cons: usize,
}

impl<F: PrimeField + Absorb, E: Pairing<ScalarField=F>> SpartanPartialVerifier<F, E>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    pub fn initialise<PC: PolyCommitmentScheme<E>>(
        proof: &CRR1CSProof<E, PC, F>,
        num_vars: usize,
        num_cons: usize,
        instance: (Vec<F>, E::G1Affine),
        evals: &(F, F, F),
        transcript: &mut Transcript<F>,
    ) -> Self {
        Transcript::append_scalars(transcript, b"input", instance.0.as_slice());
        Transcript::append_point::<E>(transcript, b"witness", &instance.1);

        let n = num_vars;

        let (num_rounds_x, num_rounds_y) = (num_cons.log_2(), (2 * num_vars).log_2());

        // derive the verifier's challenge tau
        let tau = Transcript::challenge_vector(
            transcript,
            b"challenge_tau",
            num_rounds_x,
        );

        // verify the first sum-check instance
        let claim_phase1 = F::zero();
        let (claim_post_phase1, rx) = proof
            .sc_proof_phase1
            .verify::<E>(claim_phase1, num_rounds_x, 3, transcript).unwrap();

        let sc_proof_phase1_circuit = SumcheckCircuit {
            compressed_polys: proof.sc_proof_phase1.compressed_polys.clone(),
            claim: claim_phase1,
            num_rounds: num_rounds_x,
            degree_bound: 3,
        };

        // perform the intermediate sum-check test with claimed Az, Bz, Cz, and E
        let (Az_claim, Bz_claim, Cz_claim) = proof.claims_phase2;

        Transcript::append_scalar(transcript, b"Az_claim", &Az_claim);
        Transcript::append_scalar(transcript, b"Bz_claim", &Bz_claim);
        Transcript::append_scalar(transcript, b"Cz_claim", &Cz_claim);

        let taus_bound_rx: F = (0..rx.len())
            .map(|i| rx[i] * tau[i] + (F::one() - rx[i]) * (F::one() - tau[i]))
            .product();

        let expected_claim_post_phase1 = (Az_claim * Bz_claim - Cz_claim) * taus_bound_rx;
        assert_eq!(expected_claim_post_phase1, claim_post_phase1);

        // derive three public challenges and then derive a joint claim
        let r_A = Transcript::challenge_scalar(transcript, b"challenege_Az");
        let r_B = Transcript::challenge_scalar(transcript, b"challenege_Bz");
        let r_C = Transcript::challenge_scalar(transcript, b"challenege_Cz");

        // r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;
        let claim_phase2 = r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;

        // verify the joint claim with a sum-check protocol
        let (claim_post_phase2, ry) = proof
            .sc_proof_phase2
            .verify::<E>(claim_phase2, num_rounds_y, 2, transcript).unwrap();

        let sc_proof_phase2_circuit = SumcheckCircuit {
            compressed_polys: proof.sc_proof_phase2.compressed_polys.clone(),
            claim: claim_phase2,
            num_rounds: num_rounds_y,
            degree_bound: 2,
        };

        // Compute (1,io)(r_y) so that we can use it to compute Z(r_y)
        let poly_input_eval = {
            let mut input_as_sparse_poly_entries = vec![F::ONE];
            input_as_sparse_poly_entries.extend(
                (0..instance.0.len())
                    .map(|i| instance.0[i])
                    .collect::<Vec<F>>(),
            );
            SparsePoly::new(n.log_2(), input_as_sparse_poly_entries).evaluate(&ry[1..])
        };

        // compute Z(r_y): eval_Z_at_ry = (F::one() - ry[0]) * self.eval_vars_at_ry + ry[0] * poly_input_eval
        let eval_Z_at_ry = (F::one() - ry[0]) * proof.eval_vars_at_ry + ry[0] * poly_input_eval;

        // perform the final check in the second sum-check protocol
        let (eval_A_r, eval_B_r, eval_C_r) = evals;
        let expected_claim_post_phase2 = eval_Z_at_ry * (r_A * eval_A_r + r_B * eval_B_r + r_C * eval_C_r);

        assert_eq!(expected_claim_post_phase2, claim_post_phase2);

        SpartanPartialVerifier {
            instance,
            sc_proof_phase1: sc_proof_phase1_circuit,
            claims_phase2: proof.claims_phase2,
            sc_proof_phase2: sc_proof_phase2_circuit,
            eval_vars_at_ry: proof.eval_vars_at_ry,
            evals: *evals,
            num_vars,
            num_cons,
        }
    }


    pub fn verify(&self, transcript: &mut Transcript<F>) -> (Vec<F>, Vec<F>) {
        Transcript::append_scalars(transcript, b"input", self.instance.0.as_slice());
        Transcript::append_point::<E>(transcript, b"witness", &self.instance.1);

        let n = self.num_vars;

        let (num_rounds_x, num_rounds_y) = (self.num_cons.log_2(), (2 * self.num_vars).log_2());

        // check number of round to be consistent with sc_proof_phase1 and sc_proof_phase2
        assert_eq!(self.sc_proof_phase1.num_rounds, num_rounds_x);
        assert_eq!(self.sc_proof_phase2.num_rounds, num_rounds_y);

        // derive the verifier's challenge tau
        let tau = Transcript::challenge_vector(
            transcript,
            b"challenge_tau",
            num_rounds_x,
        );

        // consistency check for sc_proof_phase1
        assert_eq!(self.sc_proof_phase1.degree_bound, 3);
        assert_eq!(self.sc_proof_phase1.claim, F::zero());
        let (claim_post_phase1, rx) = self.sc_proof_phase1.verify::<E>(transcript);

        // perform the intermediate sum-check test with claimed Az, Bz, Cz, and E
        let (Az_claim, Bz_claim, Cz_claim) = self.claims_phase2;

        Transcript::append_scalar(transcript, b"Az_claim", &Az_claim);
        Transcript::append_scalar(transcript, b"Bz_claim", &Bz_claim);
        Transcript::append_scalar(transcript, b"Cz_claim", &Cz_claim);

        let taus_bound_rx: F = (0..rx.len())
            .map(|i| rx[i] * tau[i] + (F::one() - rx[i]) * (F::one() - tau[i]))
            .product();

        let expected_claim_post_phase1 = (Az_claim * Bz_claim - Cz_claim) * taus_bound_rx;
        assert_eq!(expected_claim_post_phase1, claim_post_phase1);

        // derive three public challenges and then derive a joint claim
        let r_A = Transcript::challenge_scalar(transcript, b"challenege_Az");
        let r_B = Transcript::challenge_scalar(transcript, b"challenege_Bz");
        let r_C = Transcript::challenge_scalar(transcript, b"challenege_Cz");

        // r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;
        let claim_phase2 = r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;

        // consistency check for sc_proof_phase1
        assert_eq!(self.sc_proof_phase2.degree_bound, 2);
        assert_eq!(self.sc_proof_phase2.claim, claim_phase2);
        let (claim_post_phase2, ry) = self.sc_proof_phase2.verify::<E>(transcript);

        // Compute (1, io)(r_y) so that we can use it to compute Z(r_y)
        let poly_input_eval = {
            // constant term: one
            let mut input_as_sparse_poly_entries = vec![F::ONE];
            // remaining inputs:
            input_as_sparse_poly_entries.extend(
                (0..self.instance.0.len())
                    .map(|i| self.instance.0[i])
                    .collect::<Vec<F>>(),
            );
            SparsePoly::new(n.log_2(), input_as_sparse_poly_entries).evaluate(&ry[1..])
        };

        // compute Z(r_y): eval_Z_at_ry = (F::one() - ry[0]) * self.eval_vars_at_ry + ry[0] * poly_input_eval
        let eval_Z_at_ry = (F::one() - ry[0]) * self.eval_vars_at_ry + ry[0] * poly_input_eval;

        // perform the final check in the second sum-check protocol
        let (eval_A_r, eval_B_r, eval_C_r) = self.evals;
        let expected_claim_post_phase2 = eval_Z_at_ry * (r_A * eval_A_r + r_B * eval_B_r + r_C * eval_C_r);

        assert_eq!(expected_claim_post_phase2, claim_post_phase2);

        (rx, ry)
    }
}

#[cfg(test)]
pub mod tests {
    use crate::nexus_spartan::crr1cs::produce_synthetic_crr1cs;

    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use crate::nexus_spartan::polycommitments::ToAffine;

    #[test]
    pub fn check_verification_proof() {
        partial_verifier_test_helper::<E, MultilinearPolynomial<ScalarField>, ScalarField>();
    }

    pub fn partial_verifier_test_helper<E, PC, F>() -> (SpartanPartialVerifier<F, E>, Transcript<F>)
    where
        F: PrimeField + Absorb,
        PC: PolyCommitmentScheme<E>,
        E: Pairing<ScalarField=F>,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        let num_vars = 1024;
        let num_cons = num_vars;
        let num_inputs = 10;

        // this generates a new instance/witness for spartan as well as PCS parameters
        let (shape, instance, witness, gens) = produce_synthetic_crr1cs::<E, PC>(num_cons, num_vars, num_inputs);
        assert!(is_sat(&shape, &instance, &witness, &gens.gens_r1cs_sat).unwrap());

        let (num_cons, num_vars, _num_inputs) = (
            shape.get_num_cons(),
            shape.get_num_vars(),
            shape.get_num_inputs(),
        );

        let mut prover_transcript = Transcript::new(b"example");

        let (proof, rx, ry) = CRR1CSProof::prove(
            &shape,
            &instance,
            witness,
            &gens.gens_r1cs_sat,
            &mut prover_transcript,
        );

        let inst_evals = shape.inst.inst.evaluate(&rx, &ry);

        let mut prover_transcript = Transcript::new(b"example");
        let mut verifier_transcript = prover_transcript.clone();
        let verifier_transcript_clone = verifier_transcript.clone();
        let partial_verifier = SpartanPartialVerifier::initialise(
            &proof,
            num_vars,
            num_cons,
            (instance.input.assignment, {
                let com_w: E::G1Affine = instance.comm_W.to_affine();
                com_w
            }),
            &inst_evals,
            &mut prover_transcript,
        );

        partial_verifier.verify(&mut verifier_transcript);

        (partial_verifier, verifier_transcript_clone)
    }
}

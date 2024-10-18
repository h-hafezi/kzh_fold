use crate::math::Math;
use crate::nexus_spartan::sparse_mlpoly::{SparsePolyEntry, SparsePolynomial};
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use crate::nexus_spartan::crr1cs::is_sat;
use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
use crate::polynomial::multilinear_poly::MultilinearPolynomial;

pub struct PartialVerifierCircuit<F: PrimeField + Absorb> {
    /// io input, equivalent with
    /// let CRR1CSInstance { input: _input, comm_W, } = instance;
    /// let input = _input.assignment.as_slice();
    /// in crr1csproof.rs
    input: Vec<F>,
    /// Sumcheck proof for the polynomial g(x) = \sum eq(tau,x) * (~Az~(x) * ~Bz~(x) - u * ~Cz~(x) - ~E~(x))
    sc_proof_phase1: SumcheckCircuit<F>,
    /// Evaluation claims for ~Az~(rx), ~Bz~(rx), and ~Cz~(rx).
    claims_phase2: (F, F, F),
    /// Sumcheck proof for the polynomial F(x) = ~Z(x)~ * ~ABC~(x), where ABC(x) = \sum_t ~M~(t,x) eq(r,t)
    /// for M a random linear combination of A, B, and C.
    sc_proof_phase2: SumcheckCircuit<F>,
    /// The claimed evaluation ~Z~(ry)
    eval_vars_at_ry: F,
    /// the transcript
    transcript: Transcript<F>,
    /// matrix evaluations
    evals: (F, F, F),
    /// shape
    num_vars: usize,
    num_cons: usize,
}

impl<F: PrimeField + Absorb> PartialVerifierCircuit<F> {
    pub fn verify<E: Pairing<ScalarField=F>>(&mut self) -> (Vec<F>, Vec<F>) {
        Transcript::append_scalars(&mut self.transcript, b"input", self.input.as_slice());
        let n = self.num_vars;

        let (num_rounds_x, num_rounds_y) = (self.num_cons.log_2(), (2 * self.num_vars).log_2());

        // check number of round to be consistent with sc_proof_phase1 and sc_proof_phase2
        assert_eq!(self.sc_proof_phase1.num_rounds, num_rounds_x);
        assert_eq!(self.sc_proof_phase2.num_rounds, num_rounds_y);

        // derive the verifier's challenge tau
        let tau = Transcript::challenge_vector(
            &mut self.transcript,
            b"challenge_tau",
            num_rounds_x,
        );

        // consistency check for sc_proof_phase1
        assert_eq!(self.sc_proof_phase1.transcript.state, self.transcript.state);
        assert_eq!(self.sc_proof_phase1.degree_bound, 3);
        assert_eq!(self.sc_proof_phase1.claim, F::zero());
        let (claim_post_phase1, rx) = self.sc_proof_phase1.verify::<E>();

        // perform the intermediate sum-check test with claimed Az, Bz, Cz, and E
        let (Az_claim, Bz_claim, Cz_claim) = self.claims_phase2;

        Transcript::append_scalar(&mut self.transcript, b"Az_claim", &Az_claim);
        Transcript::append_scalar(&mut self.transcript, b"Bz_claim", &Bz_claim);
        Transcript::append_scalar(&mut self.transcript, b"Cz_claim", &Cz_claim);

        let taus_bound_rx: F = (0..rx.len())
            .map(|i| rx[i] * tau[i] + (F::one() - rx[i]) * (F::one() - tau[i]))
            .product();

        let expected_claim_post_phase1 = (Az_claim * Bz_claim - Cz_claim) * taus_bound_rx;
        assert_eq!(expected_claim_post_phase1, claim_post_phase1);

        // derive three public challenges and then derive a joint claim
        let r_A = Transcript::challenge_scalar(&mut self.transcript, b"challenege_Az");
        let r_B = Transcript::challenge_scalar(&mut self.transcript, b"challenege_Bz");
        let r_C = Transcript::challenge_scalar(&mut self.transcript, b"challenege_Cz");

        // r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;
        let claim_phase2 = r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;

        // consistency check for sc_proof_phase1
        assert_eq!(self.sc_proof_phase2.transcript.state, self.transcript.state);
        assert_eq!(self.sc_proof_phase2.degree_bound, 2);
        assert_eq!(self.sc_proof_phase2.claim, claim_phase2);
        let (claim_post_phase2, ry) = self.sc_proof_phase2.verify::<E>();

        // Compute (io,1)(r_y) so that we can use it to compute Z(r_y)
        let poly_input_eval = {
            // constant term
            let mut input_as_sparse_poly_entries = vec![SparsePolyEntry::new(0, F::ONE)];
            //remaining inputs
            input_as_sparse_poly_entries.extend(
                (0..self.input.len())
                    .map(|i| SparsePolyEntry::new(i + 1, self.input[i]))
                    .collect::<Vec<SparsePolyEntry<F>>>(),
            );
            SparsePolynomial::new(n.log_2(), input_as_sparse_poly_entries).evaluate(&ry[1..])
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

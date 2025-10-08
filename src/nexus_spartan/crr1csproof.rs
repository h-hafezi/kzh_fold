#![allow(clippy::too_many_arguments)]
#![allow(dead_code)]

use super::unipoly::unipoly::{CompressedUniPoly, UniPoly};
use ark_crypto_primitives::sponge::Absorb;

pub use super::crr1cs::*;
use super::errors::ProofVerifyError;
use super::sumcheck::SumcheckInstanceProof;
use super::timer::Timer;
use crate::math::Math;
use crate::nexus_spartan::sparse_polynomial::sparse_polynomial::SparsePoly;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use rand::thread_rng;
use crate::kzh::KZH;

#[derive(CanonicalSerialize, CanonicalDeserialize, Debug)]
pub struct CRR1CSProof<E: Pairing<ScalarField=F>, PC: KZH<E>, F: PrimeField + Absorb> {
    /// Sumcheck proof for the polynomial g(x) = \sum eq(tau,x) * (~Az~(x) * ~Bz~(x) - u * ~Cz~(x) - ~E~(x))
    pub sc_proof_phase1: SumcheckInstanceProof<F>,
    /// Evaluation claims for ~Az~(rx), ~Bz~(rx), and ~Cz~(rx).
    pub claims_phase2: (F, F, F),
    /// Sumcheck proof for the polynomial F(x) = ~Z(x)~ * ~ABC~(x), where ABC(x) = \sum_t ~M~(t,x) eq(r,t)
    /// for M a random linear combination of A, B, and C.
    pub sc_proof_phase2: SumcheckInstanceProof<F>,
    /// The claimed evaluation ~Z~(ry)
    pub eval_vars_at_ry: F,
    /// A polynomial evaluation proof of the claimed evaluation ~Z~(ry) with respect to the commitment comm_W.
    pub proof_eval_vars_at_ry: PC::Opening,
}

impl<F: PrimeField + Absorb> SumcheckInstanceProof<F> {
    pub fn prove_quad<Func, E>(
        claim: &F,
        num_rounds: usize,
        poly_A: &mut MultilinearPolynomial<F>,
        poly_B: &mut MultilinearPolynomial<F>,
        comb_func: Func,
        transcript: &mut Transcript<F>,
    ) -> (Self, Vec<F>, Vec<F>)
    where
        Func: Fn(&F, &F) -> F,
        E: Pairing<ScalarField=F>,
    {
        let mut e = *claim;
        let mut r: Vec<F> = Vec::new();
        let mut quad_polys: Vec<CompressedUniPoly<F>> = Vec::new();
        assert_eq!(poly_A.num_variables, poly_B.num_variables);
        for _j in 0..num_rounds {
            let mut eval_point_0 = F::zero();
            let mut eval_point_2 = F::zero();

            let len = poly_A.len() / 2;
            for i in 0..len {
                // eval 0: bound_func is A(low)
                eval_point_0 += comb_func(&poly_A[i], &poly_B[i]);

                // eval 2: bound_func is -A(low) + 2*A(high)
                let poly_A_bound_point = poly_A[len + i] + poly_A[len + i] - poly_A[i];
                let poly_B_bound_point = poly_B[len + i] + poly_B[len + i] - poly_B[i];

                eval_point_2 += comb_func(&poly_A_bound_point, &poly_B_bound_point);
            }

            let evals = vec![eval_point_0, e - eval_point_0, eval_point_2];
            let poly = UniPoly::from_evals(&evals);

            // append the prover's message to the transcript
            <UniPoly<F> as AppendToTranscript<F>>::append_to_transcript(&poly, b"poly", transcript);

            //derive the verifier's challenge for the next round
            let r_j = Transcript::challenge_scalar(transcript, b"challenge_nextround");

            r.push(r_j);
            // bound all tables to the verifier's challenege
            poly_A.bound_poly_var_top(&r_j);
            poly_B.bound_poly_var_top(&r_j);
            e = poly.evaluate(&r_j);
            quad_polys.push(poly.compress());
        }

        (
            SumcheckInstanceProof::new(quad_polys),
            r,
            vec![poly_A[0], poly_B[0]],
        )
    }

    pub fn prove_cubic_five_terms<Func, E>(
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
        E: Pairing<ScalarField=F>,
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

            let evals = vec![eval_point_0, e - eval_point_0, eval_point_2, eval_point_3];
            let poly = UniPoly::from_evals(&evals);

            // append the prover's message to the transcript
            <UniPoly<F> as AppendToTranscript<F>>::append_to_transcript(&poly, b"poly", transcript);

            //derive the verifier's challenge for the next round
            let r_j = Transcript::challenge_scalar(transcript, b"challenge_nextround");

            r.push(r_j);
            // bound all tables to the verifier's challenege
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

impl<E: Pairing<ScalarField=F>, PC: KZH<E>, F: PrimeField + Absorb> CRR1CSProof<E, PC, F> {
    #[allow(clippy::type_complexity)]
    /// Generates the sumcheck proof that sum_s evals_tau(s) * (evals_Az(s) * evals_Bz(s) - u * evals_Cz(s) - E(s)) == 0.
    /// Note that this proof does not use blinding factors, so this is not zero-knowledge.
    fn prove_phase_one(
        num_rounds: usize,
        evals_tau: &mut MultilinearPolynomial<F>,
        evals_Az: &mut MultilinearPolynomial<F>,
        evals_Bz: &mut MultilinearPolynomial<F>,
        evals_Cz: &mut MultilinearPolynomial<F>,
        transcript: &mut Transcript<F>,
    ) -> (
        SumcheckInstanceProof<F>,
        Vec<F>,
        Vec<F>,
    )
    {
        let relaxed_comb_func =
            |poly_tau: &F,
             poly_A: &F,
             poly_B: &F,
             poly_C: &F|
             -> F { (*poly_A * *poly_B - *poly_C) * *poly_tau };

        let (sc_proof_phase_one, r, claims) = SumcheckInstanceProof::prove_cubic_five_terms::<_, E>(
            &F::zero(), // claim is zero
            num_rounds,
            evals_tau,
            evals_Az,
            evals_Bz,
            evals_Cz,
            relaxed_comb_func,
            transcript,
        );

        (sc_proof_phase_one, r, claims)
    }

    /// Generates the sumcheck proof that `claim` = sum_{s,t} eq(r, t) evals_ABC(t, s) evals_z(s)
    #[allow(clippy::type_complexity)]
    fn prove_phase_two(
        num_rounds: usize,
        claim: &F,
        evals_z: &mut MultilinearPolynomial<F>,
        evals_ABC: &mut MultilinearPolynomial<F>,
        transcript: &mut Transcript<F>,
    ) -> (
        SumcheckInstanceProof<F>,
        Vec<F>,
        Vec<F>,
    ) {
        let comb_func = |poly_A_comp: &F,
                         poly_B_comp: &F|
                         -> F { poly_A_comp.clone() * poly_B_comp.clone() };
        let (sc_proof_phase_two, r, claims) = SumcheckInstanceProof::prove_quad::<_, E>(
            claim, num_rounds, evals_z, evals_ABC, comb_func, transcript,
        );

        (sc_proof_phase_two, r, claims)
    }

    fn protocol_name() -> &'static [u8] {
        b"CRR1CS proof"
    }

    #[allow(clippy::type_complexity)]
    pub fn prove(
        shape: &CRR1CSShape<F>,
        instance: &CRR1CSInstance<E, PC>,
        witness: CRR1CSWitness<F>,
        srs: &PC::SRS,
        transcript: &mut Transcript<F>,
    ) -> (CRR1CSProof<E, PC, F>, Vec<F>, Vec<F>) {
        let timer_prove = Timer::new("CRR1CSProof::prove");

        let _inst = &shape.inst.inst;
        let CRR1CSInstance {
            input,
            comm_W,
            aux_W,
        } = instance;

        let CRR1CSWitness { W: vars } = witness;

        let (inst, input, vars) = (&_inst, input.assignment.as_slice(), vars.assignment);

        // we currently require the number of |inputs| + 1 to be at most number of vars
        assert!(input.len() < vars.len());
        Transcript::append_scalars(transcript, b"input", input);
        AppendToTranscript::append_to_transcript(comm_W, b"witness", transcript);

        // create a multilinear polynomial using the supplied assignment for variables
        let poly_vars = MultilinearPolynomial::<F>::new(vars.clone());

        let timer_sc_proof_phase1 = Timer::new("prove_sc_phase_one");

        // append input to variables to create a single vector z
        let z = {
            let num_inputs = input.len();
            let num_vars = vars.len();
            let mut z = vars;
            z.extend(vec![F::ONE]); // add relaxed constant term in z
            z.extend(input);
            z.extend(&vec![F::zero(); num_vars - num_inputs - 1]); // we will pad with zeros
            z
        };

        // derive the verifier's challenge tau
        let (num_rounds_x, num_rounds_y) = (inst.get_num_cons().log_2(), z.len().log_2());
        let tau = Transcript::challenge_vector(
            transcript,
            b"challenge_tau",
            num_rounds_x,
        );

        // compute the initial evaluation table for R(\tau, x)
        let (mut poly_Az, mut poly_Bz, mut poly_Cz) =
            inst.multiply_vec(inst.get_num_cons(), z.len(), &z);

        let mut poly_tau = MultilinearPolynomial::new(EqPolynomial::new(tau).evals());
        let (sc_proof_phase1, rx, _claims_phase1) = CRR1CSProof::<E, PC, F>::prove_phase_one(
            num_rounds_x,
            &mut poly_tau,
            &mut poly_Az,
            &mut poly_Bz,
            &mut poly_Cz,
            transcript,
        );

        assert_eq!(poly_tau.len(), 1);
        assert_eq!(poly_Az.len(), 1);
        assert_eq!(poly_Bz.len(), 1);
        assert_eq!(poly_Cz.len(), 1);
        timer_sc_proof_phase1.stop();

        let (_, Az_claim, Bz_claim, Cz_claim) = (
            &poly_tau[0],
            &poly_Az[0],
            &poly_Bz[0],
            &poly_Cz[0],
        );

        Transcript::append_scalar(transcript, b"Az_claim", Az_claim);
        Transcript::append_scalar(transcript, b"Bz_claim", Bz_claim);
        Transcript::append_scalar(transcript, b"Cz_claim", Cz_claim);

        let timer_sc_proof_phase2 = Timer::new("prove_sc_phase_two");

        // combine the three claims into a single claim
        let r_A = Transcript::challenge_scalar(transcript, b"challenege_Az");
        let r_B = Transcript::challenge_scalar(transcript, b"challenege_Bz");
        let r_C = Transcript::challenge_scalar(transcript, b"challenege_Cz");
        let claim_phase2 = r_A * Az_claim + r_B * Bz_claim + r_C * Cz_claim;

        let evals_ABC = {
            // compute the initial evaluation table for R(\tau, x)
            let evals_rx = EqPolynomial::new(rx.clone()).evals();
            let (evals_A, evals_B, evals_C) = inst.compute_eval_table_sparse(inst.get_num_cons(), z.len(), &evals_rx);

            assert_eq!(evals_A.len(), evals_B.len());
            assert_eq!(evals_A.len(), evals_C.len());
            (0..evals_A.len())
                .map(|i| r_A * evals_A[i] + r_B * evals_B[i] + r_C * evals_C[i])
                .collect::<Vec<F>>()
        };

        assert!(z.len().is_power_of_two(), "z must be a power of two");
        assert!(evals_ABC.len().is_power_of_two(), "eval_ABC must be a power of two");
        assert_eq!(z.len(), evals_ABC.len(), "vector z and eval_ABC should have equal length");

        // another instance of the sum-check protocol
        let (sc_proof_phase2, ry, _claims_phase2) = CRR1CSProof::<E, PC, F>::prove_phase_two(
            num_rounds_y,
            &claim_phase2,
            &mut MultilinearPolynomial::new(z),
            &mut MultilinearPolynomial::new(evals_ABC),
            transcript,
        );
        timer_sc_proof_phase2.stop();

        let timer_polyeval = Timer::new("poly eval");
        let eval_vars_at_ry = poly_vars.evaluate(&ry[1..]);
        timer_polyeval.stop();

        let timer_polyevalproof = Timer::new("poly eval proof");
        let proof_eval_vars_at_ry = {
            PC::open(
                &srs,
                &ry[1..],
                comm_W,
                aux_W,
                &poly_vars,
                &mut thread_rng(),
            )
        };

        timer_polyevalproof.stop();

        timer_prove.stop();

        (
            CRR1CSProof {
                sc_proof_phase1,
                claims_phase2: (*Az_claim, *Bz_claim, *Cz_claim),
                sc_proof_phase2,
                eval_vars_at_ry,
                proof_eval_vars_at_ry,
            },
            rx,
            ry,
        )
    }

    #[allow(clippy::type_complexity)]
    pub fn verify(
        &self,
        num_vars: usize,
        num_cons: usize,
        instance: &CRR1CSInstance<E, PC>,
        evals: &(F, F, F),
        transcript: &mut Transcript<F>,
    ) -> Result<(Vec<F>, Vec<F>), ProofVerifyError> {
        let CRR1CSInstance {
            input,
            comm_W,
            aux_W,
        } = instance;

        let input = input.assignment.as_slice();

        // update the transcript
        Transcript::append_scalars(transcript, b"input", input);
        AppendToTranscript::append_to_transcript(comm_W, b"witness", transcript);

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
        let (claim_post_phase1, rx) = self
            .sc_proof_phase1
            .verify::<E>(claim_phase1, num_rounds_x, 3, transcript)?;

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

        // verify the joint claim with a sum-check protocol
        let (claim_post_phase2, ry) = self
            .sc_proof_phase2
            .verify::<E>(claim_phase2, num_rounds_y, 2, transcript)?;


        // Compute (1,io)(r_y) so that we can use it to compute Z(r_y)
        let poly_input_eval = {
            // constant term: one
            let mut input_as_sparse_poly_entries = vec![F::ONE];
            // remaining inputs:
            input_as_sparse_poly_entries.extend(
                (0..input.len())
                    .map(|i| input[i])
                    .collect::<Vec<F>>(),
            );
            SparsePoly::new(n.log_2(), input_as_sparse_poly_entries).evaluate(&ry[1..])
        };

        // compute Z(r_y): eval_Z_at_ry = (F::one() - ry[0]) * self.eval_vars_at_ry + ry[0] * poly_input_eval
        let eval_Z_at_ry = (F::one() - ry[0]) * self.eval_vars_at_ry + ry[0] * poly_input_eval;

        // perform the final check in the second sum-check protocol
        let (eval_A_r, eval_B_r, eval_C_r) = evals;
        let expected_claim_post_phase2 = eval_Z_at_ry * (r_A * eval_A_r + r_B * eval_B_r + r_C * eval_C_r);

        assert_eq!(expected_claim_post_phase2, claim_post_phase2);

        Ok((rx, ry))
    }

    pub fn decide(&self,
                  comm_W: &PC::Commitment,
                  srs: &PC::SRS,
                  ry: Vec<F>) -> Result<(), ProofVerifyError> {
        // verify Z(ry) proof against the initial commitment `comm_W`
        PC::verify(
            srs,
            &ry[1..],
            &self.eval_vars_at_ry,
            comm_W,
            &self.proof_eval_vars_at_ry,
        );

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::nexus_spartan::{crr1cs::produce_synthetic_crr1cs, r1csinstance::R1CSInstance};

    use super::*;
    use crate::constant_for_curves::{E};
    use ark_bls12_381::Fr;
    use ark_ff::PrimeField;
    use ark_std::test_rng;
    use crate::kzh::kzh2::KZH2;

    fn produce_tiny_r1cs<F: PrimeField + Absorb>() -> (R1CSInstance<F>, Vec<F>, Vec<F>) {
        // three constraints over five variables Z1, Z2, Z3, Z4, and Z5
        // rounded to the nearest power of two
        let num_cons = 128;
        let num_vars = 256;
        let num_inputs = 2;

        // encode the above constraints into three matrices
        let mut A: Vec<(usize, usize, F)> = Vec::new();
        let mut B: Vec<(usize, usize, F)> = Vec::new();
        let mut C: Vec<(usize, usize, F)> = Vec::new();

        let one = F::one();
        // constraint 0 entries
        // (Z1 + Z2) * I0 - Z3 = 0;
        A.push((0, 0, one));
        A.push((0, 1, one));
        B.push((0, num_vars + 1, one));
        C.push((0, 2, one));

        // constraint 1 entries
        // (Z1 + I1) * (Z3) - Z4 = 0
        A.push((1, 0, one));
        A.push((1, num_vars + 2, one));
        B.push((1, 2, one));
        C.push((1, 3, one));
        // constraint 3 entries
        // Z5 * 1 - 0 = 0
        A.push((2, 4, one));
        B.push((2, num_vars, one));

        let inst = R1CSInstance::new(num_cons, num_vars, num_inputs, &A, &B, &C);

        // compute a satisfying assignment
        let mut prng = test_rng();
        let i0 = F::rand(&mut prng);
        let i1 = F::rand(&mut prng);
        let z1 = F::rand(&mut prng);
        let z2 = F::rand(&mut prng);
        let z3 = (z1 + z2) * i0; // constraint 1: (Z1 + Z2) * I0 - Z3 = 0;
        let z4 = (z1 + i1) * z3; // constraint 2: (Z1 + I1) * (Z3) - Z4 = 0
        let z5 = F::zero(); //constraint 3

        let mut vars = vec![F::zero(); num_vars];
        vars[0] = z1;
        vars[1] = z2;
        vars[2] = z3;
        vars[3] = z4;
        vars[4] = z5;

        let mut input = vec![F::zero(); num_inputs];
        input[0] = i0;
        input[1] = i1;

        (inst, vars, input)
    }

    #[test]
    fn test_tiny_r1cs() {
        test_tiny_r1cs_helper::<Fr>()
    }

    fn test_tiny_r1cs_helper<F: PrimeField + Absorb>() {
        let (inst, vars, input) = produce_tiny_r1cs::<F>();
        let is_sat = inst.is_sat(&vars, &input);
        assert!(is_sat);
    }

    #[test]
    fn test_synthetic_r1cs() {
        test_synthetic_r1cs_helper::<Fr>()
    }

    fn test_synthetic_r1cs_helper<F: PrimeField + Absorb>() {
        let (inst, vars, input) = R1CSInstance::<F>::produce_synthetic_r1cs(1024, 1024, 10);
        let is_sat = inst.is_sat(&vars, &input);
        assert!(is_sat);
    }

    #[test]
    pub fn check_crr1cs_proof() {
        check_crr1cs_proof_helper::<E, KZH2<E>>()
    }

    fn check_crr1cs_proof_helper<E: Pairing, PC: KZH<E>>()
    where
        <E as Pairing>::ScalarField: Absorb,
    {
        let num_vars = 1024;
        let num_cons = num_vars;
        let num_inputs = 10;
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

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(
                num_vars,
                num_cons,
                &instance,
                &inst_evals,
                &mut verifier_transcript,
            )
            .is_ok());

        assert!(proof.decide(&instance.comm_W, &gens.gens_r1cs_sat, ry).is_ok());
    }
}

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use rand::RngCore;

use crate::accumulation;
use crate::accumulation::accumulator::{AccInstance, AccWitness, Accumulator};
use crate::constant_for_curves::ScalarField;
use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::pcs::multilinear_pcs::{Commitment, PolyCommit};
use crate::polynomial::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::Transcript;
use ark_ff::Zero;

// XXX move to mod.rs or somewhere neutral
#[derive(Clone, Debug)]
pub struct SRS<E: Pairing> {
    pub acc_srs: accumulation::accumulator::AccSRS<E>,
}

impl<E: Pairing> SRS<E>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    <E as Pairing>::ScalarField: Absorb,
{
    fn new<T: RngCore>(degree_x: usize, degree_y: usize, rng: &mut T) -> Self {
        let pcs_srs = PolyCommit::setup(degree_x, degree_y, rng);

        SRS {
            acc_srs: Accumulator::setup(pcs_srs, rng),
        }
    }
}

#[derive(Clone, Debug)]
pub struct SignatureAggrData<E: Pairing>
where
    <E as Pairing>::ScalarField: Absorb,
{
    // TODO comment this out for now. we will figure out the BLS stuff later.
    //pk: E::G1Affine,
    //sig: E::G2Affine,
    // Commitments to b_1(x) and b_2(x)
    B_1_commitment: Option<Commitment<E>>,
    B_2_commitment: Option<Commitment<E>>,
    // c(x): the union poly
    bitfield_poly: MultilinearPolynomial<E::ScalarField>,
    // Commitment to c(x)
    bitfield_commitment: Commitment<E>,
    sumcheck_proof: Option<SumcheckInstanceProof<E::ScalarField>>,
    // Accumulator witness for random evaluation of p(x) at rho:
    // p(rho) = b_1(rho) + c_1 * b_2(rho) + c_2 * c(rho)
    sumcheck_eval_acc_witness: Option<AccWitness<E>>,
    // Evaluations of the inner polynomials at rho:
    b_1_at_rho: Option<E::ScalarField>, // b_1(rho)
    b_2_at_rho: Option<E::ScalarField>, // b_2(rho)
    c_at_rho: Option<E::ScalarField>, // c(rho)
    // ivc_proof: Option<IVCProof<E>>
    // state_acc_witness: Option<AccWitness<E>>
}

impl<E: Pairing> SignatureAggrData<E>
where
    <E as Pairing>::ScalarField: Absorb,
{
    pub fn new(bitfield_poly: MultilinearPolynomial<E::ScalarField>, _sumcheck_proof: Option<SumcheckInstanceProof<E::ScalarField>>, srs: &SRS<E>) -> Self {
        // XXX this PolyCommit is not very ergonomic
        let poly_commit = PolyCommit { srs: srs.acc_srs.pc_srs.clone() }; // XXX no clone
        let bitfield_commitment = poly_commit.commit(&bitfield_poly);
        SignatureAggrData {
            B_1_commitment: None,
            B_2_commitment: None,
            bitfield_poly,
            bitfield_commitment,
            sumcheck_proof: None,
            sumcheck_eval_acc_witness: None,
            b_1_at_rho: None,
            b_2_at_rho: None,
            c_at_rho: None,
        }
    }
}

/// This struct represents an aggregator on the network that receives data from two parties and needs to aggregate them
/// into one. After aggregation, the aggregator forwards the aggregated data.
pub struct Aggregator<E, F>
where
    E: Pairing<ScalarField=F>,
    F: PrimeField + Absorb,
{
    srs: SRS<E>,
    A_1: SignatureAggrData<E>,
    A_2: SignatureAggrData<E>,
}

impl<E, F> Aggregator<E, F>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    /// Return (A.X, A.W) given f(x), and z and y such that f(z) = y
    fn get_accumulator_from_evaluation(&self,
                                       bitfield_poly: &MultilinearPolynomial<F>,
                                       bitfield_commitment: &Commitment<E>,
                                       eval_result: &F,
                                       eval_point: &Vec<F>,
    ) -> Accumulator<E> {
        let poly_commit = PolyCommit { srs: self.srs.acc_srs.pc_srs.clone() }; // XXX no clone. bad ergonomics

        // Split the evaluation point in half since open() just needs the first half
        // XXX ergonomics
        assert_eq!(eval_point.len() % 2, 0);
        let mid = eval_point.len() / 2;
        let (eval_point_first_half, eval_point_second_half) = eval_point.split_at(mid);

        let opening_proof = poly_commit.open(bitfield_poly, bitfield_commitment.clone(), eval_point_first_half); // XXX needless clone

        // XXX bad name for function
        let acc_instance = Accumulator::new_accumulator_instance_from_proof(
            &self.srs.acc_srs,
            &bitfield_commitment.C,
            eval_point_first_half,
            eval_point_second_half,
            eval_result,
        );
        let acc_witness = Accumulator::new_accumulator_witness_from_proof(
            &self.srs.acc_srs,
            opening_proof,
            eval_point_first_half,
            eval_point_second_half,
        );
        Accumulator {
            witness: acc_witness,
            instance: acc_instance,
        }
    }

    pub fn aggregate(&self, transcript: &mut Transcript<F>) -> SignatureAggrData<E> {
        let poly_commit = PolyCommit { srs: self.srs.acc_srs.pc_srs.clone() }; // XXX no clone. bad ergonomics
        // Step 1:
        // let pk = self.A_1.pk + self.A_2.pk;
        // let sk = self.A_1.sig + self.A_2.sig;

        // Step 2: Compute c(x)
        let b_1_poly = &self.A_1.bitfield_poly;
        let b_2_poly = &self.A_2.bitfield_poly;

        let c_poly = b_1_poly.get_bitfield_union_poly(&b_2_poly);
        let C_commitment = poly_commit.commit(&c_poly);

        // Step 3: Get r from verifier: it's the evaluation point challenge (for the zerocheck)
        transcript.append_point::<E>(b"poly", &C_commitment.C);
        let vec_r = transcript.challenge_vector(b"vec_r", b_1_poly.num_variables);

        // Step 4: Do the sumcheck for the following polynomial:
        // eq(r,x) * (b_1 + b_2 - b_1 * b_2 - c)
        let union_comb_func =
            |eq_poly: &F, b_1_poly: &F, b_2_poly: &F, c_poly: &F|
             -> F { *eq_poly * (*b_1_poly + *b_2_poly - *b_1_poly * *b_2_poly - *c_poly) };

        // Start preparing for the sumcheck
        let num_rounds = c_poly.num_variables;
        let eq_at_r = MultilinearPolynomial::new(EqPolynomial::new(vec_r).evals());

        // Sanity check: This is not true in general, but it's true for our tests
        assert_eq!(b_1_poly.len, b_2_poly.len);
        assert_eq!(b_1_poly.len, c_poly.len);
        assert_eq!(b_1_poly.len, eq_at_r.len);

        // Run the sumcheck and get back the verifier's challenge (random eval point rho)
        let (sumcheck_proof, sumcheck_challenges, _) =
            SumcheckInstanceProof::prove_cubic_four_terms::<_, E::G1>(&F::zero(),
                                                                      num_rounds,
                                                                      &mut eq_at_r.clone(), // eq(r, x)
                                                                      &mut b_1_poly.clone(), // b_1(x)
                                                                      &mut b_2_poly.clone(), // b_2(x)
                                                                      &mut c_poly.clone(), // c(x)
                                                                      union_comb_func,
                                                                      transcript);

        // Step 5: Send evaluations to verifier
        // The verifier needs the following evaluation to verify the sumcheck:
        // y_1 = b_1(rho), y_2 = b_2(rho), and y_3 = c(rho)
        // where rho are the sumcheck challenges.
        //
        // Instead of sending three KZH proofs to the verifier, we ask the verifier for challenges c_1 and c_2
        // then we combine three polys into a single polynomial using a random linear combination, and send a
        // proof for the resulting polynomial p(x) where p(x) = b_1(x) + c_1 * b_2(x) + c_2 * c(x)

        // Get c_1 and c_2 (XXX could also get just c and then compute c^2)
        let vec_c: Vec<E::ScalarField> = transcript.challenge_vector(b"vec_c", 2);

        // Step 5.1: First compute p(x):
        // Get c_1 * b_2(x)
        let mut c_1_times_b_2_poly = b_2_poly.clone();
        c_1_times_b_2_poly.scalar_mul(&vec_c[0]);

        // Get c_2 * c(x)
        let mut c_2_times_c_poly = c_poly.clone();
        c_2_times_c_poly.scalar_mul(&vec_c[1]);

        // Now combine everything to p(x)
        let p_x = b_1_poly.clone() + c_1_times_b_2_poly + c_2_times_c_poly;

        // Step 5.2: Now compute P commitment to p(x)
        // First compute c_1 * B_2
        // let b_1_eval_accumulator = self.get_accumulator_from_evaluation(
        //     &self.A_1.bitfield_poly,
        //     &self.A_1.bitfield_commitment,
        //     &sumcheck_challenges,
        // );
        let mut c_1_times_B_2 = self.A_2.bitfield_commitment.clone();
        c_1_times_B_2.scale_by_r(&vec_c[0]);
        // Now compute c_2 * C
        let mut c_2_times_C = C_commitment.clone();
        c_2_times_C.scale_by_r(&vec_c[1]);
        let P_commitment = self.A_1.bitfield_commitment.clone() + c_1_times_B_2 + c_2_times_C;

        // Step 5.3: Compute b_1(rho), b_2(rho), c(rho) to send it to verifier
        let b_1_at_rho = b_1_poly.evaluate(&sumcheck_challenges);
        let b_2_at_rho = b_2_poly.evaluate(&sumcheck_challenges);
        let c_at_rho = c_poly.evaluate(&sumcheck_challenges);
        let p_at_rho = b_1_at_rho + vec_c[0] * b_2_at_rho + vec_c[1] * c_at_rho;

        // Step 5.4: Compute accumulator for opening of p(rho)
        let sumcheck_eval_accumulator = self.get_accumulator_from_evaluation(
            &p_x,
            &P_commitment,
            &p_at_rho,
            &sumcheck_challenges,
        );

        // Step 6: Aggregate accumulators 5-to-1:
        // Hossein: At this point we will also have four more accumulators from the IVC proofs of Bob and Charlie
        // Hossein: Accumulate the five accumulators into one
        // let bob_accumulator_ivc = self.A_1.ivc_proof.acc_witness;
        // let bob_accumulator_state = self.A_1.state_acc_witness;
        // let charlie_accumulator_ivc = self.A_2.ivc_proof.acc_witness;
        // let charlie_accumulator_state = self.A_2.state_acc_witness;
        // let (ivc_proof, state_accumulator) = self.accumulate_everything(sumcheck_eval_accumulator, bob_accumulator_ivc, bob_accumulator_state,
        //                                                                 charlie_accumulator_ivc, charlie_accumulator_state);

        // Step 7: Accumulate forward the evaluations of A, B, C from the IVC proofs
        // let A_eval_1 = self.A_1.ivc_proof.A_eval;
        // let A_eval_2 = self.A_2.ivc_proof.A_eval;
        // let A_eval_3 = ivc_proof.A_eval;

        SignatureAggrData {
            B_1_commitment: Some(self.A_1.bitfield_commitment.clone()),
            B_2_commitment: Some(self.A_2.bitfield_commitment.clone()),
            bitfield_poly: c_poly,
            bitfield_commitment: C_commitment,
            sumcheck_proof: Some(sumcheck_proof),
            sumcheck_eval_acc_witness: Some(sumcheck_eval_accumulator.witness),
            b_1_at_rho: Some(b_1_at_rho),
            b_2_at_rho: Some(b_2_at_rho),
            c_at_rho: Some(c_at_rho),
            // ivc_proof: ivc_proof
            // state_acc_witness: state_acc_witness
        }
    }
}

/// This struct represents a network node that just received an aggregate signature. The verifier needs to verify the
/// aggregate signature (and later aggregate it with more signatures herself).
/// For the purposes of this module, we will only do the verification.
pub struct Verifier<E, F>
where
    F: PrimeField + Absorb,
    E: Pairing<ScalarField=F>,
{
    pub srs: SRS<E>,
    pub A: SignatureAggrData<E>,
}

impl<E, F> Verifier<E, F>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
    E: Pairing<ScalarField=F>,
{
    fn get_acc_instance_from_evaluation(&self,
                                        bitfield_commitment: &Commitment<E>,
                                        eval_result: &F,
                                        eval_point: &Vec<F>,
    ) -> AccInstance<E> {
        // Split the evaluation point in half since open() just needs the first half
        // XXX ergonomics
        assert_eq!(eval_point.len() % 2, 0);
        let mid = eval_point.len() / 2;
        let (eval_point_first_half, eval_point_second_half) = eval_point.split_at(mid);

        // XXX bad name for function
        Accumulator::new_accumulator_instance_from_proof(
            &self.srs.acc_srs,
            &bitfield_commitment.C,
            eval_point_first_half,
            eval_point_second_half,
            eval_result,
        )
    }

    pub fn verify(&self, transcript: &mut Transcript<F>) -> (bool, Vec<F>) {
        // Step 1: Get r challenge from verifier
        transcript.append_point::<E>(b"poly", &self.A.bitfield_commitment.C);
        let vec_r = transcript.challenge_vector(b"vec_r", self.A.bitfield_poly.num_variables);

        // Step 2: Verify the sumcheck proof
        let zero = F::zero();
        let num_rounds = self.A.bitfield_poly.num_variables;
        let (tensorcheck_claim, sumcheck_challenges) =
            self.A.sumcheck_proof.clone()
                .unwrap()
                .verify::<E>(
                    zero,
                    num_rounds,
                    3,
                    transcript,
                ).unwrap();

        // Step 3: Verify the sumcheck tensor check (the random evaluation at the end of the protocol)
        // We need to check: p(rho) = tensor check_claim
        // where rho are the sumcheck challenges and
        // where p(x) = eq(r,x) (b_1(x) + b_2(x) - b_1(x) * b_2(x) - c(x))
        let eq_at_r = MultilinearPolynomial::new(EqPolynomial::new(vec_r).evals());
        let eq_at_r_rho = eq_at_r.evaluate(&sumcheck_challenges);
        let b_1_at_rho = self.A.b_1_at_rho.unwrap();
        let b_2_at_rho = self.A.b_2_at_rho.unwrap();
        let c_at_rho = self.A.c_at_rho.unwrap();
        assert_eq!(tensorcheck_claim, eq_at_r_rho * (b_1_at_rho + b_2_at_rho - b_1_at_rho * b_2_at_rho - c_at_rho));

        // Step 4: Verify the IVC proof

        // Step 5: Verify the BLS signature

        (true, sumcheck_challenges)
    }

    pub fn decide(&self, transcript: &mut Transcript<E::ScalarField>, sumcheck_challenges: Vec<E::ScalarField>) -> bool {
        let b_1_at_rho = self.A.b_1_at_rho.unwrap();
        let b_2_at_rho = self.A.b_2_at_rho.unwrap();
        let c_at_rho = self.A.c_at_rho.unwrap();

        // Verify the accumulator
        // Get c_1 and c_2 (XXX could also get just c and then compute c^2)
        let vec_c: Vec<F> = transcript.challenge_vector(b"vec_c", 2);

        // Now compute commitment to P using B_1, B_2, and C
        let mut c_1_times_B_2 = self.A.B_2_commitment.clone().unwrap(); // XXX stop the cloning!!!
        c_1_times_B_2.scale_by_r(&vec_c[0]);
        let mut c_2_times_C = self.A.bitfield_commitment.clone();
        c_2_times_C.scale_by_r(&vec_c[1]);
        let P_commitment = self.A.B_1_commitment.clone().unwrap() + c_1_times_B_2 + c_2_times_C;

        // Now compute p(rho)
        let p_at_rho = b_1_at_rho + vec_c[0] * b_2_at_rho + vec_c[1] * c_at_rho;

        // Compute the decider's accumulator instance
        let acc_instance = self.get_acc_instance_from_evaluation(
            &P_commitment,
            &p_at_rho,
            &sumcheck_challenges);

        // Compute the full accumulator using the witness from the prover
        let accumulator = Accumulator {
            instance: acc_instance,
            witness: self.A.sumcheck_eval_acc_witness.clone().unwrap(),
        };

        // Decide the accumulator!
        assert!(Accumulator::decide(&self.srs.acc_srs, &accumulator));

        true
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use crate::constant_for_curves::{ScalarField, E};

    type F = ScalarField;

    #[test]
    fn test_aggregate() {
        let rng = &mut rand::thread_rng();
        let mut transcript_p = Transcript::<F>::new(b"aggr");
        let mut transcript_v = Transcript::<F>::new(b"aggr");

        // num_vars = log(degree_x) + log(degree_y)
        let degree_x = 8usize;
        let degree_y = 8usize;
        let num_vars = 6usize;

        let srs = SRS::<E>::new(degree_x, degree_y, rng);

        let b_1 = MultilinearPolynomial::random_binary(num_vars, rng);
        let sig_aggr_data_1 = SignatureAggrData::new(b_1, None, &srs);

        let b_2 = MultilinearPolynomial::random_binary(num_vars, rng);
        let sig_aggr_data_2 = SignatureAggrData::new(b_2, None, &srs);

        let aggregator = Aggregator {
            srs: srs.clone(),
            A_1: sig_aggr_data_1,
            A_2: sig_aggr_data_2,
        };

        let agg_data = aggregator.aggregate(&mut transcript_p);

        // Now let's do verification
        let verifier = Verifier {
            srs,
            A: agg_data,
        };

        let (is_valid, sumcheck_challenges) = verifier.verify(&mut transcript_v);
        assert!(is_valid);

        assert!(verifier.decide(&mut transcript_v, sumcheck_challenges));
    }
}


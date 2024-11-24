use crate::kzh_fold::kzh2_fold::{Acc2Instance, Acc2SRS, Accumulator2 as KZHAccumulator, Accumulator2};
use crate::kzh::kzh2::{split_between_x_and_y, KZH2Commitment, KZH2};
use crate::kzh::KZH;
use crate::math::Math;
use crate::nexus_spartan::sumcheck::SumcheckInstanceProof;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_ff::UniformRand;
use rand::RngCore;

#[derive(Clone, Debug)]
pub struct SignatureAggrSRS<E: Pairing> {
    pub acc_srs: Acc2SRS<E>,
}


/// This is the data sent by an aggregator node to the next aggregator node on the network
#[derive(Clone, Debug)]
pub struct SignatureAggrData<E, F>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    /////////////// Signature aggregation data ///////////////

    /// Commitments to b_1(x) and b_2(x)
    pub B_1_commitment: KZH2Commitment<E>,
    pub B_2_commitment: KZH2Commitment<E>,
    /// c(x): the union poly
    pub bitfield_poly: MultilinearPolynomial<F>,

    /// aggregated signature and message
    sig: E::G2Affine,

    /////////////// All the z_i data that go into SignatureVerifierCircuit //////////////

    /// aggregated public key
    pub pk: E::G1Affine,

    /// Commitment to c(x)
    pub bitfield_commitment: KZH2Commitment<E>,

    /// proof that bitfield poly has been computed correctly in fact
    pub sumcheck_proof: SumcheckInstanceProof<F>,

    /// Evaluations of the inner polynomials at rho:
    pub b_1_at_rho: F, // b_1(rho)
    pub b_2_at_rho: F, // b_2(rho)
    pub c_at_rho: F, // c(rho)

    /////////////// KZH accumulator for the sig aggr sumcheck (goes into 3-to-1) //////////////

    /// Accumulator for random evaluation of p(x) at rho:
    /// p(rho) = b_1(rho) + c_1 * b_2(rho) + c_2 * c(rho)
    pub sumcheck_eval_KZH_accumulator: KZHAccumulator<E>,

    /////////////// IVC proof for all the previous steps `\pi_i` //////////////////

    // The IVC proof contains a KZH accumulator and an A,B,C accumulator

    // ivc_proof: Option<CRR1CSProof<E>>
}

impl<E, F> SignatureAggrData<E, F>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    /// Generate a random SignatureAggrData for a network participant
    pub fn rand<R: RngCore>(num_vars: usize, srs: &Acc2SRS<E>, rng: &mut R) -> Self {
        let mut transcript = Transcript::<F>::new(b"aggr");

        // Random polynomials
        let b_1_poly = MultilinearPolynomial::random_binary(num_vars, rng);
        let b_2_poly = MultilinearPolynomial::random_binary(num_vars, rng);
        let c_poly = b_1_poly.get_bitfield_union_poly(&b_2_poly);

        let B_1_commitment = KZH2::commit(&srs.pc_srs, &b_1_poly);
        let B_2_commitment = KZH2::commit(&srs.pc_srs, &b_2_poly);
        let C_commitment = KZH2::commit(&srs.pc_srs, &c_poly);

        // Some random stuff for the sig/pk
        let sig = E::G2Affine::rand(rng);
        let pk = E::G1Affine::rand(rng);
        // let message = E::G2Affine::rand(rng);

        // Perform the sig aggr sumcheck
        let (sumcheck_proof, rho) = perform_sig_aggr_sumcheck::<E, F>(&b_1_poly, &b_2_poly, &c_poly,
                                                                      &mut transcript);

        // Get the accumulator for the sumcheck
        let (b_1_at_rho, b_2_at_rho, c_at_rho, sumcheck_eval_KZH_accumulator) =
            compute_signature_aggr_KZH_accumulator(srs, &b_1_poly, &b_2_poly, &c_poly, &rho, &mut transcript);

        Self {
            B_1_commitment,
            B_2_commitment,
            bitfield_poly: c_poly,
            sig,
            pk,
            bitfield_commitment: C_commitment,
            sumcheck_proof,
            b_1_at_rho,
            b_2_at_rho,
            c_at_rho,
            sumcheck_eval_KZH_accumulator,
        }
    }
}

/// This struct represents an IVC aggregator, Alice, that receives network data from a single party, Bob, and also has
/// her own running accumulator. It aggregates the received data with the running accumulator and produces her own
/// `SignatureAggrData` that can be forwarded to the next node.
pub struct AggregatorIVC<E, F>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    pub srs: SignatureAggrSRS<E>,

    // Alice's running bitfield (the result of previous aggregations)
    pub running_bitfield_poly: MultilinearPolynomial<E::ScalarField>,
    // Commitment to Alice's running bitfield
    pub running_bitfield_commitment: KZH2Commitment<E>,
    // Alice's running accumulator
    pub running_accumulator: Accumulator2<E>,
    // running signature
    pub running_signature: E::G2Affine,
    // running public key
    pub running_public_key: E::G1Affine,
    // the message, this is supposed to be constant during the IVC/PCD
    // pub message: E::G2Affine,

    // Data received from Bob
    pub bob_data: SignatureAggrData<E, F>,
}

// Perform sumcheck for the following polynomial: eq(r,x) * (b_1 + b_2 - b_1 * b_2 - c)
pub fn perform_sig_aggr_sumcheck<E, F>(
    b_1_poly: &MultilinearPolynomial<F>,
    b_2_poly: &MultilinearPolynomial<F>,
    c_poly: &MultilinearPolynomial<F>,
    transcript: &mut Transcript<F>,
) -> (SumcheckInstanceProof<F>, Vec<F>)
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    let vec_r = transcript.challenge_vector(b"vec_r", b_1_poly.num_variables);
    let union_comb_func =
        |eq_poly: &F, b_1_poly: &F, b_2_poly: &F, c_poly: &F|
         -> F { *eq_poly * (*b_1_poly + *b_2_poly - *b_1_poly * *b_2_poly - *c_poly) };

    // Start preparing for the sumcheck
    let num_rounds = c_poly.num_variables;
    let eq_at_r = MultilinearPolynomial::new(EqPolynomial::new(vec_r).evals());

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
    let rho = sumcheck_challenges;

    (sumcheck_proof, rho)
}

/// Return (A.X, A.W) given f(x), and z and y such that f(z) = y
fn get_accumulator_from_evaluation<E, F>(acc_srs: &Acc2SRS<E>,
                                         bitfield_poly: &MultilinearPolynomial<F>,
                                         eval_result: &F,
                                         eval_point: &Vec<F>,
) -> KZHAccumulator<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    let bitfield_commitment = KZH2::commit(
        &acc_srs.pc_srs,
        bitfield_poly,
    );

    let opening_proof = KZH2::open(
        &acc_srs.pc_srs,
        eval_point,
        &bitfield_commitment,
        &bitfield_poly,
    );


    let length_x = acc_srs.pc_srs.degree_x.log_2();
    let length_y = acc_srs.pc_srs.degree_y.log_2();
    let (eval_point_first_half, eval_point_second_half) = split_between_x_and_y::<F>(length_x, length_y, eval_point, F::ZERO);

    let acc_instance = KZHAccumulator::new_accumulator_instance_from_fresh_kzh_instance(
        acc_srs,
        &bitfield_commitment.C,
        eval_point_first_half.as_slice(),
        eval_point_second_half.as_slice(),
        eval_result,
    );

    let acc_witness = KZHAccumulator::new_accumulator_witness_from_fresh_kzh_witness(
        acc_srs,
        opening_proof,
        eval_point_first_half.as_slice(),
        eval_point_second_half.as_slice(),
    );

    KZHAccumulator {
        witness: acc_witness,
        instance: acc_instance,
    }
}

/// Signature aggregation verifier needs the following evaluation to verify the sumcheck:
/// y_1 = b_1(rho), y_2 = b_2(rho), and y_3 = c(rho)
/// where rho are the sumcheck challenges.
///
/// Instead of sending three KZH proofs to the verifier, we ask the verifier for challenges c_1 and c_2
/// then we combine three polys into a single polynomial using a random linear combination, and send a
/// proof for the resulting polynomial p(x) where p(x) = b_1(x) + c_1 * b_2(x) + c_2 * c(x)
///
/// Return b_1(rho), b_2(rho), c(rho) and an accumulator for the proof on p(x).
pub fn compute_signature_aggr_KZH_accumulator<E, F>(
    acc_srs: &Acc2SRS<E>,
    b_1_poly: &MultilinearPolynomial<F>,
    b_2_poly: &MultilinearPolynomial<F>,
    c_poly: &MultilinearPolynomial<F>,
    rho: &Vec<F>,
    transcript: &mut Transcript<F>,
) -> (F, F, F, KZHAccumulator<E>)
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    // Get c_1 and c_2 (one could also get just c and then compute c^2)
    let vec_c: Vec<F> = transcript.challenge_vector(b"vec_c", 2);

    // Step 5.1: First compute p(x):
    // Get c_1 * b_2(x)
    let mut c_1_times_b_2_poly = b_2_poly.clone();
    c_1_times_b_2_poly.scalar_mul(&vec_c[0]);

    // Get c_2 * c(x)
    let mut c_2_times_c_poly = c_poly.clone();
    c_2_times_c_poly.scalar_mul(&vec_c[1]);

    // Now combine everything to p(x)
    let p_x = b_1_poly.clone() + c_1_times_b_2_poly + c_2_times_c_poly;

    // Step 5.2: Compute b_1(rho), b_2(rho), c(rho) to send it to verifier
    let b_1_at_rho = b_1_poly.evaluate(&rho);
    let b_2_at_rho = b_2_poly.evaluate(&rho);
    let c_at_rho = c_poly.evaluate(&rho);
    let p_at_rho = b_1_at_rho + vec_c[0] * b_2_at_rho + vec_c[1] * c_at_rho;

    // Step 5.4: Compute accumulator for opening of p(rho)
    let sumcheck_eval_KZH_accumulator = get_accumulator_from_evaluation(
        &acc_srs,
        &p_x,
        &p_at_rho,
        &rho,
    );

    (b_1_at_rho, b_2_at_rho, c_at_rho, sumcheck_eval_KZH_accumulator)
}

impl<E, F> AggregatorIVC<E, F>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
{
    pub fn aggregate(&self, transcript: &mut Transcript<F>) -> SignatureAggrData<E, F> {

        // Step 1:
        let pk = self.running_public_key + self.bob_data.pk;
        let sig = self.running_signature + self.bob_data.sig;

        // assert_eq!(self.message, self.bob_data.message, "two messages should be equal");

        // Step 2: Compute c(x)
        let b_1_poly = &self.running_bitfield_poly;
        let b_2_poly = &self.bob_data.bitfield_poly;

        let c_poly = b_1_poly.get_bitfield_union_poly(&b_2_poly);
        let C_commitment = KZH2::commit(&self.srs.acc_srs.pc_srs, &c_poly);

        // Step 3: Get r from verifier: it's the evaluation point challenge (for the zerocheck)
        transcript.append_point::<E>(
            b"poly",
            &C_commitment.C,
        );

        // Step 4: Do the sumcheck for the following polynomial:
        // eq(r,x) * (b_1 + b_2 - b_1 * b_2 - c)
        let (sumcheck_proof, rho) = perform_sig_aggr_sumcheck::<E, F>(b_1_poly, b_2_poly, &c_poly,
                                                                      transcript);

        // Step 5: Send KZH accumulator to verifier
        let (b_1_at_rho, b_2_at_rho, c_at_rho, sumcheck_eval_KZH_accumulator) =
            compute_signature_aggr_KZH_accumulator(&self.srs.acc_srs,
                                                   &b_1_poly, &b_2_poly, &c_poly,
                                                   &rho, transcript);


        // Step 6: Aggregate accumulators 3-to-1:
        // At this point we will also have two more KZH accumulators: one from our running accumulator, and another one from Bob
        // Accumulate thet three accumulators into one
        // let bob_KZH_accumulator = self.bob_data.ivc_proof.KZH_accumulator;
        // let running_KZH_accumulator = self.running_accumulator.KZH_accumulator

        // We will also have two A,B,C accumulators: one from our running accumulator, and one from Bob
        // let bob_A_B_C_eval_accumulator = self.bob_data.ivc_proof.A_B_C_eval_accumulator;
        // let running_A_B_C_eval_accumulator=  self.running_accumulator.A_B_C_eval_accumulator;

        // Accumulate everything!
        // let (ivc_proof, state_accumulator) = self.accumulate_everything(sumcheck_eval_accumulator, bob_accumulator, running_accumulator,
        //                                                                 bob_A_B_C_eval_accumulator, running_A_B_C_eval_accumulator);

        SignatureAggrData {
            B_1_commitment: self.running_bitfield_commitment.clone(),
            B_2_commitment: self.bob_data.bitfield_commitment.clone(),
            bitfield_poly: c_poly,
            sig: sig.into(),
            // message: self.message,
            pk: pk.into(),
            bitfield_commitment: C_commitment,
            sumcheck_proof,
            b_1_at_rho,
            b_2_at_rho,
            c_at_rho,
            sumcheck_eval_KZH_accumulator,
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
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    E: Pairing<ScalarField=F>,
{
    pub srs: SignatureAggrSRS<E>,
    pub A: SignatureAggrData<E, F>,
}

impl<E, F> Verifier<E, F>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    F: PrimeField + Absorb,
    E: Pairing<ScalarField=F>,
{
    fn get_acc_instance_from_evaluation(&self,
                                        bitfield_commitment: &KZH2Commitment<E>,
                                        eval_result: &F,
                                        eval_point: &Vec<F>,
    ) -> Acc2Instance<E> {
        // Split the evaluation point in half since open() just needs the first half
        let length_x = self.srs.acc_srs.pc_srs.degree_x.log_2();
        let length_y = self.srs.acc_srs.pc_srs.degree_y.log_2();
        let (eval_point_first_half, eval_point_second_half) = split_between_x_and_y::<F>(length_x, length_y, eval_point, F::ZERO);
        KZHAccumulator::new_accumulator_instance_from_fresh_kzh_instance(
            &self.srs.acc_srs,
            &bitfield_commitment.C,
            eval_point_first_half.as_slice(),
            eval_point_second_half.as_slice(),
            eval_result,
        )
    }

    pub fn verify(&self, transcript: &mut Transcript<F>) -> (bool, Vec<F>, KZH2Commitment<E>) {
        // Step 1: Get r challenge from verifier
        transcript.append_point::<E>(
            b"poly",
            &self.A.bitfield_commitment.C,
        );
        let vec_r = transcript.challenge_vector(b"vec_r", self.A.bitfield_poly.num_variables);

        // Step 2: Verify the sumcheck proof
        let zero = F::zero();
        let num_rounds = self.A.bitfield_poly.num_variables;
        let (tensor_check_claim, sumcheck_challenges) =
            self.A.sumcheck_proof.clone()
                .verify::<E>(
                    zero,
                    num_rounds,
                    3,
                    transcript,
                ).unwrap();
        let rho = sumcheck_challenges;

        // Step 3: Verify the sumcheck tensor check (the random evaluation at the end of the protocol)
        // We need to check: p(rho) = tensor check_claim
        // where rho are the sumcheck challenges and
        // where p(x) = eq(r,x) (b_1(x) + b_2(x) - b_1(x) * b_2(x) - c(x))
        let eq_at_r = MultilinearPolynomial::new(EqPolynomial::new(vec_r).evals());
        let eq_at_r_rho = eq_at_r.evaluate(&rho);
        let b_1_at_rho = self.A.b_1_at_rho;
        let b_2_at_rho = self.A.b_2_at_rho;
        let c_at_rho = self.A.c_at_rho;
        assert_eq!(tensor_check_claim, eq_at_r_rho * (b_1_at_rho + b_2_at_rho - b_1_at_rho * b_2_at_rho - c_at_rho));

        // Step 4: Compute aggregated commitment P to check against accumulator
        let vec_c: Vec<F> = transcript.challenge_vector(b"vec_c", 2);
        let mut c_1_times_B_2 = self.A.B_2_commitment.clone();
        c_1_times_B_2.scale_by_r(&vec_c[0]);
        let mut c_2_times_C = self.A.bitfield_commitment.clone();
        c_2_times_C.scale_by_r(&vec_c[1]);
        let P_commitment = self.A.B_1_commitment.clone() + c_1_times_B_2 + c_2_times_C;

        (true, rho, P_commitment)
    }

    pub fn decide(&self, P_commitment: KZH2Commitment<E>, vec_c: Vec<F>, sumcheck_challenges: Vec<F>) -> bool {
        let rho = sumcheck_challenges;
        let b_1_at_rho = self.A.b_1_at_rho;
        let b_2_at_rho = self.A.b_2_at_rho;
        let c_at_rho = self.A.c_at_rho;

        // Verify the accumulator
        let p_at_rho = b_1_at_rho + vec_c[0] * b_2_at_rho + vec_c[1] * c_at_rho;

        // Compute the decider's accumulator instance
        let _acc_instance = self.get_acc_instance_from_evaluation(
            &P_commitment,
            &p_at_rho,
            &rho);

        // Do the cross-check that the accumulator is the right one using _acc_instance

        // Decide the accumulator!
        assert!(KZHAccumulator::decide(&self.srs.acc_srs, &self.A.sumcheck_eval_KZH_accumulator.clone()));

        // XXX Verify the IVC proof

        true
    }
}


#[cfg(test)]
pub mod test {
    use crate::kzh_fold::kzh2_fold::Accumulator2;
    use crate::constant_for_curves::{G1Affine, G2Affine, ScalarField, E};
    use crate::kzh::kzh2::KZH2;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use crate::signature_aggregation::signature_aggregation::{AggregatorIVC, SignatureAggrData, SignatureAggrSRS};
    use crate::transcript::transcript::Transcript;
    use ark_std::UniformRand;
    use crate::kzh::KZH;

    type F = ScalarField;

    /// Bob sends signature data to Alice. Alice aggregates it and sends it forward.
    #[test]
    fn test_signature_aggregation_IVC_end_to_end() {
        // Setup:
        let rng = &mut rand::thread_rng();
        let mut transcript_p = Transcript::<F>::new(b"aggr");
        // let mut transcript_v = Transcript::<F>::new(b"aggr");

        // num_vars = log(degree_x) + log(degree_y)
        let degree_x = 64usize;
        let degree_y = 64usize;
        let num_vars = 12usize;
        let srs = SignatureAggrSRS::<E>::new(degree_x, degree_y, rng);

        // Generate signature aggregation payload from Bob
        let bob_data = SignatureAggrData::rand(num_vars, &srs.acc_srs, rng);

        // Generate random running data for Alice
        let alice_bitfield = MultilinearPolynomial::random_binary(num_vars, rng);
        let alice_bitfield_commitment = KZH2::commit(&srs.acc_srs.pc_srs, &alice_bitfield);
        let alice_running_accumulator = Accumulator2::rand(&srs.acc_srs, rng);
        let alice_running_sig = G2Affine::rand(rng);
        let alice_running_pk = G1Affine::rand(rng);

        ////////////// Aggregation ////////////////

        let alice = AggregatorIVC {
            srs: srs.clone(),
            running_bitfield_poly: alice_bitfield,
            running_bitfield_commitment: alice_bitfield_commitment,
            running_accumulator: alice_running_accumulator,
            running_signature: alice_running_sig,
            running_public_key: alice_running_pk,
            bob_data,
        };

        let _aggregated_data = alice.aggregate(&mut transcript_p);

        //////////// Verification //////////////////
    }
}

      

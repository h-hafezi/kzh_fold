use crate::accumulation::accumulator::AccInstance;
use crate::accumulation_circuit::affine_to_projective;
use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit_var::SumcheckCircuitVar;
use crate::nova::cycle_fold::coprocessor_constraints::{OvaInstanceVar, RelaxedOvaInstanceVar};
use crate::polynomial::eq_poly::eq_poly_var::EqPolynomialVar;
use crate::polynomial::multilinear_poly::multilinear_poly_var::MultilinearPolynomialVar;
use crate::signature_aggregation::verifier_circuit::verifier_circuit::SignatureVerifierCircuit;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::CurveConfig;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::ToBitsGadget;
use ark_relations::ns;
use ark_relations::r1cs::{Namespace, SynthesisError};
use std::borrow::Borrow;
use std::marker::PhantomData;

pub struct SignatureVerifierCircuitVar<F, G1, G2, C2>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<ScalarField=F> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G2: SWCurveConfig<
        BaseField=G1::ScalarField,
        ScalarField=G1::BaseField
    >,
{
    /// public keys pk_t = pk_1 + pk_2
    pk_1_var: NonNativeAffineVar<G1>,
    pk_2_var: NonNativeAffineVar<G1>,
    pk_t_var: NonNativeAffineVar<G1>,

    /// bitfield commitment
    pub com_bitfield: (NonNativeFieldVar<G1::BaseField, F>, NonNativeFieldVar<G1::BaseField, F>),

    /// beta
    pub beta: NonNativeFieldVar<G1::BaseField, F>,

    /// cross term error
    pub com_pk_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub cycle_fold_fresh_instance: OvaInstanceVar<G2, C2>,
    pub cycle_fold_running_instance: RelaxedOvaInstanceVar<G2, C2>,
    pub cycle_fold_new_running_instance: RelaxedOvaInstanceVar<G2, C2>,

    /// the bitfield polynomial
    bitfield_poly_var: MultilinearPolynomialVar<F>,

    /// the sumcheck proof
    sumcheck_proof_var: SumcheckCircuitVar<F>,

    /// Evaluations of the inner polynomials at rho:
    b_1_at_rho: FpVar<F>,
    b_2_at_rho: FpVar<F>,
    c_at_rho: FpVar<F>,

}

impl<E, F, G1, G2, C2> AllocVar<SignatureVerifierCircuit<E, F, G1, G2, C2>, F> for SignatureVerifierCircuitVar<F, G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<
        BaseField=G2::ScalarField,
        ScalarField=G2::BaseField
    >,
    F: PrimeField + Absorb,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
{
    fn new_variable<T: Borrow<SignatureVerifierCircuit<E, F, G1, G2, C2>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let verifier_circuit = res.as_ref()
            .map(|e| e.borrow())
            .map_err(|err| *err);

        let pk_1_var = NonNativeAffineVar::new_variable(
            ns!(cs, "pk_1_var"),
            || verifier_circuit.map(|e| affine_to_projective(e.pk_1)),
            mode,
        ).unwrap();

        let pk_2_var = NonNativeAffineVar::new_variable(
            ns!(cs, "pk_2_var"),
            || verifier_circuit.map(|e| affine_to_projective(e.pk_2)),
            mode,
        ).unwrap();

        let pk_t_var = NonNativeAffineVar::new_variable(
            ns!(cs, "pk_t_var"),
            || verifier_circuit.map(|e| affine_to_projective(e.pk_t)),
            mode,
        ).unwrap();

        let com_pk_var = ProjectiveVar::new_variable(
            ns!(cs, "cycle fold running instance"),
            || verifier_circuit.map(|e| e.com_pk.clone()),
            mode,
        ).unwrap();

        let bitfield_poly_var = MultilinearPolynomialVar::new_variable(
            ns!(cs, "bitfield poly var"),
            || verifier_circuit.map(|e| e.bitfield_poly.clone()),
            mode,
        ).unwrap();

        let sumcheck_proof_var = SumcheckCircuitVar::new_variable(
            ns!(cs, "sumcheck proof var"),
            || verifier_circuit.map(|e| e.sumcheck_proof.clone()),
            mode,
        ).unwrap();

        let b_1_at_rho = FpVar::new_variable(
            ns!(cs, "b_1 at rho"),
            || verifier_circuit.map(|e| e.b_1_at_rho.clone()),
            mode,
        ).unwrap();

        let b_2_at_rho = FpVar::new_variable(
            ns!(cs, "b_2 at rho"),
            || verifier_circuit.map(|e| e.b_2_at_rho.clone()),
            mode,
        ).unwrap();

        let c_at_rho = FpVar::new_variable(
            ns!(cs, "c at rho"),
            || verifier_circuit.map(|e| e.c_at_rho.clone()),
            mode,
        ).unwrap();

        // auxiliary inputs
        let cycle_fold_fresh_instance = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary input pk var"),
            || Ok(verifier_circuit.map(|e| e.cycle_fold_fresh_instance.clone()).unwrap()),
            mode,
        ).unwrap();

        // cycle fold instances
        let cycle_fold_running_instance = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "running auxiliary input pk_var"),
            || verifier_circuit.map(|e| e.cycle_fold_running_instance.clone()),
            mode,
        ).unwrap();

        let cycle_fold_new_running_instance = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "final auxiliary input pk var"),
            || verifier_circuit.map(|e| e.cycle_fold_final_instance.clone()),
            mode,
        ).unwrap();

        let com_bitfield_x = NonNativeFieldVar::new_variable(
            ns!(cs, "com bitfield x"),
            || verifier_circuit.map(|e| e.com_bitfield.0.clone()),
            mode,
        ).unwrap();

        let com_bitfield_y = NonNativeFieldVar::new_variable(
            ns!(cs, "com bitfield y"),
            || verifier_circuit.map(|e| e.com_bitfield.1.clone()),
            mode,
        ).unwrap();

        let beta = NonNativeFieldVar::new_variable(
            ns!(cs, "beta"),
            || verifier_circuit.map(|e| e.beta.clone()),
            mode,
        ).unwrap();

        Ok(SignatureVerifierCircuitVar {
            pk_1_var,
            pk_2_var,
            pk_t_var,
            cycle_fold_fresh_instance,
            cycle_fold_running_instance,
            cycle_fold_new_running_instance,
            com_pk_var,
            bitfield_poly_var,
            sumcheck_proof_var,
            b_1_at_rho,
            b_2_at_rho,
            c_at_rho,
            com_bitfield: (com_bitfield_x, com_bitfield_y),
            beta,
        })
    }
}

impl<F, G1, G2, C2> SignatureVerifierCircuitVar<F, G1, G2, C2>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<ScalarField=F> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G2: SWCurveConfig<
        BaseField=G1::ScalarField,
        ScalarField=G1::BaseField
    >,
{
    pub fn verify(&self, transcript: &mut TranscriptVar<F>) {
        // Step 1: Get challenge
        transcript.append_scalar_non_native(b"poly", &self.com_bitfield.0);
        transcript.append_scalar_non_native(b"poly", &self.com_bitfield.1);

        let vec_r = transcript.challenge_vector(b"vec_r", self.bitfield_poly_var.num_variables);

        // Step 2: Verify the sumcheck proof
        let zero: FpVar<F> = FpVar::zero();

        // assert the sumcheck proof is indeed well-formatted
        self.sumcheck_proof_var.claim.enforce_equal(&zero).expect("equality error");
        assert_eq!(self.sumcheck_proof_var.num_rounds, self.bitfield_poly_var.num_variables);
        assert_eq!(self.sumcheck_proof_var.degree_bound, 3);

        let (tensor_check_claim, sumcheck_challenges) = self.sumcheck_proof_var.verify(transcript);

        // Step 3: Verify the sumcheck tensor check (the random evaluation at the end of the protocol)
        // We need to check: p(rho) = tensor check_claim
        // where rho are the sumcheck challenges and
        // where p(x) = eq(r,x) (b_1(x) + b_2(x) - b_1(x) * b_2(x) - c(x))
        let eq_at_r_rho = MultilinearPolynomialVar::new(EqPolynomialVar::new(vec_r).evals()).evaluate(&sumcheck_challenges);
        FpVar::enforce_equal(
            &tensor_check_claim,
            &(eq_at_r_rho * (self.b_1_at_rho.clone() + self.b_2_at_rho.clone() - self.b_1_at_rho.clone() * self.b_2_at_rho.clone() - self.c_at_rho.clone())),
        ).expect("equality error");

        // Step 4: Do the cycle fold math
        // Non-native scalar multiplication: linear combination pk_t = pk_1 + pk_2
        let (flag,
            r,
            pk_1,
            pk_2,
            pk_t
        ) = self.cycle_fold_fresh_instance.parse_secondary_io::<G1>().unwrap();
        // g1 == pk_1
        pk_1.enforce_equal(&self.pk_1_var).expect("error while enforcing equality");
        // g2 == pk_2
        pk_2.enforce_equal(&self.pk_2_var).expect("error while enforcing equality");
        // enforce flag to be true
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        // check r to be equal to one
        r.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        // check out the result pk_t is consistent with input pk_t
        pk_t.enforce_equal(&self.pk_t_var).expect("error while enforcing equality");

        // Step 5: fold the cycle fold instance
        let one_bits = <NonNativeFieldVar<G1::BaseField, F> as ToBitsGadget<F>>::to_bits_le(&self.beta).unwrap();
        let final_instance = self.cycle_fold_running_instance.fold(
            &[((&self.cycle_fold_fresh_instance, None), &self.com_pk_var, &self.beta, &one_bits)]
        ).unwrap();

        self.cycle_fold_new_running_instance.X.enforce_equal(&final_instance.X).expect("XXX: panic message");
        self.cycle_fold_new_running_instance.commitment.enforce_equal(&final_instance.commitment).expect("XXX: panic message");
    }
}

#[cfg(test)]
mod test {
    use crate::commitment::CommitmentScheme;
    use crate::constant_for_curves::{BaseField, G1Affine, ScalarField, E, G1, G2};
    use crate::hash::pederson::PedersenCommitment;
    use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
    use crate::nova::cycle_fold::coprocessor::setup_shape;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use crate::signature_aggregation::signature_aggregation::{Aggregator, SignatureAggrData, Verifier, SRS};
    use crate::signature_aggregation::verifier_circuit::prover::SignatureVerifierProver;
    use crate::signature_aggregation::verifier_circuit::verifier_circuit::SignatureVerifierCircuit;
    use crate::signature_aggregation::verifier_circuit::verifier_circuit_var::SignatureVerifierCircuitVar;
    use crate::transcript::transcript::Transcript;
    use crate::transcript::transcript_var::TranscriptVar;
    use ark_ec::short_weierstrass::{Affine, Projective};
    use ark_ec::AffineRepr;
    use ark_ff::AdditiveGroup;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_relations::r1cs::ConstraintSystem;
    use rand::thread_rng;

    type GrumpkinCurveGroup = ark_grumpkin::Projective;
    type C2 = PedersenCommitment<GrumpkinCurveGroup>;
    type F = ScalarField;
    type Q = BaseField;

    #[test]
    pub fn test() {
        // get a random prover
        let rng = &mut thread_rng();

        // the randomness used to take linear combination for cycle fold
        let beta = Q::from(2u8);

        // get a random prover
        let prover = SignatureVerifierProver::<G1, G2, C2, E>::rand(rng, beta);

        // get a random verifier
        let verifier = {
            let mut transcript_p = Transcript::<F>::new(b"aggr");

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

            verifier
        };

        // get a set of random public keys
        let (pk_1, pk_2, _pk_t): (G1Affine, G1Affine, G1Affine) = SignatureVerifierProver::<G1, G2, C2, E>::get_satisfying_public_keys(rng);

        let shape = setup_shape::<G1, G2>().unwrap();
        let commitment_pp: Vec<Affine<G2>> = PedersenCommitment::<Projective<G2>>::setup(shape.num_vars + shape.num_constraints, b"test", &());

        let (com_pk, cycle_fold_fresh_instance, cycle_fold_running_instance, cycle_fold_final_instance) = {
            // fold it with the prover's running instance/witness
            let (instance, witness) = SignatureVerifierProver::<G1, G2, C2, E>::get_auxiliary_input_for_public_keys(&shape, &commitment_pp, pk_1, pk_2);
            let (com_T, new_running_instance, _) = prover.compute_cycle_fold_proofs_and_final_instance(&instance, &witness);

            (com_T, instance, prover.cycle_fold_running_instance, new_running_instance)
        };

        let signature_verifier_circuit = SignatureVerifierCircuit::<E, F, G1, G2, C2> {
            pk_1,
            pk_2,
            pk_t: (pk_1 + pk_2).into(),
            com_bitfield: (verifier.A.bitfield_commitment.C.x().unwrap(),
                           verifier.A.bitfield_commitment.C.y().unwrap()
            ),
            beta,
            com_pk,
            cycle_fold_fresh_instance,
            cycle_fold_running_instance,
            cycle_fold_final_instance,
            bitfield_poly: verifier.A.bitfield_poly.clone(),
            sumcheck_proof: SumcheckCircuit {
                compressed_polys: verifier.A.sumcheck_proof.unwrap().compressed_polys,
                claim: F::ZERO,
                num_rounds: verifier.A.bitfield_poly.num_variables,
                degree_bound: 3,
            },
            b_1_at_rho: verifier.A.b_1_at_rho.unwrap(),
            b_2_at_rho: verifier.A.b_2_at_rho.unwrap(),
            c_at_rho: verifier.A.c_at_rho.unwrap(),
        };

        let cs = ConstraintSystem::<F>::new_ref();
        let signature_verifier_circuit_var = SignatureVerifierCircuitVar::new_variable(
            cs.clone(),
            || Ok(signature_verifier_circuit),
            AllocationMode::Input,
        ).unwrap();

        let mut transcript_var = TranscriptVar::<F>::new(cs.clone(), b"aggr");
        signature_verifier_circuit_var.verify(&mut transcript_var);

        assert!(cs.is_satisfied().unwrap());
        println!("{} {}", cs.num_constraints(), cs.borrow().unwrap().num_instance_variables);
    }
}

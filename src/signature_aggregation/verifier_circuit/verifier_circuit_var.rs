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

    /// cross term error
    pub com_pk_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub auxiliary_input_pk_var: OvaInstanceVar<G2, C2>,
    pub running_auxiliary_input_pk_var: RelaxedOvaInstanceVar<G2, C2>,
    pub final_auxiliary_input_pk_var: RelaxedOvaInstanceVar<G2, C2>,

    /// the bitfield polynomial
    bitfield_poly_var: MultilinearPolynomialVar<F>,
    /// commitment to bitfield (only for Fiat-Shamir)
    pub bitfield_poly_commitment: NonNativeAffineVar<G1>,

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

        let bitfield_poly_commitment = NonNativeAffineVar::new_variable(
            ns!(cs, "bitfield poly commitment"),
            || verifier_circuit.map(|e| affine_to_projective(e.bitfield_poly_commitment)),
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
        let auxiliary_input_pk_var = OvaInstanceVar::new_variable(
            ns!(cs, "auxiliary input pk var"),
            || Ok(verifier_circuit.map(|e| e.auxiliary_input_pk.clone()).unwrap()),
            mode,
        ).unwrap();

        // cycle fold instances
        let running_auxiliary_input_pk_var = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "running auxiliary input pk_var"),
            || verifier_circuit.map(|e| e.running_auxiliary_input_pk.clone()),
            mode,
        ).unwrap();

        let final_auxiliary_input_pk_var = RelaxedOvaInstanceVar::new_variable(
            ns!(cs, "final auxiliary input pk var"),
            || verifier_circuit.map(|e| e.final_auxiliary_input_pk.clone()),
            mode,
        ).unwrap();

        Ok(SignatureVerifierCircuitVar {
            pk_1_var,
            pk_2_var,
            pk_t_var,
            auxiliary_input_pk_var,
            running_auxiliary_input_pk_var,
            final_auxiliary_input_pk_var,
            com_pk_var,
            bitfield_poly_var,
            bitfield_poly_commitment,
            sumcheck_proof_var,
            b_1_at_rho,
            b_2_at_rho,
            c_at_rho,
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
        let eq_at_r = MultilinearPolynomialVar::new(
            EqPolynomialVar::new(vec_r).evals()
        );
        let eq_at_r_rho = eq_at_r.evaluate(&sumcheck_challenges);
        FpVar::enforce_equal(
            &tensor_check_claim,
            &(eq_at_r_rho * (self.b_1_at_rho.clone() + self.b_2_at_rho.clone() - self.b_1_at_rho.clone() * self.b_2_at_rho.clone() - self.c_at_rho.clone())),
        ).expect("equality error");

        // Step 4: Do the cycle fold math
        // Non-native scalar multiplication: linear combination E'' = E_{temp} + (1-beta) * beta * Q
        let (flag,
            r,
            pk_1,
            pk_2,
            pk_t
        ) = self.auxiliary_input_pk_var.parse_secondary_io::<G1>().unwrap();
        // g1 == Q
        pk_1.enforce_equal(&self.pk_1_var).expect("error while enforcing equality");
        // g2 == E_temp
        pk_2.enforce_equal(&self.pk_2_var).expect("error while enforcing equality");
        // enforce flag to be true
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        // check r to be equal to one
        r.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        // check out the result E_var is consistent with result_acc
        pk_t.enforce_equal(&self.pk_t_var).expect("error while enforcing equality");

        // Step 5: fold the cycle fold instance
        let one_bits = <NonNativeFieldVar<G1::BaseField, F> as ToBitsGadget<F>>::to_bits_le(&NonNativeFieldVar::one()).unwrap();
        let final_instance = self.running_auxiliary_input_pk_var.fold(
            &[((&self.auxiliary_input_pk_var, None), &self.com_pk_var, &NonNativeFieldVar::one(), &one_bits)]
        ).unwrap();

        self.final_auxiliary_input_pk_var.X.enforce_equal(&final_instance.X).expect("XXX: panic message");
        self.final_auxiliary_input_pk_var.commitment.enforce_equal(&final_instance.commitment).expect("XXX: panic message");
    }
}
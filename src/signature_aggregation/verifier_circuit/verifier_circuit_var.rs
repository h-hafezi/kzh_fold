use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit_var::SumcheckCircuitVar;
use crate::nova::cycle_fold::coprocessor_constraints::{OvaInstanceVar, RelaxedOvaInstanceVar};
use crate::polynomial::eq_poly::eq_poly_var::EqPolynomialVar;
use crate::polynomial::multilinear_poly::multilinear_poly_var::MultilinearPolynomialVar;
use crate::signature_aggregation::verifier_circuit::verifier_circuit::SignatureVerifierCircuit;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::ToBitsGadget;
use ark_relations::ns;
use ark_relations::r1cs::{Namespace, SynthesisError};
use std::borrow::Borrow;

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
    running_pk: NonNativeAffineVar<G1>,
    current_pk: NonNativeAffineVar<G1>,
    final_pk: NonNativeAffineVar<G1>,

    /// randomness for homomorphically combining bitfield commitments
    pub c_0: FpVar<F>,
    pub c_0_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,
    pub c_1: FpVar<F>,
    pub c_1_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// bitfield commitment, non-native only used for fiat shamir and computing random combination P
    pub com_bitfield_C: NonNativeAffineVar<G1>,
    pub com_bitfield_B_1: NonNativeAffineVar<G1>,
    pub com_bitfield_B_2: NonNativeAffineVar<G1>,

    /// com_homomorphic_bitfield = com_bitfield_B_1 + c_0 * com_bitfield_B_2 + c1 * com_bitfield_C
    pub com_homomorphic_bitfield: NonNativeAffineVar<G1>,

    /// beta
    pub beta: FpVar<F>,
    pub beta_non_native: NonNativeFieldVar<G1::BaseField, F>,

    /// adding pk
    pub ova_cross_term_error_pk: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub ova_auxiliary_input_pk: OvaInstanceVar<G2, C2>,

    /// adding temp = B1 + c_0 * B2
    pub ova_cross_term_error_bitfield_1: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub ova_auxiliary_input_bitfield_1: OvaInstanceVar<G2, C2>,

    /// adding res = temp + c1 * C
    pub ova_cross_term_error_bitfield_2: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub ova_auxiliary_input_bitfield_2: OvaInstanceVar<G2, C2>,

    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub ova_running_instance: RelaxedOvaInstanceVar<G2, C2>,
    pub ova_final_instance: RelaxedOvaInstanceVar<G2, C2>,

    /// the sumcheck proof
    sumcheck_proof: SumcheckCircuitVar<F>,

    /// Evaluations of the inner polynomials at rho:
    b_1_at_rho: FpVar<F>,
    b_2_at_rho: FpVar<F>,
    c_at_rho: FpVar<F>,

    /// size of the bitfield
    pub bitfield_num_variables: usize,
}

impl<F, G1, G2, C2> AllocVar<SignatureVerifierCircuit<F, G1, G2, C2>, F> for SignatureVerifierCircuitVar<F, G1, G2, C2>
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
{
    fn new_variable<T: Borrow<SignatureVerifierCircuit<F, G1, G2, C2>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let res = f();
        let circuit = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        // allocate public keys
        let running_pk = NonNativeAffineVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.running_pk.clone()),
            mode,
        )?;
        let current_pk = NonNativeAffineVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.current_pk.clone()),
            mode,
        )?;
        let final_pk = NonNativeAffineVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.final_pk.clone()),
            mode,
        )?;


        // allocate bitfield commitments
        let com_bitfield_C = NonNativeAffineVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.com_bitfield_C.clone()),
            mode,
        )?;
        let com_bitfield_B_1 = NonNativeAffineVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.com_bitfield_B_1.clone()),
            mode,
        )?;
        let com_bitfield_B_2 = NonNativeAffineVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.com_bitfield_B_2.clone()),
            mode,
        )?;
        let com_homomorphic_bitfield = NonNativeAffineVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.com_homomorphic_bitfield.clone()),
            mode,
        )?;


        // allocate beta and beta_non_native
        let beta = FpVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.beta.clone()),
            mode,
        )?;

        let beta_non_native = NonNativeFieldVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.beta_non_native.clone()),
            mode,
        )?;


        // allocate ova instance for public keys
        let ova_cross_term_error_pk = ProjectiveVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_cross_term_error_pk.clone()),
            mode,
        )?;
        let ova_auxiliary_input_pk = OvaInstanceVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_auxiliary_input_pk.clone()),
            mode,
        )?;


        // allocate ova instance for bitfield 1
        let ova_cross_term_error_bitfield_1 = ProjectiveVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_cross_term_error_bitfield_1.clone()),
            mode,
        )?;
        let ova_auxiliary_input_bitfield_1 = OvaInstanceVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_auxiliary_input_bitfield_1.clone()),
            mode,
        )?;

        // allocate ova instance for bitfield 2
        let ova_cross_term_error_bitfield_2 = ProjectiveVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_cross_term_error_bitfield_2.clone()),
            mode,
        )?;
        let ova_auxiliary_input_bitfield_2 = OvaInstanceVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_auxiliary_input_bitfield_2.clone()),
            mode,
        )?;

        // allocate ova running and final instance
        let ova_running_instance = RelaxedOvaInstanceVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_running_instance.clone()),
            mode,
        )?;
        let ova_final_instance = RelaxedOvaInstanceVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.ova_final_instance.clone()),
            mode,
        )?;


        // allocate coefficients c0 and c1 used to compute homomorphic commitments
        let c_0 = FpVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.c_0.clone()),
            mode,
        )?;

        let c_0_non_native = NonNativeFieldVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.c_0_non_native.clone()),
            mode,
        )?;

        let c_1 = FpVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.c_1.clone()),
            mode,
        )?;

        let c_1_non_native = NonNativeFieldVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.c_1_non_native.clone()),
            mode,
        )?;


        // allocate sumcheck proof
        let sumcheck_proof = SumcheckCircuitVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.sumcheck_proof.clone()),
            mode,
        )?;


        // allocate b_1_at_rho, b_2_at_rho, c_at_rho
        let b_1_at_rho = FpVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.b_1_at_rho.clone()),
            mode,
        )?;
        let b_2_at_rho = FpVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.b_2_at_rho.clone()),
            mode,
        )?;
        let c_at_rho = FpVar::new_variable(
            cs.clone(),
            || circuit.map(|e| e.c_at_rho.clone()),
            mode,
        )?;

        let bitfield_num_variables = circuit.map(|e| e.bitfield_num_variables).unwrap();

        // Pass each variable into the struct
        Ok(SignatureVerifierCircuitVar {
            running_pk,
            current_pk,
            final_pk,
            com_bitfield_C,
            com_bitfield_B_1,
            com_bitfield_B_2,
            com_homomorphic_bitfield,
            beta,
            beta_non_native,
            ova_cross_term_error_pk,
            ova_auxiliary_input_pk,
            ova_cross_term_error_bitfield_1,
            ova_auxiliary_input_bitfield_1,
            ova_cross_term_error_bitfield_2,
            ova_auxiliary_input_bitfield_2,
            ova_running_instance,
            ova_final_instance,
            c_0,
            c_0_non_native,
            c_1,
            c_1_non_native,
            sumcheck_proof,
            b_1_at_rho,
            b_2_at_rho,
            c_at_rho,
            bitfield_num_variables,
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
        transcript.append_scalars(b"poly", self.com_bitfield_C.to_sponge_field_elements().unwrap().as_slice());
        transcript.append_scalars(b"poly", self.com_bitfield_B_1.to_sponge_field_elements().unwrap().as_slice());
        transcript.append_scalars(b"poly", self.com_bitfield_B_2.to_sponge_field_elements().unwrap().as_slice());

        let vec_r = transcript.challenge_vector(b"vec_r", self.bitfield_num_variables);

        // Step 2: Verify the sumcheck proof
        let zero: FpVar<F> = FpVar::zero();

        // assert the sumcheck proof is indeed well-formatted
        self.sumcheck_proof.claim.enforce_equal(&zero).expect("equality error");
        assert_eq!(self.sumcheck_proof.num_rounds, self.bitfield_num_variables);
        assert_eq!(self.sumcheck_proof.degree_bound, 3);

        let (tensor_check_claim, sumcheck_challenges) = self.sumcheck_proof.verify(transcript);

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
            g1,
            g2,
            g_out,
        ) = self.ova_auxiliary_input_pk.parse_secondary_io::<G1>().unwrap();
        g1.enforce_equal(&self.running_pk).expect("error while enforcing equality");
        g2.enforce_equal(&self.current_pk).expect("error while enforcing equality");
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        r.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        g_out.enforce_equal(&self.final_pk).expect("error while enforcing equality");

        let vec_c = transcript.challenge_vector(b"vec_r", 2);

        // enforce it's equal to c_0 and c_1
        vec_c[0].enforce_equal(&self.c_0).expect("error while enforcing equality");
        vec_c[1].enforce_equal(&self.c_1).expect("error while enforcing equality");

        // now enforce that the non-native version of c0 and c1 are correct
        self.c_0.enforce_equal(&{
            let bits = self.c_0.to_bits_le().unwrap();
            Boolean::le_bits_to_fp_var(bits.as_slice()).unwrap()
        }).expect("error while enforcing equality");
        self.c_1.enforce_equal(&{
            let bits = self.c_1.to_bits_le().unwrap();
            Boolean::le_bits_to_fp_var(bits.as_slice()).unwrap()
        }).expect("error while enforcing equality");

        // derive beta
        let beta = transcript.challenge_scalar(b"beta");

        // enforce it's consistent with the original challenge beta
        beta.enforce_equal(&self.beta).expect("error while enforcing equality");

        // Non-native scalar multiplication: linear combination temp = B1 + c_0 * B2
        let (flag,
            r,
            g1,
            g2,
            temp
        ) = self.ova_auxiliary_input_bitfield_1.parse_secondary_io::<G1>().unwrap();
        g1.enforce_equal(&self.com_bitfield_B_2).expect("error while enforcing equality");
        g2.enforce_equal(&self.com_bitfield_B_1).expect("error while enforcing equality");
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        r.enforce_equal(&self.c_0_non_native).expect("error while enforcing equality");

        let (flag,
            r,
            g1,
            g2,
            g_out
        ) = self.ova_auxiliary_input_bitfield_2.parse_secondary_io::<G1>().unwrap();
        g1.enforce_equal(&self.com_bitfield_C).expect("error while enforcing equality");
        g2.enforce_equal(&temp).expect("error while enforcing equality");
        flag.enforce_equal(&NonNativeFieldVar::one()).expect("error while enforcing equality");
        r.enforce_equal(&self.c_1_non_native).expect("error while enforcing equality");
        g_out.enforce_equal(&self.com_homomorphic_bitfield).expect("error while enforcing equality");


        // Step 5: fold the cycle fold instance
        let beta_bits = <NonNativeFieldVar<G1::BaseField, F> as ToBitsGadget<F>>::to_bits_le(&self.beta_non_native).unwrap();
        let final_instance = self.ova_running_instance.fold(
            &[
                (
                    (&self.ova_auxiliary_input_pk, None),
                    &self.ova_cross_term_error_pk,
                    &self.beta_non_native,
                    &beta_bits
                ),
                (
                    (&self.ova_auxiliary_input_bitfield_1, None),
                    &self.ova_cross_term_error_bitfield_1,
                    &self.beta_non_native,
                    &beta_bits
                ),
                (
                    (&self.ova_auxiliary_input_bitfield_2, None),
                    &self.ova_cross_term_error_bitfield_2,
                    &self.beta_non_native,
                    &beta_bits
                ),
            ]
        ).unwrap();

        self.ova_final_instance.X.enforce_equal(&final_instance.X).expect("XXX: panic message");
        self.ova_final_instance.commitment.enforce_equal(&final_instance.commitment).expect("XXX: panic message");
    }
}


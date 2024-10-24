use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::gadgets::r1cs::OvaInstance;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit_var::SumcheckCircuitVar;
use crate::nova::cycle_fold::coprocessor_constraints::OvaInstanceVar;
use crate::polynomial::eq_poly::eq_poly_var::EqPolynomialVar;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::polynomial::multilinear_poly::multilinear_poly_var::MultilinearPolynomialVar;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::CurveConfig;
use ark_ff::PrimeField;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;

pub struct SignatureVerifierCircuit<E, F, Primary, Secondary, TwistedPrimary, C2>
where
// primary is the primary curve and G1
    Primary: SWCurveConfig + Clone,
    Primary::BaseField: PrimeField,
    Primary::ScalarField: PrimeField,
// secondary and primary form a cycle of curves
    Secondary: SWCurveConfig,
    Secondary::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<Secondary>>,
// primary and secondary base/scalar field must match
    Primary: SWCurveConfig<
        BaseField=Secondary::ScalarField,
        ScalarField=Secondary::BaseField
    >,
// Pairing curves are Primary and Twisted Primary
    TwistedPrimary: SWCurveConfig,
    E: Pairing<
        G1Affine=Affine<Primary>,
        ScalarField=<Primary as CurveConfig>::ScalarField,
        G2Affine=Affine<TwistedPrimary>,
    >,
// the scalar field of E must satisfy Absorb trait
    F: PrimeField + Absorb,
{
    /// public keys pk_t = pk_1 + pk_2
    pk_1: E::G1Affine,
    pk_2: E::G1Affine,
    pk_t: E::G1Affine,

    /// signatures sig_t = sig_1 + sig_2
    sig_1: E::G2Affine,
    sig_2: E::G2Affine,
    sig_t: E::G2Affine,

    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub auxiliary_input_pk: OvaInstance<Secondary, C2>,
    /// auxiliary input which helps to have sig_t = sig_2 + sig_1
    pub auxiliary_input_sig: OvaInstance<Secondary, C2>,

    /// the bitfield polynomial
    bitfield_poly: MultilinearPolynomial<F>,

    /// the sumcheck proof
    sumcheck_proof: SumcheckCircuit<F>,

    /// Evaluations of the inner polynomials at rho:
    b_1_at_rho: F,
    b_2_at_rho: F,
    c_at_rho: F,
}

pub struct SignatureVerifierCircuitVar<F, Primary, Secondary, TwistedPrimary, C2>
where
    F: PrimeField + Absorb,
    Primary: SWCurveConfig<ScalarField=F> + Clone,
    Primary::BaseField: PrimeField,
    Primary::ScalarField: PrimeField,
    Secondary: SWCurveConfig,
    Secondary::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<Secondary>>,
// the secondary and primary must form a cycle of curves
    Secondary: SWCurveConfig<
        BaseField=Primary::ScalarField,
        ScalarField=Primary::BaseField
    >,
// the primary and twisted primary must have identical base / scalar fields
    TwistedPrimary: SWCurveConfig<
        BaseField=Primary::BaseField,
        ScalarField=Primary::ScalarField,
    >,
{
    /// public keys pk_t = pk_1 + pk_2
    pk_1_var: NonNativeAffineVar<Primary>,
    pk_2: NonNativeAffineVar<Primary>,
    pk_t: NonNativeAffineVar<Primary>,

    /// signatures sig_t = sig_1 + sig_2
    sig_1: NonNativeAffineVar<TwistedPrimary>,
    sig_2: NonNativeAffineVar<TwistedPrimary>,
    sig_t: NonNativeAffineVar<TwistedPrimary>,

    /// auxiliary input which helps to have pk_t = pk_2 + pk_1
    pub auxiliary_input_pk_var: OvaInstanceVar<Secondary, C2>,
    /// auxiliary input which helps to have sig_t = sig_2 + sig_1
    pub auxiliary_input_sig_var: OvaInstanceVar<Secondary, C2>,

    /// the bitfield polynomial
    bitfield_poly: MultilinearPolynomialVar<F>,

    /// the sumcheck proof
    sumcheck_proof: SumcheckCircuitVar<F>,

    /// Evaluations of the inner polynomials at rho:
    b_1_at_rho: F,
    b_2_at_rho: F,
    c_at_rho: F,
}

impl<F, Primary, Secondary, TwistedPrimary, C2> SignatureVerifierCircuitVar<F, Primary, Secondary, TwistedPrimary, C2>
where
    F: PrimeField + Absorb,
    Primary: SWCurveConfig<ScalarField=F> + Clone,
    Primary::BaseField: PrimeField,
    Primary::ScalarField: PrimeField,
    Secondary: SWCurveConfig,
    Secondary::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<Secondary>>,
// the secondary and primary must form a cycle of curves
    Secondary: SWCurveConfig<
        BaseField=Primary::ScalarField,
        ScalarField=Primary::BaseField
    >,
// the primary and twisted primary must have identical base / scalar fields
    TwistedPrimary: SWCurveConfig<
        BaseField=Primary::BaseField,
        ScalarField=Primary::ScalarField,
    >,
{
    pub fn verify(&self, transcript: &mut TranscriptVar<F>) {
        // Step 1: Get challenge
        let vec_r = transcript.challenge_vector(b"vec_r", self.bitfield_poly.num_variables);

        // Step 2: Verify the sumcheck proof
        let zero: FpVar<F> = FpVar::zero();
        let num_rounds = self.bitfield_poly.num_variables;

        // assert the sumcheck proof is indeed well-formatted
        self.sumcheck_proof.claim.enforce_equal(&zero).expect("equality error");
        assert_eq!(self.sumcheck_proof.num_rounds, self.bitfield_poly.num_variables);
        assert_eq!(self.sumcheck_proof.degree_bound, 3);

        let (tensor_check_claim, sumcheck_challenges) = self.sumcheck_proof.verify(transcript);

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
            &(eq_at_r_rho * (self.b_1_at_rho + self.b_2_at_rho - self.b_1_at_rho * self.b_2_at_rho - self.c_at_rho))
        ).expect("equality error");
    }
}
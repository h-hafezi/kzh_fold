/// Accumulation verifier implementation
use std::marker::PhantomData;

use ark_crypto_primitives::crh::{
    poseidon::constraints::{CRHGadget, CRHParametersVar},
    CRHSchemeGadget,
};
use ark_ec::{pairing::Pairing, short_weierstrass::{Projective, SWCurveConfig}};
use ark_ec::AffineRepr;
use ark_ec::VariableBaseMSM;
use ark_ff::PrimeField;
use ark_relations::{r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError}};
use ark_r1cs_std::{
    alloc::AllocVar,
    eq::EqGadget,
    fields::fp::FpVar,
    groups::{curves::short_weierstrass::ProjectiveVar, CurveVar},
    ToBitsGadget,
};

use crate::{pcs::{Commitment, OpeningProof, SRS}, univariate_poly::UnivariatePolynomial, utils::compute_powers};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccInstance<E: Pairing> {
    C: Commitment<E>,
    T: E::G1Affine,
    b: E::ScalarField,
    c: E::ScalarField,
    y: E::ScalarField,
    E: E::G1Affine,
}

/// AccVerifierSNARKCircuit takes as input A.X, A_1.X, Q and folds them.
/// We assume that we are using a cycle of curves or cyclefold here.
pub struct AccVerifierSNARKCircuit<G1: SWCurveConfig> {
    // XXX should this be BaseField? That's what
    // https://github.com/nexus-xyz/nexus-zkvm/blob/7f7789a271a4d7950ad6b4347cb6190b589c15c2/nova/src/folding/cyclefold/secondary/mod.rs#L39
    // is doing
    // We assume that the fold challenge is being taken as public input from the main curve (?)
    pub beta: G1::BaseField,
    pub E_1: Projective<G1>,
    // pub E_2: Projective<G1>,
    // pub Q: Projective<G1>,
    // All the rest...
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccVerifier<E: Pairing> {
    _phantom: PhantomData<E>,
}

impl<G1: SWCurveConfig> ConstraintSynthesizer<G1::BaseField> for AccVerifierSNARKCircuit<G1>
where
    G1::BaseField: PrimeField,
{
    /// Generate the constraints for this circuit
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<G1::BaseField>
    ) -> Result<(), SynthesisError> {
        let E_1 = ProjectiveVar::<G1, FpVar<G1::BaseField>>::new_input(cs.clone(), || Ok(self.E_1))?;
        let beta = FpVar::<G1::BaseField>::new_input(cs.clone(), || Ok(self.beta))?;
        let beta_bits = beta.to_bits_le()?;

        // XXX Just a demonstration of how to do a **native** scalar mul
        let _tmp = E_1.scalar_mul_le(beta_bits.iter())?;
        // out.enforce_equal(g_out);

        unimplemented!();
    }
}

impl<E: Pairing> AccVerifier<E> {
    fn new_acc_instance_from_proof(C: &Commitment<E>, proof: &OpeningProof<E>, b: &E::ScalarField, c: &E::ScalarField, y: &E::ScalarField, srs: &SRS<E>) -> AccInstance<E> {
        let n = proof.vec_D_i.len();
        let vec_b_powers = compute_powers(b, n);
        let vec_c_powers = compute_powers(c, n);

        let mut vec_scalars = vec_b_powers.clone();
        vec_scalars.extend(vec_c_powers);
        vec_scalars.push(*y);

        let T = E::G1::msm_unchecked(&srs.vec_G_accumulation, &vec_scalars);

        AccInstance {
            C: C.clone(),
            T: T.clone().into(),
            b: *b,
            c: *c,
            y: *y,
            E: E::G1Affine::zero(),
        }
    }

    /// Accumulation verifier: Given (A_1.X, A_2.X) compute A_3.X
    fn verify(_instance_1: &AccInstance<E>, _instance_2: &AccInstance<E>, _quotient: E::G1Affine) -> AccInstance<E> {
        unimplemented!();
    }
}

#[cfg(test)]
pub mod tests {
    use std::ops::Mul;

    use super::*;
    use ark_r1cs_std::{alloc::AllocVar, fields::{fp::FpVar}};
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;

    /// Test that the accumulation verifier circuit correctly folds the instances
    #[test]
    fn test_circuit() {
        let _rng = &mut rand::thread_rng();

        unimplemented!();
    }
}

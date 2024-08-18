//! Helper definitions for a conditional elliptic curve operation circuit.
//!
//! This circuit accepts `g1`, `g2`, `g_out`, `r`, and `flag` (in this exact order) as its public input, where
//! `g1`, `g2`, `g_out` are points on the elliptic curve G, `r` is an element from the scalar field, and `flag`
//! is a boolean value. The circuit performs one of two possible operations based on the value of `flag`:
//!
//! - If `flag` is `true`, the circuit enforces `g_out = r * g1 + g2`.
//! - If `flag` is `false`, the circuit enforces `g_out = r * (g1 - g2) + g2 = r * g1 + (1-r) * g2`.
//!
//! The circuit efficiently computes the result by conditionally selecting between `g1` and `g1 - g2` before
//! performing a single scalar multiplication and adding `g2`. This ensures that the circuit minimizes the
//! number of expensive scalar multiplication operations.

use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::{AdditiveGroup, PrimeField};
use ark_r1cs_std::{
    alloc::AllocVar,
    eq::EqGadget,
    fields::fp::FpVar,
    groups::{curves::short_weierstrass::ProjectiveVar, CurveVar},
    ToBitsGadget,
};
use ark_r1cs_std::boolean::Boolean;
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError, SynthesisMode,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::Zero;

use crate::gadgets::r1cs::*;
use crate::commitment::CommitmentScheme;

/// Leading One + 3 curve points + 1 scalar + 1 flag.
const SECONDARY_NUM_IO: usize = 12;

/// Public input of secondary circuit.
#[derive(Debug, PartialEq, Eq)]
pub struct SecondaryCircuit<G1: SWCurveConfig> {
    pub(crate) g1: Projective<G1>,
    pub(crate) g2: Projective<G1>,
    pub(crate) g_out: Projective<G1>,
    /// Scalar for elliptic curve points multiplication is part of the public
    /// input and hence should fit into the base field of G1.
    ///
    /// See [`super::nimfs::SQUEEZE_ELEMENTS_BIT_SIZE`].
    pub(crate) r: G1::BaseField,
    pub(crate) flag: bool,
}

impl<G1: SWCurveConfig> SecondaryCircuit<G1> {
    pub const NUM_IO: usize = SECONDARY_NUM_IO;
}

impl<G: SWCurveConfig> Default for SecondaryCircuit<G> {
    fn default() -> Self {
        Self {
            g1: Projective::zero(),
            g2: Projective::zero(),
            g_out: Projective::zero(),
            r: G::BaseField::ZERO,
            flag: false,
        }
    }
}

impl<G: SWCurveConfig> Clone for SecondaryCircuit<G> {
    fn clone(&self) -> Self {
        Self {
            g1: self.g1,
            g2: self.g2,
            g_out: self.g_out,
            r: self.r,
            flag: self.flag,
        }
    }
}

impl<G1: SWCurveConfig> ConstraintSynthesizer<G1::BaseField> for SecondaryCircuit<G1>
where
    G1::BaseField: PrimeField,
{
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<G1::BaseField>,
    ) -> Result<(), SynthesisError> {
        // Allocate inputs
        let g1 = ProjectiveVar::<G1, FpVar<G1::BaseField>>::new_input(cs.clone(), || Ok(self.g1))?;
        let g2 = ProjectiveVar::<G1, FpVar<G1::BaseField>>::new_input(cs.clone(), || Ok(self.g2))?;
        let g_out = ProjectiveVar::<G1, FpVar<G1::BaseField>>::new_input(cs.clone(), || Ok(self.g_out))?;
        let r = FpVar::<G1::BaseField>::new_input(cs.clone(), || Ok(self.r))?;
        // It's automatically converted to either 0 or 1 on G1::BaseField
        let flag = Boolean::new_input(cs.clone(), || Ok(self.flag))?;

        // Compute A = g1 and A' = g1 - g2
        let a = g1.clone();
        let a_prime = g1 - &g2;

        // Conditionally select B = A if flag is true, else B = A'
        let b = flag.select(&a, &a_prime)?;

        // Compute B * r + g2
        let r_bits = r.to_bits_le()?;
        let result = b.scalar_mul_le(r_bits.iter())? + &g2;

        // Enforce the result to be equal to g_out
        result.enforce_equal(&g_out)?;

        Ok(())
    }
}

/// Setup [`R1CSShape`] for a secondary circuit, defined over `G2::BaseField`.
pub fn setup_shape<G1, G2>() -> Result<R1CSShape<G2>, SynthesisError>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
    G2: SWCurveConfig<BaseField=G1::ScalarField, ScalarField=G1::BaseField>,
{
    let cs = ConstraintSystem::<G1::BaseField>::new_ref();
    cs.set_mode(SynthesisMode::Setup);

    SecondaryCircuit::<G1>::default().generate_constraints(cs.clone())?;

    cs.finalize();
    Ok(R1CSShape::from(cs.clone()))
}

/// Synthesize public input and a witness-trace.
pub fn synthesize<G1, G2, C2>(
    circuit: SecondaryCircuit<G1>,
    pp_secondary: &C2::PP,
) -> Result<(R1CSInstance<G2, C2>, R1CSWitness<G2>), SynthesisError>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
    G2: SWCurveConfig<BaseField=G1::ScalarField, ScalarField=G1::BaseField>,
    C2: CommitmentScheme<Projective<G2>>,
{
    let cs = ConstraintSystem::<G1::BaseField>::new_ref();
    cs.set_mode(SynthesisMode::Prove { construct_matrices: false });

    circuit.generate_constraints(cs.clone())?;

    cs.finalize();
    let cs_borrow = cs.borrow().unwrap();

    let witness = cs_borrow.witness_assignment.clone();
    let pub_io = cs_borrow.instance_assignment.clone();

    let W = R1CSWitness::<G2> { W: witness };

    let commitment_W = W.commit::<C2>(pp_secondary);
    let U = R1CSInstance::<G2, C2> { commitment_W, X: pub_io };

    Ok((U, W))
}

#[cfg(any(test, feature = "spartan"))]
macro_rules! parse_projective {
    ($X:expr) => {
        match &$X[..3] {
            &[x, y, z, ..] => {
                let point = ark_ec::CurveGroup::into_affine(Projective::<G1> { x, y, z });
                if !point.is_on_curve() || !point.is_in_correct_subgroup_assuming_on_curve() {
                    return None;
                }
                $X = &$X[3..];
                point.into()
            }
            _ => return None,
        }
    };
}

impl<G2, C2> R1CSInstance<G2, C2>
where
    G2: SWCurveConfig,
    C2: CommitmentScheme<Projective<G2>>,
{
    #[cfg(any(test, feature = "spartan"))]
    pub(crate) fn parse_secondary_io<G1>(&self) -> Option<SecondaryCircuit<G1>>
    where
        G2::BaseField: PrimeField,
        G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    {
        let mut X = &self.X[1..];

        let g1 = parse_projective!(X);
        let g2 = parse_projective!(X);
        let g_out = parse_projective!(X);

        if X.len() < 2 {
            return None;
        }

        let r = X[0];
        let flag = X[1];

        Some(SecondaryCircuit { g1, g2, g_out, r, flag: !flag.is_zero() })
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{Field, PrimeField};
    use ark_pallas::{Fq, Fr, PallasConfig, Projective};
    use ark_std::UniformRand;
    use ark_vesta::VestaConfig;

    use crate::hash::pederson::PedersenCommitment;
    use crate::utils::cast_field_element;

    use super::*;

    #[test]
    pub fn parse_pub_input() {
        let mut rng = ark_std::test_rng();
        let g1 = Projective::rand(&mut rng);
        let g2 = Projective::rand(&mut rng);

        let val = u64::rand(&mut rng);
        let r = <Fq as PrimeField>::BigInt::from(val).into();
        let r_scalar = unsafe { cast_field_element::<Fq, Fr>(&r) };

        let g_out = g1 * r_scalar + g2;

        let expected_pub_io = SecondaryCircuit::<PallasConfig> { g1, g2, g_out, r, flag: true };
        let X = [
            Fq::ONE,
            g1.x,
            g1.y,
            g1.z,
            g2.x,
            g2.y,
            g2.z,
            g_out.x,
            g_out.y,
            g_out.z,
            unsafe { cast_field_element(&r) },
            Fq::ONE,
        ];

        assert_eq!(X.len(), SECONDARY_NUM_IO);

        let r1cs = R1CSInstance::<VestaConfig, PedersenCommitment<ark_vesta::Projective>> {
            commitment_W: Default::default(),
            X: X.into(),
        };

        let pub_io = r1cs.parse_secondary_io().unwrap();

        assert_eq!(pub_io.g1, expected_pub_io.g1);
        assert_eq!(pub_io.g2, expected_pub_io.g2);
        assert_eq!(pub_io.g_out, expected_pub_io.g_out);
        assert_eq!(pub_io.r, expected_pub_io.r);
        assert_eq!(pub_io.flag, expected_pub_io.flag);

        // incorrect length
        let _X = &X[..10];
        let r1cs = R1CSInstance::<VestaConfig, PedersenCommitment<ark_vesta::Projective>> {
            commitment_W: Default::default(),
            X: _X.into(),
        };
        assert!(r1cs.parse_secondary_io::<PallasConfig>().is_none());

        // not on curve
        let mut _X = X.to_vec();
        _X[1] -= Fq::ONE;
        let r1cs = R1CSInstance::<VestaConfig, PedersenCommitment<ark_vesta::Projective>> {
            commitment_W: Default::default(),
            X: _X,
        };
        assert!(r1cs.parse_secondary_io::<PallasConfig>().is_none());
    }

    #[test]
    pub fn parse_synthesized() {
        let shape = setup_shape::<PallasConfig, VestaConfig>().unwrap();
        let mut rng = ark_std::test_rng();
        let g1 = Projective::rand(&mut rng);
        let g2 = Projective::rand(&mut rng);

        let val = u64::rand(&mut rng);
        let r = <Fq as PrimeField>::BigInt::from(val).into();
        let r_scalar = unsafe { cast_field_element::<Fq, Fr>(&r) };

        let g_out = g1 * r_scalar + g2;

        let pp = PedersenCommitment::<ark_vesta::Projective>::setup(shape.num_vars, b"test", &());
        let (U, _) = synthesize::<
            PallasConfig,
            VestaConfig,
            PedersenCommitment<ark_vesta::Projective>,
        >(SecondaryCircuit { g1, g2, g_out, r, flag: false }, &pp).unwrap();

        let pub_io = U.parse_secondary_io::<PallasConfig>().unwrap();

        assert_eq!(pub_io.g1, g1);
        assert_eq!(pub_io.g2, g2);
        assert_eq!(pub_io.g_out, g_out);
        assert_eq!(pub_io.r, r);
        assert_eq!(pub_io.flag, false);
    }

    #[test]
    pub fn satisfiability_test_for_true_flag() {
        let mut rng = ark_std::test_rng();
        let g1 = Projective::rand(&mut rng);
        let g2 = Projective::rand(&mut rng);

        let r = <Fq as PrimeField>::BigInt::from(u64::rand(&mut rng)).into();
        let r_scalar = unsafe { cast_field_element::<Fq, Fr>(&r) };

        let g_out = g1 * r_scalar + g2;

        let c = SecondaryCircuit { g1, g2, g_out, r, flag: true };

        let cs = ConstraintSystem::new_ref();
        c.generate_constraints(cs.clone()).unwrap();

        println!("{}", cs.num_constraints());
        println!("{}", cs.is_satisfied().unwrap());
    }

    #[test]
    pub fn satisfiability_test_for_false_flag() {
        let mut rng = ark_std::test_rng();
        let g1 = Projective::rand(&mut rng);
        let g2 = Projective::rand(&mut rng);

        let r = <Fq as PrimeField>::BigInt::from(u64::rand(&mut rng)).into();
        let r_scalar = unsafe { cast_field_element::<Fq, Fr>(&r) };

        let g_out = g1 * r_scalar + g2 * (Fr::ONE - r_scalar);

        let c = SecondaryCircuit { g1, g2, g_out, r, flag: false };

        let cs = ConstraintSystem::new_ref();
        c.generate_constraints(cs.clone()).unwrap();

        println!("{}", cs.num_constraints());
        println!("{}", cs.is_satisfied().unwrap());
    }
}
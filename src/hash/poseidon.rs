use ark_crypto_primitives::crh::poseidon::constraints::CRHParametersVar;
/// Poseidon config stolen from sonobe

use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge, poseidon::{PoseidonConfig}};
use ark_crypto_primitives::sponge::constraints::{AbsorbGadget, CryptographicSpongeVar};
use ark_crypto_primitives::sponge::poseidon::{find_poseidon_ark_and_mds, PoseidonSponge};
use ark_crypto_primitives::sponge::poseidon::constraints::PoseidonSpongeVar;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_relations::r1cs::ConstraintSystemRef;

pub struct PoseidonHash<F: Absorb + PrimeField> {
    pub(crate) poseidon_params: PoseidonConfig<F>,
    sponge: PoseidonSponge<F>,
}

pub fn get_poseidon_config<F: PrimeField>() -> PoseidonConfig<F> {
    // 120 bit security target as in
    // https://eprint.iacr.org/2019/458.pdf
    // t = rate + 1

    let full_rounds = 8;
    let partial_rounds = 60;
    let alpha = 5;
    let rate = 4;

    let (ark, mds) = find_poseidon_ark_and_mds::<F>(
        F::MODULUS_BIT_SIZE as u64,
        rate,
        full_rounds,
        partial_rounds,
        0,
        );
    PoseidonConfig::new(
        full_rounds as usize,
        partial_rounds as usize,
        alpha,
        mds,
        ark,
        rate,
        1,
    )
}

impl<F: Absorb + PrimeField> PoseidonHash<F> {
    /// This Poseidon configuration generator agrees with Circom's Poseidon(4) in the case of BN254's scalar field
    pub fn new(poseidon_params: &PoseidonConfig<F>) -> Self {
        Self {
            poseidon_params: poseidon_params.clone(),
            sponge: PoseidonSponge::new(&poseidon_params),
        }
    }

    pub fn update_sponge<A: Absorb>(&mut self, field_vector: Vec<A>) -> () {
        for field_element in field_vector {
            self.sponge.absorb(&field_element);
        }
    }

    pub fn output(&mut self) -> F {
        let squeezed_field_element: Vec<F> = self.sponge.squeeze_field_elements(1);
        squeezed_field_element[0]
    }
}

pub struct PoseidonHashVar<F: Absorb + PrimeField> {
    poseidon_params: CRHParametersVar<F>,
    sponge: PoseidonSpongeVar<F>,
}

impl<F: Absorb + PrimeField> PoseidonHashVar<F> {
    pub fn new(cs: ConstraintSystemRef<F>) -> Self {
        let poseidon_config = get_poseidon_config();
        let hash = PoseidonHash::new(&poseidon_config);
        // XXX: later don't clone
        let poseidon_params = CRHParametersVar::<F>::new_witness(cs.clone(), || Ok(hash.poseidon_params.clone())).unwrap();
        let sponge = PoseidonSpongeVar::new(cs, &hash.poseidon_params);
        PoseidonHashVar {
            poseidon_params,
            sponge,
        }
    }

    pub fn update_sponge<A: AbsorbGadget<F>>(&mut self, field_vector: Vec<A>) -> () {
        for field_element in field_vector {
            self.sponge.absorb(&field_element).expect("Error while sponge absorbing");
        }
    }

    pub fn output(&mut self) -> FpVar<F> {
        let squeezed_field_element: Vec<FpVar<F>> = self.sponge.squeeze_field_elements(1).unwrap();
        squeezed_field_element[0].clone()
    }
}


#[cfg(test)]
mod tests {
    use std::ops::Mul;
    use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
    use ark_crypto_primitives::sponge::poseidon::{PoseidonConfig, PoseidonSponge};
    use ark_ec::CurveGroup;
    use ark_ec::short_weierstrass::{SWCurveConfig, Projective};
    use ark_ff::{Field, PrimeField};
    use ark_r1cs_std::alloc::{AllocationMode, AllocVar};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
    use ark_r1cs_std::{R1CSVar, ToConstraintFieldGadget};
    use ark_r1cs_std::fields::FieldVar;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
    use ark_std::UniformRand;
    use rand::rngs::OsRng;
    use rand::thread_rng;

    use super::*;

    use crate::constant_for_curves::{BaseField, ScalarField};
    use crate::gadgets::non_native::util::{convert_field_one_to_field_two, non_native_to_fpvar};


    #[test]
    fn hash_test() {
        let poseidon_config = get_poseidon_config();
        let mut hash_object: PoseidonHash<ScalarField> = PoseidonHash::new(&poseidon_config);
        hash_object.update_sponge(vec![ScalarField::ONE, ScalarField::ONE]);
        hash_object.update_sponge(vec![convert_field_one_to_field_two::<BaseField, ScalarField>(BaseField::ONE)]);

        let cs = ConstraintSystem::new_ref();

        let mut hash_object_var: PoseidonHashVar<ScalarField> = PoseidonHashVar::new(cs.clone());

        let one_on_first_curve_var = FpVar::Constant(ScalarField::ONE);
        hash_object_var.update_sponge(vec![one_on_first_curve_var.clone(), one_on_first_curve_var.clone()]);
        let one_on_second_curve_var = NonNativeFieldVar::Constant(BaseField::ONE);

        hash_object_var.update_sponge(vec![non_native_to_fpvar(&one_on_second_curve_var)]);

        assert_eq!(hash_object_var.output().value().unwrap(), hash_object.output());
    }
}



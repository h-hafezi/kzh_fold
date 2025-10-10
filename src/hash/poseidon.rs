use ark_crypto_primitives::sponge::constraints::{AbsorbGadget, CryptographicSpongeVar};
use ark_crypto_primitives::sponge::poseidon::constraints::PoseidonSpongeVar;
use ark_crypto_primitives::sponge::poseidon::{find_poseidon_ark_and_mds, PoseidonSponge};
use ark_crypto_primitives::sponge::{poseidon::PoseidonConfig, Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_relations::r1cs::ConstraintSystemRef;


#[derive(Clone)]
pub struct PoseidonHash<F: Absorb + PrimeField> {
    pub sponge: PoseidonSponge<F>,
}

pub fn get_poseidon_config<F: PrimeField>() -> PoseidonConfig<F> {
    // 120 bit security target as in
    // https://eprint.iacr.org/2019/458.pdf
    // t = rate + 1

    let full_rounds = 8;
    let partial_rounds = 60;
    let alpha = 5;
    let rate = 32;

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
    pub fn new() -> Self {
        let poseidon_params = get_poseidon_config::<F>();
        Self {
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
    sponge: PoseidonSpongeVar<F>,
}

impl<F: Absorb + PrimeField> PoseidonHashVar<F> {
    pub fn new(cs: ConstraintSystemRef<F>) -> Self {
        let poseidon_params = get_poseidon_config::<F>();

        // get the SpongeVar
        let sponge = PoseidonSpongeVar::new(cs, &poseidon_params);

        PoseidonHashVar {
            sponge,
        }
    }

    /// the function takes a PoseidonHash object and converts it into a
    /// PoseidonVar object by transforming it states from F into FpVar<F>,
    /// but the rest of arguments e.g. parameters and mode, do not change
    pub fn from_poseidon_hash(cs: ConstraintSystemRef<F>, poseidon_hash: PoseidonHash<F>) -> Self {
        // convert state from F into FpVar
        let state = {
            let mut res = Vec::new();
            for i in poseidon_hash.sponge.state {
                res.push(FpVar::new_input(cs.clone(), ||Ok(i)).unwrap());
            }
            res
        };

        PoseidonHashVar{
            sponge: PoseidonSpongeVar {
                cs,
                parameters: poseidon_hash.sponge.parameters.clone(),
                state,
                mode: poseidon_hash.sponge.mode.clone(),
            },
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
    use ark_ff::Field;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;

    use super::*;

    use crate::constant_for_curves::{BaseField, ScalarField};
    use crate::gadgets::non_native::util::{cast_field, non_native_to_fpvar};


    // In this test, we have a Poseidon and PoseidonVar, we feed them on native and
    // non-native data on both, at the end we check if the results are equal

    #[test]
    fn hash_test() {
        let mut hash_object: PoseidonHash<ScalarField> = PoseidonHash::new();
        hash_object.update_sponge(vec![ScalarField::ONE, ScalarField::ONE]);
        hash_object.update_sponge(vec![cast_field::<BaseField, ScalarField>(BaseField::ONE)]);

        let cs = ConstraintSystem::new_ref();

        let mut hash_object_var: PoseidonHashVar<ScalarField> = PoseidonHashVar::new(cs.clone());

        let one_on_first_curve_var = FpVar::Constant(ScalarField::ONE);
        hash_object_var.update_sponge(vec![one_on_first_curve_var.clone(), one_on_first_curve_var.clone()]);
        let one_on_second_curve_var = NonNativeFieldVar::Constant(BaseField::ONE);

        hash_object_var.update_sponge(vec![non_native_to_fpvar(&one_on_second_curve_var)]);

        assert_eq!(hash_object_var.output().value().unwrap(), hash_object.output());

        // test from_poseidon_function
        let mut hash_object_var_new = PoseidonHashVar::from_poseidon_hash(
            ConstraintSystem::new_ref(),
            hash_object.clone(),
        );

        assert_eq!(
            hash_object_var_new.output().value().unwrap(),
            hash_object_var.output().value().unwrap()
        );
    }
}



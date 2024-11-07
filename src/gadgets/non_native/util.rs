use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ff::{AdditiveGroup, BigInteger, Field, PrimeField};
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::ToBitsGadget;

pub fn non_native_to_fpvar<ScalarField, BaseField>(
    non_native_var: &NonNativeFieldVar<BaseField, ScalarField>,
) -> FpVar<ScalarField>
where
    ScalarField: PrimeField,
    BaseField: PrimeField,
{
    // Convert the non-native field variable to a bit vector
    let bits = non_native_var.to_bits_le().unwrap();

    // Convert the bit vector into an FpVar
    let x = Boolean::le_bits_to_fp_var(bits.as_slice()).unwrap();

    x
}

/// given an E::G1Affine, it converts it into a tuple of two scalars by converting each coordinate into a scalar point
/// finally, in case the point is zero is returns (1, 0) which is standard representation of the infinity point
/// it's supposed to be consistent with NonNativeAffineVar::to_sponge_field_elements
pub fn convert_affine_to_scalars<E: Pairing>(point: E::G1Affine) -> (E::ScalarField, E::ScalarField)
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    if point.is_zero() {
        (E::ScalarField::ONE, E::ScalarField::ZERO)
    } else {
        // Extract x and y coordinates and convert them
        let x = cast_field::<
            <<E as Pairing>::G1Affine as AffineRepr>::BaseField,
            <E as Pairing>::ScalarField,
        >(point.x().unwrap());
        let y = cast_field::<
            <<E as Pairing>::G1Affine as AffineRepr>::BaseField,
            <E as Pairing>::ScalarField,
        >(point.y().unwrap());

        // return the converted tuple
        (x, y)
    }
}


pub fn cast_field<Fr, Fq>(first_field: Fr) -> Fq
where
    Fr: PrimeField,
    Fq: PrimeField,
{
    // Convert the Fr element to its big integer representation
    let bytes = first_field.into_bigint().to_bytes_le();

    // Convert the big integer representation to an Fq element
    let fq_element = Fq::from_le_bytes_mod_order(bytes.as_slice());

    fq_element
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{BaseField, ScalarField};
    use crate::gadgets::non_native::util::{cast_field, non_native_to_fpvar};

    // This test makes sure non_native_to_fpvar works correctly but generating a
    // random non-native value and then converting it into FpVar, it's the zk version
    #[test]
    fn test_conversion1() {
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        let g = BaseField::rand(&mut thread_rng());

        // Create the non-native field variable outside the function
        let non_native_var = NonNativeFieldVar::new_variable(
            cs.clone(),
            || Ok(g),
            AllocationMode::Input,
        ).unwrap();

        println!("constraint counts before: {}", cs.num_constraints());

        let x = non_native_to_fpvar(&non_native_var);

        println!("constraint counts after: {}", cs.num_constraints());

        // make sure the non-native to fpvar in fact works correctly
        assert_eq!(g.into_bigint(), x.value().unwrap().into_bigint());
    }

    // equivalent of test above, in non-zk

    #[test]
    fn test_conversion2() {
        let g = BaseField::rand(&mut thread_rng());
        assert_eq!(g.into_bigint(), cast_field::<BaseField, ScalarField>(g).into_bigint());
    }
}
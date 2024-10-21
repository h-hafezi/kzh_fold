use ark_ff::PrimeField;
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::ToBitsGadget;

pub fn get_bits_canonical_order<F: PrimeField>(fp: &FpVar<F>, n: usize) -> Vec<Boolean<F>> {
    let new_vector: Vec<Boolean<F>> = fp.to_bits_be().unwrap();
    new_vector[new_vector.len() - n..].to_vec()
}

#[cfg(test)]
mod tests {
    use crate::constant_for_curves::ScalarField;
    use crate::math::Math;
    use crate::nexus_spartan::sparse_mlpoly::util::get_bits_canonical_order;
    use ark_bls12_381::Fr;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::FieldVar;
    use ark_r1cs_std::{R1CSVar, ToBitsGadget};
    use ark_relations::ns;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::{random, thread_rng};

    #[test]
    fn test() {
        let cs = ConstraintSystem::<ScalarField>::new_ref();
        let n = u16::rand(&mut thread_rng());
        let fp = FpVar::<ScalarField>::new_variable(ns!(cs, "num"),
                                                    || Ok(ScalarField::from(n)),
                                                    AllocationMode::Witness).unwrap();
        let new_vector: Vec<bool> = get_bits_canonical_order(&fp, 32)
            .iter()
            .map(|element| element.value().unwrap())
            .collect();
        assert_eq!(new_vector, (n as usize).get_bits_canonical_order(32));
    }
}
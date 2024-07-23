/// Poseidon config stolen from sonobe

use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge, poseidon::{PoseidonConfig}};
use ark_crypto_primitives::sponge::poseidon::{find_poseidon_ark_and_mds, PoseidonSponge};
use ark_ff::PrimeField;


/// This Poseidon configuration generator agrees with Circom's Poseidon(4) in the case of BN254's scalar field
pub(crate) fn poseidon_canonical_config<F: PrimeField>() -> PoseidonConfig<F> {
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

fn hash<F: PrimeField + Absorb>(elements: Vec<F>) -> F {
    let poseidon_params: PoseidonConfig<F> = poseidon_canonical_config();
    // apply the hash function
    let mut sponge = PoseidonSponge::new(&poseidon_params);
    for field_element in elements {
        sponge.absorb(&field_element);
    }
    let squeezed_field_element: Vec<F> = sponge.squeeze_field_elements(1);
    squeezed_field_element[0]
}

#[cfg(test)]
mod tests {
    use std::ops::Mul;
    use ark_bn254::{Bn254, Fq, Fr, G1Projective, G2Projective};
    use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
    use ark_crypto_primitives::sponge::poseidon::{PoseidonConfig, PoseidonSponge};
    use ark_ff::Field;
    use ark_std::UniformRand;
    use rand::rngs::OsRng;
    use crate::accumulation::poseidon::{hash, poseidon_canonical_config};

    type FirstCurve = Fr;
    type SecondCurve = Fq;

    #[test]
    fn lagrange_test() {
        let field_vector: Vec<FirstCurve> = vec![FirstCurve::ONE, FirstCurve::ONE];
        let h = hash(field_vector);
        println!("{}", h);
    }
}



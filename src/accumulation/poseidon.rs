/// Poseidon config stolen from sonobe

use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge, poseidon::{PoseidonConfig}};
use ark_crypto_primitives::sponge::poseidon::{find_poseidon_ark_and_mds, PoseidonSponge};
use ark_ff::PrimeField;

pub struct PoseidonHash<F: Absorb + PrimeField> {
    poseidon_params: PoseidonConfig<F>,
    sponge: PoseidonSponge<F>,
}

pub trait PoseidonHashTrait<F: Absorb + PrimeField> {
    fn new() -> Self;

    fn update_sponge<A: Absorb>(&mut self, field_vector: Vec<A>) -> ();

    fn output(&mut self) -> F;
}

impl<F: Absorb + PrimeField> PoseidonHashTrait<F> for PoseidonHash<F> {

    /// This Poseidon configuration generator agrees with Circom's Poseidon(4) in the case of BN254's scalar field
    fn new() -> Self {
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
        let poseidon_params = PoseidonConfig::new(
            full_rounds as usize,
            partial_rounds as usize,
            alpha,
            mds,
            ark,
            rate,
            1,
        );

        Self {
            poseidon_params: poseidon_params.clone(),
            sponge: PoseidonSponge::new(&poseidon_params),
        }
    }

    fn update_sponge<A: Absorb>(&mut self, field_vector: Vec<A>) -> () {
        for field_element in field_vector {
            self.sponge.absorb(&field_element);
        }
    }

    fn output(&mut self) -> F {
        let squeezed_field_element: Vec<F> = self.sponge.squeeze_field_elements(1);
        squeezed_field_element[0]
    }
}


#[cfg(test)]
mod tests {
    use std::ops::Mul;
    use ark_bn254::{Bn254, Fq, Fr, G1Projective, G2Projective};
    use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
    use ark_crypto_primitives::sponge::poseidon::{PoseidonConfig, PoseidonSponge};
    use ark_ec::CurveGroup;
    use ark_ec::short_weierstrass::{SWCurveConfig, Projective};
    use ark_ff::{Field, PrimeField};
    use ark_std::UniformRand;
    use rand::rngs::OsRng;
    use rand::thread_rng;
    use crate::accumulation::poseidon::{PoseidonHash, PoseidonHashTrait};

    type FirstCurve = Fr;
    type SecondCurve = Fq;

    #[test]
    fn lagrange_test() {
        let mut hash_object: PoseidonHash<Fr> = PoseidonHash::new();
        let field_vector: Vec<FirstCurve> = vec![FirstCurve::ONE, FirstCurve::ONE];
        hash_object.update_sponge(field_vector);
        let field_vector: Vec<SecondCurve> = vec![SecondCurve::ONE, SecondCurve::ONE];
        hash_object.update_sponge(field_vector);
        println!("{}", hash_object.output());
    }
}



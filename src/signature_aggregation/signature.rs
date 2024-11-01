use std::ops::Mul;
use ark_ec::AffineRepr;
use ark_ec::pairing::Pairing;
use ark_std::UniformRand;
use rand::Rng;

pub struct BLS;

impl BLS {
    // Generate a new secret key and public key pair
    pub fn generate_key_pair<R: Rng, E: Pairing>(rng: &mut R) -> (E::ScalarField, E::G1Affine) {
        let sk = E::ScalarField::rand(rng);
        let pk = E::G1Affine::generator().mul(sk);
        (sk, pk.into())
    }

    // Sign a message using a provided secret key
    pub fn sign<E: Pairing>(sk: E::ScalarField, message: E::G2Affine) -> E::G2Affine {
        message.mul(&sk).into()
    }

    // Verify the signature against a provided public key and message
    pub fn verify<E: Pairing>(pk: E::G1Affine, sig: E::G2Affine, message: E::G2Affine) {
        let e1 = E::pairing(pk, message);
        let e2 = E::pairing(E::G1Affine::generator(), sig);
        assert_eq!(e1, e2, "invalid signature")
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::AffineRepr;
    use ark_std::{test_rng, UniformRand};
    use rand::thread_rng;
    use crate::constant_for_curves::{G2Affine, E};
    use crate::signature_aggregation::signature::BLS;

    #[test]
    fn test_bls_signature() {
        let mut rng = test_rng();

        // Generate a key pair
        let (sk, pk) = BLS::generate_key_pair::<_, E>(&mut rng);

        // Define a message as an example G2Affine point
        let message = G2Affine::rand(&mut thread_rng());

        // Sign the message
        let signature = BLS::sign::<E>(sk, message);

        // Verify the signature
        BLS::verify::<E>(pk, signature, message);
    }
}

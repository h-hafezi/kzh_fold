/// This is basically KZG with a modified Eval function

use crate::kzg::{KZG10, KZGPowers, KZGUniversalParams, KZGVerifierKey, trim};

/// Prover for PCS:
/// Step 1) Run the HPI protocol
/// Step 2) Run a KZG eval on the log verification scalars

/// Verifier for PCS:
/// Step 1) Verify the HPI protocol
/// Step 2) Compute C'
/// Step 3) Verify KZG eval using C'

#[cfg(test)]
mod tests {
    use ark_bn254::{Bn254, Fr};
    use ark_ec::pairing::Pairing;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    type F = Fr;
    type E = Bn254;
    type Poly = DensePolynomial<<E as Pairing>::ScalarField>;

    #[test]
    pub fn halo_infinite_pcs() {
        // Set up public parameters
        let rng = &mut test_rng();
        let degree = 128 * 128;
        let params = KZG10::<E, Poly>::setup(degree, false, rng).expect("Setup failed");
        let (ck, vk) = trim(&params, degree);

        // Generate commitment
        let polynomial = Poly::rand(degree, rng);

        let (comm, r) = KZG10::<E, Poly>::commit(&ck, &polynomial, None, None).expect("Commitment failed");

        // Open commitment to get proof
        let evaluation_point = F::rand(rng);
        let proof = KZG10::<E, Poly>::open(&ck, &polynomial, evaluation_point, &r).expect("Proof generation failed");

        // Verify proof
        let evaluation_result = polynomial.evaluate(&evaluation_point);
        let is_valid = KZG10::<E, Poly>::check(&vk, &comm, evaluation_point, evaluation_result, &proof).expect("Verification failed");

        assert!(is_valid, "Proof verification failed");
    }
}

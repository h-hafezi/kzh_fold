#[cfg(test)]
mod tests {
    use ark_crypto_primitives::sponge::Absorb;
    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use crate::transcript::transcript::Transcript;
    use crate::transcript::transcript_var::TranscriptVar;
    use ark_ff::{PrimeField, UniformRand};
    use ark_r1cs_std::prelude::*;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
    use rand::{thread_rng, Rng};
    use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit::SumcheckCircuit;
    use crate::nexus_spartan::sumcheck_circuit::sumcheck_circuit_var::SumcheckCircuitVar;
    use crate::nexus_spartan::unipoly::unipoly::CompressedUniPoly;

    // Mock implementation of CompressedUniPoly for test purposes
    impl<F: PrimeField + Absorb> CompressedUniPoly<F> {
        pub fn random<R: Rng>(rng: &mut R, degree_bound: usize) -> Self {
            // Generate a random univariate polynomial with the given degree bound
            let coeffs_except_linear_term: Vec<F> = (0..degree_bound).map(|_| F::rand(rng)).collect();
            CompressedUniPoly { coeffs_except_linear_term }
        }
    }

    type F = ScalarField;

    #[test]
    fn test_sumcheck_circuit() {
        test_sumcheck_circuit_helper::<F>();
    }


    pub fn test_sumcheck_circuit_helper<F: PrimeField + Absorb>() -> SumcheckCircuit<F> {
        // Set up the random number generator
        let mut rng = thread_rng();

        // Parameters for the test
        let num_rounds = 5;
        let degree_bound = 3;

        // Create a claim as a random field element
        let claim = F::rand(&mut rng);

        // Create random compressed univariate polynomials for the Sumcheck proof
        let compressed_polys: Vec<CompressedUniPoly<F>> = (0..num_rounds)
            .map(|_| CompressedUniPoly::random(&mut rng, degree_bound))
            .collect();

        // Initialize the SumcheckCircuit with random values
        let sumcheck_circuit = SumcheckCircuit {
            compressed_polys: compressed_polys.clone(),
            claim: claim.clone(),
            num_rounds,
            degree_bound,
        };

        sumcheck_circuit
    }

    #[test]
    fn test_sumcheck_circuit_var() {
        let label = b"test label";
        let mut transcript = Transcript::<F>::new(label);
        // Initialize the SumcheckCircuit
        let sumcheck_circuit = test_sumcheck_circuit_helper();

        // Run the verification on the normal SumcheckCircuit
        let (e1, r1) = sumcheck_circuit.verify::<E>(&mut transcript);

        // Initialize the constraint system
        let cs: ConstraintSystemRef<F> = ConstraintSystem::new_ref();

        // Initialize a fresh variable transcript for the variable circuit
        let mut transcript_var = TranscriptVar::<F>::new(cs.clone(), label);

        // Initialize the SumCheckCircuitVar with the same values
        let mut sumcheck_circuit_var = SumcheckCircuitVar::new_variable(
            cs.clone(),
            || Ok(sumcheck_circuit.clone()),
            AllocationMode::Witness,
        ).unwrap();

        assert_eq!(sumcheck_circuit, sumcheck_circuit_var.value().unwrap());

        // Run the verification on the SumCheckCircuitVar
        let (e2, r2) = sumcheck_circuit_var.verify(&mut transcript_var);

        // Compare the results between the normal and variable circuits
        assert_eq!(e1, e2.value().unwrap());
        assert_eq!(r1.len(), r2.len());

        // check equality element wise
        for i in 0..r1.len() {
            assert_eq!(r1[i], r2[i].value().unwrap());
        }

        // Ensure that the constraint system is satisfied
        assert!(cs.is_satisfied().unwrap());
        println!("constraint count: {}", cs.num_constraints());
    }
}

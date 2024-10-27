use crate::nexus_spartan::crr1cs::CRR1CSShape;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;

pub fn fold_matrices_evaluations<F: PrimeField + Absorb>(
    shape: &CRR1CSShape<F>,
    r: (Vec<F>, Vec<F>),
    r_prime: (Vec<F>, Vec<F>),
    beta: F,
    evaluations: ((F, F, F), (F, F, F)),
    check_evaluations: bool,
) -> (F, F, F) {
    let (r_x, r_y) = r;
    let (r_x_prime, r_y_prime) = r_prime;
    let (z, z_prime) = evaluations;

    if check_evaluations {
        let z_new = shape.inst.inst.evaluate(&r_x, &r_y);
        assert_eq!(z, z_new);

        let z_prime_new = shape.inst.inst.evaluate(&r_x_prime, &r_y_prime);
        assert_eq!(z_prime, z_prime_new);
    }

    // Perform the random combination for r_x_folded and r_y_folded
    let r_x_folded: Vec<F> = r_x.iter()
        .zip(r_x_prime.iter())
        .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
        .collect();

    let r_y_folded: Vec<F> = r_y.iter()
        .zip(r_y_prime.iter())
        .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
        .collect();

    // Evaluate the folded r_x_folded and r_y_folded
    let new_evaluation = shape.inst.inst.evaluate(&r_x_folded, &r_y_folded);

    // Calculate the proof values
    let proof_A = (new_evaluation.0 - ((F::one() - beta) * z.0 + beta * z_prime.0)) / (beta * (F::one() - beta));
    let proof_B = (new_evaluation.1 - ((F::one() - beta) * z.1 + beta * z_prime.1)) / (beta * (F::one() - beta));
    let proof_C = (new_evaluation.2 - ((F::one() - beta) * z.2 + beta * z_prime.2)) / (beta * (F::one() - beta));

    (proof_A, proof_B, proof_C)
}

#[cfg(test)]
pub mod tests {
    use ark_crypto_primitives::sponge::Absorb;
    use ark_ec::CurveConfig;
    use ark_ec::pairing::Pairing;
    use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
    use ark_ff::{One, PrimeField};
    use crate::constant_for_curves::{ScalarField, E, G1};
    use crate::nexus_spartan::conversion::tests::TrivialCircuit;
    use crate::nexus_spartan::crr1cs::{is_sat, CRR1CSInstance, CRR1CSKey, CRR1CSShape, CRR1CSWitness};
    use crate::nexus_spartan::crr1csproof::CRR1CSProof;
    use crate::nexus_spartan::matrix_evaluation_accumulation::prover::{fold_matrices_evaluations};
    use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
    use crate::pcs::multilinear_pcs::SRS;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use crate::transcript::transcript::Transcript;
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisMode};
    use ark_std::UniformRand;
    use rand::thread_rng;

    type F = ScalarField;

    pub fn matrix_evaluation_setup<F, E, G1>() -> (CRR1CSShape<F>, Vec<F>, Vec<F>)
    where
        F: PrimeField + Absorb,
        E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
        G1: SWCurveConfig<ScalarField=F> + Clone,
        <G1 as CurveConfig>::BaseField: PrimeField
    {
        // Create a new constraint system for: a + b == a^2
        let cs = ConstraintSystem::<F>::new_ref();

        // Example public inputs
        // 4 + 12 == 16 == 4^2
        let a = F::from(4u32);
        let b = F::from(12u32);

        // Instantiate the trivial circuit with inputs
        let circuit = TrivialCircuit { a, b };

        // Generate the constraints
        circuit.generate_constraints(cs.clone()).unwrap();
        assert!(cs.is_satisfied().unwrap());
        println!("{}", cs.num_constraints());

        cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
        cs.finalize();

        // convert to the corresponding Spartan types
        let shape = CRR1CSShape::<F>::convert::<G1>(cs.clone());
        let SRS: SRS<E> = MultilinearPolynomial::setup(4, &mut thread_rng()).unwrap();
        let key: CRR1CSKey<E, MultilinearPolynomial<F>> = CRR1CSKey::new(&SRS, shape.get_num_cons(), shape.get_num_vars());
        let instance: CRR1CSInstance<E, MultilinearPolynomial<F>> = CRR1CSInstance::convert(cs.clone(), &key.keys.ck);
        let witness = CRR1CSWitness::<F>::convert(cs.clone());
        // check that the Spartan instance-witness pair is still satisfying
        assert!(is_sat(&shape, &instance, &witness, &key).unwrap());

        let (num_cons, num_vars, _num_inputs) = (
            shape.get_num_cons(),
            shape.get_num_vars(),
            shape.get_num_inputs(),
        );

        let mut prover_transcript = Transcript::new(b"example");

        let (proof, r_x, r_y) = CRR1CSProof::prove(
            &shape,
            &instance,
            witness,
            &key,
            &mut prover_transcript,
        );

        (shape, r_x, r_y)
    }

    #[test]
    fn test() {
        // ************************************************ set up the spartan proof ************************************************
        let (shape, r_x, r_y) = matrix_evaluation_setup::<F, E, G1>();

        // ************************************************ setup params for folding prover ************************************************
        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let r_x_prime: Vec<F> = (0..r_x.len()).map(|_| F::rand(&mut rng)).collect();
        let r_y_prime: Vec<F> = (0..r_y.len()).map(|_| F::rand(&mut rng)).collect();

        let z: (F, F, F) = shape.inst.inst.evaluate(&r_x, &r_y);
        let z_prime: (F, F, F) = shape.inst.inst.evaluate(&r_x_prime, &r_y_prime);

        let beta = F::rand(&mut rng);

        let (proof_A, proof_B, proof_C) = fold_matrices_evaluations(
            &shape,
            (r_x.clone(), r_y.clone()),
            (r_x_prime.clone(), r_y_prime.clone()),
            beta,
            (z, z_prime),
            true,
        );

        let expected_eval = |z: (F, F, F), z_prime: (F, F, F), proof: (F, F, F)| -> (F, F, F) {
            (
                (F::one() - beta) * z.0 + beta * z_prime.0 + (F::one() - beta) * beta * proof.0,
                (F::one() - beta) * z.1 + beta * z_prime.1 + (F::one() - beta) * beta * proof.1,
                (F::one() - beta) * z.2 + beta * z_prime.2 + (F::one() - beta) * beta * proof.2,
            )
        };

        // Compute the expected evaluation tuple
        let expected_z_folded = expected_eval(z, z_prime, (proof_A, proof_B, proof_C));


        // Perform the random combination for r_x_folded and r_y_folded
        let r_x_folded: Vec<F> = r_x.iter()
            .zip(r_x_prime.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        let r_y_folded: Vec<F> = r_y.iter()
            .zip(r_y_prime.iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        let z_folded: (F, F, F) = shape.inst.inst.evaluate(&r_x_folded, &r_y_folded);

        assert_eq!(expected_z_folded, z_folded);
    }
}
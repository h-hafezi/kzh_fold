use crate::nexus_spartan::crr1cs::CRR1CSShape;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use crate::transcript::transcript::Transcript;

pub fn fold_matrices_evaluations<F: PrimeField + Absorb>(
    shape: &CRR1CSShape<F>,
    running_input: (Vec<F>, Vec<F>),
    current_input: (Vec<F>, Vec<F>),
    transcript: &mut Transcript<F>,
    running_evaluations: (F, F, F),
    current_evaluations: (F, F, F),
    check_evaluations: bool,
) -> (F, (F, F, F)) {
    let (running_input_x, running_input_y) = running_input.clone();
    let (current_input_x, current_input_y) = current_input.clone();

    // append struct to transcript to generate challenge beta
    transcript.append_scalars(b"matrix evaluations",
                              to_sponge_vector(
                                  &running_input,
                                  &current_input,
                                  running_evaluations,
                                  current_evaluations,
                              ).as_slice(),
    );
    let beta = transcript.challenge_scalar(b"beta");


    if check_evaluations {
        let expected_running_evaluations = shape.inst.inst.evaluate(&running_input_x, &running_input_y);
        assert_eq!(running_evaluations, expected_running_evaluations);

        let expected_current_evaluations = shape.inst.inst.evaluate(&current_input_x, &current_input_y);
        assert_eq!(current_evaluations, expected_current_evaluations);
    }

    // Perform the random combination for r_x_folded and r_y_folded
    let folded_input_x: Vec<F> = running_input_x.iter()
        .zip(current_input_x.iter())
        .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
        .collect();

    let folded_input_y: Vec<F> = running_input_y.iter()
        .zip(current_input_y.iter())
        .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
        .collect();

    // Evaluate the folded r_x_folded and r_y_folded
    let new_evaluation = shape.inst.inst.evaluate(&folded_input_x, &folded_input_y);

    // Calculate the proof values
    let proof_A = (new_evaluation.0 - ((F::one() - beta) * running_evaluations.0 + beta * current_evaluations.0)) / (beta * (F::one() - beta));
    let proof_B = (new_evaluation.1 - ((F::one() - beta) * running_evaluations.1 + beta * current_evaluations.1)) / (beta * (F::one() - beta));
    let proof_C = (new_evaluation.2 - ((F::one() - beta) * running_evaluations.2 + beta * current_evaluations.2)) / (beta * (F::one() - beta));

    (beta, (proof_A, proof_B, proof_C))
}

pub fn to_sponge_vector<T: Clone>(running_input: &(Vec<T>, Vec<T>),
                                  new_input: &(Vec<T>, Vec<T>),
                                  running_evaluation: (T, T, T),
                                  new_evaluation: (T, T, T),
) -> Vec<T> {
    let mut res = Vec::new();
    res.extend(running_input.0.clone());
    res.extend(running_input.1.clone());
    res.extend(new_input.0.clone());
    res.extend(new_input.1.clone());
    res.extend(vec![new_evaluation.0, new_evaluation.1, new_evaluation.2]);
    res.extend(vec![running_evaluation.0, running_evaluation.1, running_evaluation.2]);

    res
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
        let (shape, running_input_x, running_input_y) = matrix_evaluation_setup::<F, E, G1>();

        // ************************************************ setup params for folding prover ************************************************
        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let current_input_x: Vec<F> = (0..running_input_x.len()).map(|_| F::rand(&mut rng)).collect();
        let current_input_y: Vec<F> = (0..running_input_y.len()).map(|_| F::rand(&mut rng)).collect();

        let running_evaluations: (F, F, F) = shape.inst.inst.evaluate(&running_input_x, &running_input_y);
        let current_evaluations: (F, F, F) = shape.inst.inst.evaluate(&current_input_x, &current_input_y);

        let (beta, (proof_A, proof_B, proof_C)) = fold_matrices_evaluations(
            &shape,
            (running_input_x.clone(), running_input_y.clone()),
            (current_input_x.clone(), current_input_y.clone()),
            &mut Transcript::new(b"new transcript"),
            running_evaluations,
            current_evaluations,
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
        let expected_folded_evaluations = expected_eval(running_evaluations, current_evaluations, (proof_A, proof_B, proof_C));


        // Perform the random combination for r_x_folded and r_y_folded
        let folded_input_x: Vec<F> = running_input_x.iter()
            .zip(current_input_x.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        let folded_input_y: Vec<F> = running_input_y.iter()
            .zip(current_input_y.iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        let folded_evaluations: (F, F, F) = shape.inst.inst.evaluate(&folded_input_x, &folded_input_y);

        assert_eq!(expected_folded_evaluations, folded_evaluations);
    }
}
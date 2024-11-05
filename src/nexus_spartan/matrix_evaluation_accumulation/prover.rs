use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use rand::RngCore;
use rayon::prelude::IntoParallelRefIterator;
use crate::math::Math;
use crate::nexus_spartan::crr1cs::CRR1CSShape;
use crate::polynomial::univariate::univariate::PolynomialInterpolator;
use crate::transcript::transcript::Transcript;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixEvaluationAccumulator<F: PrimeField + Absorb> {
    // (r_x, r_y)
    pub evaluation_point: (Vec<F>, Vec<F>),
    // A(r_x, r_y), B(r_x, r_y), C(r_x, r_y)
    pub evaluations: (F, F, F),
}

// todo: I think it's weird, it should take a shape not generate everything randomly
impl<F: PrimeField + Absorb> MatrixEvaluationAccumulator<F> {
    pub fn rand<R: RngCore>(x_len: usize, y_len: usize, rng: &mut R) -> Self {
        let eval_point_x: Vec<F> = (0..x_len).map(|_| F::rand(rng)).collect();
        let eval_point_y: Vec<F> = (0..y_len).map(|_| F::rand(rng)).collect();
        let evaluations = (F::rand(rng), F::rand(rng), F::rand(rng));

        Self {
            evaluation_point: (eval_point_x, eval_point_y),
            evaluations,
        }
    }
}

pub fn fold_matrices_evaluations<F: PrimeField + Absorb>(
    shape: &CRR1CSShape<F>,
    eval_point_1: (Vec<F>, Vec<F>),
    eval_point_2: (Vec<F>, Vec<F>),
    transcript: &mut Transcript<F>,
    evals_1: (F, F, F),
    evals_2: (F, F, F),
    check_evaluations: bool,
) -> (F, (PolynomialInterpolator<F>, PolynomialInterpolator<F>, PolynomialInterpolator<F>)) {
    let (eval_point_1_x, eval_point_1_y) = eval_point_1.clone();
    let (eval_point_2_x, eval_point_2_y) = eval_point_2.clone();

    let (q_A, q_B, q_C) = compute_q(
        &shape,
        eval_point_1.clone(),
        eval_point_2.clone(),
    );

    // append struct to transcript to generate challenge beta
    transcript.append_scalars(b"matrix evaluations",
                              to_sponge_vector(
                                  &eval_point_1,
                                  &eval_point_2,
                                  evals_1,
                                  evals_2,
                              ).as_slice(),
    );
    transcript.append_scalars(b"quotient polynomial", q_A.coefficients.as_slice());
    transcript.append_scalars(b"quotient polynomial", q_B.coefficients.as_slice());
    transcript.append_scalars(b"quotient polynomial", q_C.coefficients.as_slice());
    let beta = transcript.challenge_scalar(b"beta");

    // also check evaluation points
    if check_evaluations {
        let expected_evals_1 = shape.inst.inst.evaluate(&eval_point_1_x, &eval_point_1_y);
        assert_eq!(evals_1, expected_evals_1);

        let expected_evals_2 = shape.inst.inst.evaluate(&eval_point_2_x, &eval_point_2_y);
        assert_eq!(evals_2, expected_evals_2);
    }

    (beta, (q_A, q_B, q_C))
}

pub fn compute_q<F: PrimeField + Absorb>(shape: &CRR1CSShape<F>,
                                         eval_point_1: (Vec<F>, Vec<F>),
                                         eval_point_2: (Vec<F>, Vec<F>),
) -> (PolynomialInterpolator<F>, PolynomialInterpolator<F>, PolynomialInterpolator<F>) {
    let (eval_point_1_x, eval_point_1_y) = eval_point_1.clone();
    let (eval_point_2_x, eval_point_2_y) = eval_point_2.clone();

    assert_eq!(eval_point_1_x.len(), eval_point_2_x.len(), "arrays with unequal length");
    assert_eq!(eval_point_1_y.len(), eval_point_2_y.len(), "arrays with unequal length");

    let evals_1 = shape.inst.inst.evaluate(&eval_point_1_x, &eval_point_1_y);
    let evals_2 = shape.inst.inst.evaluate(&eval_point_2_x, &eval_point_2_y);

    let mut Q_A = Vec::new();
    let mut Q_B = Vec::new();
    let mut Q_C = Vec::new();

    // we start from 2 because the term (1 - x) * x is zero at 0 and 1
    for i in 2..(64 *  shape.get_num_cons()).log_2() {
        let beta = F::from(i as u128);

        // Perform the random combination for r_x_folded and r_y_folded
        let folded_input_x: Vec<F> = eval_point_1_x.par_iter()
            .zip(eval_point_2_x.par_iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        let folded_input_y: Vec<F> = eval_point_1_y.par_iter()
            .zip(eval_point_2_y.par_iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        // Evaluate the folded r_x_folded and r_y_folded
        let new_evaluation = shape.inst.inst.evaluate(&folded_input_x, &folded_input_y);

        // Calculate the proof values
        let q_A = (new_evaluation.0 - ((F::one() - beta) * evals_1.0 + beta * evals_2.0)) / (beta * (F::one() - beta));
        let q_B = (new_evaluation.1 - ((F::one() - beta) * evals_1.1 + beta * evals_2.1)) / (beta * (F::one() - beta));
        let q_C = (new_evaluation.2 - ((F::one() - beta) * evals_1.2 + beta * evals_2.2)) / (beta * (F::one() - beta));

        Q_A.push((beta, q_A));
        Q_B.push((beta, q_B));
        Q_C.push((beta, q_C));
    }

    let mut poly_Q_A = PolynomialInterpolator::new();
    poly_Q_A.interpolate(Q_A.as_slice());

    let mut poly_Q_B = PolynomialInterpolator::new();
    poly_Q_B.interpolate(Q_B.as_slice());

    let mut poly_Q_C = PolynomialInterpolator::new();
    poly_Q_C.interpolate(Q_C.as_slice());

    (poly_Q_A, poly_Q_B, poly_Q_C)
}

pub fn to_sponge_vector<T: Clone>(eval_point_1: &(Vec<T>, Vec<T>),
                                  new_input: &(Vec<T>, Vec<T>),
                                  running_evaluation: (T, T, T),
                                  new_evaluation: (T, T, T),
) -> Vec<T> {
    let mut res = Vec::new();
    res.extend(eval_point_1.0.clone());
    res.extend(eval_point_1.1.clone());
    res.extend(new_input.0.clone());
    res.extend(new_input.1.clone());
    res.extend(vec![new_evaluation.0, new_evaluation.1, new_evaluation.2]);
    res.extend(vec![running_evaluation.0, running_evaluation.1, running_evaluation.2]);

    res
}

#[cfg(test)]
pub mod tests {
    use crate::constant_for_curves::{ScalarField, E, G1};
    use crate::nexus_spartan::conversion::tests::TrivialCircuit;
    use crate::nexus_spartan::crr1cs::{is_sat, CRR1CSInstance, CRR1CSKey, CRR1CSShape, CRR1CSWitness};
    use crate::nexus_spartan::crr1csproof::CRR1CSProof;
    use crate::nexus_spartan::matrix_evaluation_accumulation::prover::fold_matrices_evaluations;
    use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
    use crate::pcs::multilinear_pcs::PolynomialCommitmentSRS;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use crate::transcript::transcript::Transcript;
    use ark_crypto_primitives::sponge::Absorb;
    use ark_ec::pairing::Pairing;
    use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
    use ark_ec::CurveConfig;
    use ark_ff::{One, PrimeField};
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisMode};
    use ark_std::UniformRand;
    use rand::thread_rng;

    type F = ScalarField;

    pub fn matrix_evaluation_setup<F, E, G1>() -> (CRR1CSShape<F>, Vec<F>, Vec<F>)
    where
        F: PrimeField + Absorb,
        E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
        G1: SWCurveConfig<ScalarField=F> + Clone,
        <G1 as CurveConfig>::BaseField: PrimeField,
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
        let SRS: PolynomialCommitmentSRS<E> = MultilinearPolynomial::setup(4, &mut thread_rng()).unwrap();
        let key: CRR1CSKey<E, MultilinearPolynomial<F>> = CRR1CSKey::new(&SRS, shape.get_num_cons(), shape.get_num_vars());
        let instance: CRR1CSInstance<E, MultilinearPolynomial<F>> = CRR1CSInstance::convert(cs.clone(), &key.keys.ck);
        let witness = CRR1CSWitness::<F>::convert(cs.clone());
        // check that the Spartan instance-witness pair is still satisfying
        assert!(is_sat(&shape, &instance, &witness, &key).unwrap());

        let mut prover_transcript = Transcript::new(b"example");

        let (_, r_x, r_y) = CRR1CSProof::prove(
            &shape,
            &instance,
            witness,
            &key,
            &mut prover_transcript,
        );

        (shape, r_x, r_y)
    }

    #[test]
    fn test_matrix_evaluation_prover() {
        // ************************************************ set up the spartan proof ************************************************
        let (shape, eval_point_1_x, eval_point_1_y) = matrix_evaluation_setup::<F, E, G1>();

        // ************************************************ setup params for folding prover ************************************************
        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let eval_point_2_x: Vec<F> = (0..eval_point_1_x.len()).map(|_| F::rand(&mut rng)).collect();
        let eval_point_2_y: Vec<F> = (0..eval_point_1_y.len()).map(|_| F::rand(&mut rng)).collect();

        let evals_1: (F, F, F) = shape.inst.inst.evaluate(&eval_point_1_x, &eval_point_1_y);
        let evals_2: (F, F, F) = shape.inst.inst.evaluate(&eval_point_2_x, &eval_point_2_y);

        let (beta, (proof_A, proof_B, proof_C)) = fold_matrices_evaluations(
            &shape,
            (eval_point_1_x.clone(), eval_point_1_y.clone()),
            (eval_point_2_x.clone(), eval_point_2_y.clone()),
            &mut Transcript::new(b"new transcript"),
            evals_1,
            evals_2,
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
        let expected_folded_evaluations = expected_eval(
            evals_1,
            evals_2,
            (
                proof_A.evaluate(beta),
                proof_B.evaluate(beta),
                proof_C.evaluate(beta),
            ),
        );


        // Perform the random combination for r_x_folded and r_y_folded
        let folded_input_x: Vec<F> = eval_point_1_x.iter()
            .zip(eval_point_2_x.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        let folded_input_y: Vec<F> = eval_point_1_y.iter()
            .zip(eval_point_2_y.iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        let folded_evaluations: (F, F, F) = shape.inst.inst.evaluate(&folded_input_x, &folded_input_y);

        assert_eq!(expected_folded_evaluations, folded_evaluations);
    }
}

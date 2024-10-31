use crate::nexus_spartan::crr1cs::CRR1CSShape;
use crate::nexus_spartan::matrix_evaluation_accumulation::prover::{fold_matrices_evaluations, to_sponge_vector};
use crate::transcript::transcript::Transcript;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use rand::Rng;
use std::borrow::Borrow;
use std::ops::Mul;

/// XXX Hossein this should contain two `MatrixEvaluationAccumulator` instead of the raw data
/// Verify the accumulation of (r_x, r_y, z) and (r_x', r_y', z') into (r_x'', r_y'', z'')
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixEvaluationAccVerifier<F: PrimeField + Absorb> {
    // (r_x, r_y)
    pub eval_point_1: (Vec<F>, Vec<F>),
    // z = A(r_x, r_y)
    pub evals_1: (F, F, F),

    // (r_x', r_y')
    pub eval_point_2: (Vec<F>, Vec<F>),
    // z' = A(r_x', r_y')
    pub evals_2: (F, F, F),

    // XXX Hossein This should likely be a polynomial q(x). Not its evaluations
    pub proof: (F, F, F),
}

impl<F: PrimeField + Absorb> MatrixEvaluationAccVerifier<F> {
    pub fn accumulate(&self, transcript: &mut Transcript<F>) -> ((Vec<F>, Vec<F>), (F, F, F)) {
        // add the whole struct to transcript
        transcript.append_scalars(
            b"whole struct",
            to_sponge_vector(
                &self.eval_point_1,
                &self.eval_point_2,
                self.evals_1,
                self.evals_2)
                .as_slice(),
        );
        let beta = transcript.challenge_scalar(b"beta");

        // parse variables
        let (eval_point_1_x, eval_point_1_y) = self.eval_point_1.clone();
        let (eval_point_2_x, eval_point_2_y) = self.eval_point_2.clone();

        // Perform the random combination for r_x_folded and r_y_folded
        let folded_input_x: Vec<F> = eval_point_1_x.iter()
            .zip(eval_point_2_x.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        let folded_input_y: Vec<F> = eval_point_1_y.iter()
            .zip(eval_point_2_y.iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        let expected_eval = |z: (F, F, F), z_prime: (F, F, F), proof: (F, F, F)| -> (F, F, F) {
            (
                (F::one() - beta) * z.0 + beta * z_prime.0 + (F::one() - beta) * beta * proof.0,
                (F::one() - beta) * z.1 + beta * z_prime.1 + (F::one() - beta) * beta * proof.1,
                (F::one() - beta) * z.2 + beta * z_prime.2 + (F::one() - beta) * beta * proof.2,
            )
        };

        // Compute the expected evaluation tuple
        let next_evaluation = expected_eval(self.evals_1, self.evals_2, self.proof);

        ((folded_input_x, folded_input_y), next_evaluation)
    }

    pub fn random_from_eval_point<R: Rng>(shape: &CRR1CSShape<F>, rx: Vec<F>, ry: Vec<F>, mut transcript: Transcript<F>, mut rng: R) -> Self {
        // Generate random elements in the field for r_x_prime and r_y_prime
        let r_x_prime: Vec<F> = (0..rx.len()).map(|_| F::rand(&mut rng)).collect();
        let r_y_prime: Vec<F> = (0..ry.len()).map(|_| F::rand(&mut rng)).collect();

        // Generate A(r_x', r_y'), B(r_x', r_y'), C(r_x', r_y')
        let z_prime: (F, F, F) = shape.inst.inst.evaluate(&r_x_prime, &r_y_prime);

        let current_A_B_C_evaluations = shape.inst.inst.evaluate(&rx, &ry);

        let (_, matrix_evaluation_proof) = fold_matrices_evaluations(
            &shape,
            (rx.clone(), ry.clone()),
            (r_x_prime.clone(), r_y_prime.clone()),
            &mut transcript,
            current_A_B_C_evaluations,
            z_prime,
            true,
        );

        let verifier = MatrixEvaluationAccVerifier {
            eval_point_1: (rx.clone(), ry.clone()),
            eval_point_2: (r_x_prime.clone(), r_y_prime.clone()),
            evals_1: current_A_B_C_evaluations,
            evals_2: z_prime,
            proof: matrix_evaluation_proof,
        };

        verifier
    }
}

pub struct MatrixEvaluationAccVerifierVar<F: PrimeField + Absorb> {
    pub eval_point_1: (Vec<FpVar<F>>, Vec<FpVar<F>>),
    pub eval_point_2: (Vec<FpVar<F>>, Vec<FpVar<F>>),

    pub evals_1: (FpVar<F>, FpVar<F>, FpVar<F>),
    pub evals_2: (FpVar<F>, FpVar<F>, FpVar<F>),

    pub proof: (FpVar<F>, FpVar<F>, FpVar<F>),
}

// Implement the `AllocVar` and `R1CSVar` traits to integrate with R1CS constraint systems
impl<F: PrimeField + Absorb> AllocVar<MatrixEvaluationAccVerifier<F>, F> for MatrixEvaluationAccVerifierVar<F> {
    fn new_variable<T: Borrow<MatrixEvaluationAccVerifier<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let binding = f()?;
        let instance = binding.borrow();

        let eval_point_1_x = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.eval_point_1.0.clone()), mode)?;
        let eval_point_1_y = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.eval_point_1.1.clone()), mode)?;
        let eval_point_2_x = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.eval_point_2.0.clone()), mode)?;
        let eval_point_2_y = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.eval_point_2.1.clone()), mode)?;

        let running_eval_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.evals_1.0), mode)?;
        let running_eval_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.evals_1.1), mode)?;
        let running_eval_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.evals_1.2), mode)?;

        let current_eval_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.evals_2.0), mode)?;
        let current_eval_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.evals_2.1), mode)?;
        let current_eval_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.evals_2.2), mode)?;

        let proof_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.0), mode)?;
        let proof_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.1), mode)?;
        let proof_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.2), mode)?;

        Ok(MatrixEvaluationAccVerifierVar {
            eval_point_1: (eval_point_1_x, eval_point_1_y),
            eval_point_2: (eval_point_2_x, eval_point_2_y),
            evals_1: (running_eval_0, running_eval_1, running_eval_2),
            evals_2: (current_eval_0, current_eval_1, current_eval_2),
            proof: (proof_0, proof_1, proof_2),
        })
    }
}

impl<F: PrimeField + Absorb> R1CSVar<F> for MatrixEvaluationAccVerifierVar<F> {
    type Value = MatrixEvaluationAccVerifier<F>;

    /// Returns a reference to the constraint system.
    fn cs(&self) -> ConstraintSystemRef<F> {
        unreachable!()
    }

    /// Retrieves the underlying values from the variable if they exist.
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        Ok(MatrixEvaluationAccVerifier {
            eval_point_1: (
                self.eval_point_1.0.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?,
                self.eval_point_1.1.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?
            ),
            eval_point_2: (
                self.eval_point_2.0.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?,
                self.eval_point_2.1.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?
            ),
            evals_1: (
                self.evals_1.0.value()?,
                self.evals_1.1.value()?,
                self.evals_1.2.value()?,
            ),
            evals_2: (
                self.evals_2.0.value()?,
                self.evals_2.1.value()?,
                self.evals_2.2.value()?,
            ),
            proof: (
                self.proof.0.value()?,
                self.proof.1.value()?,
                self.proof.2.value()?,
            ),
        })
    }
}


impl<F: PrimeField + Absorb> MatrixEvaluationAccVerifierVar<F> {
    /// Constraint-friendly accumulate function.
    pub fn accumulate(&self, transcript: &mut TranscriptVar<F>) -> ((Vec<FpVar<F>>, Vec<FpVar<F>>), (FpVar<F>, FpVar<F>, FpVar<F>)) {
        // add the whole struct to transcript
        transcript.append_scalars(
            b"whole struct",
            to_sponge_vector(
                &self.eval_point_1,
                &self.eval_point_2,
                self.evals_1.clone(),
                self.evals_2.clone(),
            ).as_slice());
        let beta = transcript.challenge_scalar(b"beta");

        // Parse variables
        let (eval_point_1_x, eval_point_1_y) = &self.eval_point_1;
        let (eval_point_2_x, eval_point_2_y) = &self.eval_point_2;

        // Compute r_x_folded and r_y_folded
        let folded_input_x: Vec<FpVar<F>> = eval_point_1_x.iter()
            .zip(eval_point_2_x.iter())
            .map(|(rx, rx_prime)| {
                rx.mul(&(FpVar::one() - &beta)) + rx_prime.mul(&beta)
            })
            .collect::<Vec<_>>();

        let folded_input_y: Vec<FpVar<F>> = eval_point_1_y.iter()
            .zip(eval_point_2_y.iter())
            .map(|(ry, ry_prime)| {
                ry.mul(&(FpVar::one() - &beta)) + ry_prime.mul(&beta)
            })
            .collect::<Vec<_>>();

        // Compute expected evaluation
        let expected_eval = |z: &(FpVar<F>, FpVar<F>, FpVar<F>),
                             z_prime: &(FpVar<F>, FpVar<F>, FpVar<F>),
                             proof: &(FpVar<F>, FpVar<F>, FpVar<F>)| -> (FpVar<F>, FpVar<F>, FpVar<F>) {
            (
                z.0.clone() * (FpVar::one() - beta.clone()) + z_prime.0.clone() * beta.clone() + proof.0.clone() * (FpVar::one() - beta.clone()) * beta.clone(),
                z.1.clone() * (FpVar::one() - beta.clone()) + z_prime.1.clone() * beta.clone() + proof.1.clone() * (FpVar::one() - beta.clone()) * beta.clone(),
                z.2.clone() * (FpVar::one() - beta.clone()) + z_prime.2.clone() * beta.clone() + proof.2.clone() * (FpVar::one() - beta.clone()) * beta.clone(),
            )
        };

        let next_evaluation = expected_eval(
            &self.evals_1,
            &self.evals_2,
            &self.proof,
        );

        ((folded_input_x, folded_input_y), next_evaluation)
    }
}

#[cfg(test)]
pub mod tests {
    use crate::constant_for_curves::{ScalarField, E, G1};
    use crate::nexus_spartan::matrix_evaluation_accumulation::prover::fold_matrices_evaluations;
    use crate::nexus_spartan::matrix_evaluation_accumulation::prover::tests::matrix_evaluation_setup;
    use crate::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationAccVerifier, MatrixEvaluationAccVerifierVar};
    use crate::transcript::transcript::Transcript;
    use crate::transcript::transcript_var::TranscriptVar;
    use ark_ff::One;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;

    type F = ScalarField;

    #[test]
    fn test_A_B_C_eval_non_zk() {
        let (shape, eval_point_1_x, eval_point_1_y) = matrix_evaluation_setup::<F, E, G1>();

        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let eval_point_2_x: Vec<F> = (0..eval_point_1_x.len()).map(|_| F::rand(&mut rng)).collect();
        let eval_point_2_y: Vec<F> = (0..eval_point_1_y.len()).map(|_| F::rand(&mut rng)).collect();

        let evals_1: (F, F, F) = shape.inst.inst.evaluate(&eval_point_1_x, &eval_point_1_y);
        let evals_2: (F, F, F) = shape.inst.inst.evaluate(&eval_point_2_x, &eval_point_2_y);

        let mut prover_transcript = Transcript::<F>::new(b"new transcript");
        let mut verifier_transcript = prover_transcript.clone();

        let (beta, matrix_evaluation_proof) = fold_matrices_evaluations(
            &shape,
            (eval_point_1_x.clone(), eval_point_1_y.clone()),
            (eval_point_2_x.clone(), eval_point_2_y.clone()),
            &mut prover_transcript,
            evals_1, evals_2,
            true,
        );

        let verifier = MatrixEvaluationAccVerifier {
            eval_point_1: (eval_point_1_x.clone(), eval_point_1_y.clone()),
            eval_point_2: (eval_point_2_x.clone(), eval_point_2_y.clone()),
            evals_1: evals_1,
            evals_2: evals_2,
            proof: matrix_evaluation_proof,
        };

        let ((folded_input_x, folded_input_y), folded_evaluations) = verifier.accumulate(&mut verifier_transcript);

        // Perform the random combination for r_x_folded and r_y_folded
        let folded_input_x_expected: Vec<F> = eval_point_1_x.iter()
            .zip(eval_point_2_x.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        assert_eq!(folded_input_x_expected, folded_input_x);

        let folded_input_y_expected: Vec<F> = eval_point_1_y.iter()
            .zip(eval_point_2_y.iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        assert_eq!(folded_input_y_expected, folded_input_y);

        let new_evaluations_expected: (F, F, F) = shape.inst.inst.evaluate(&folded_input_x, &folded_input_y);

        assert_eq!(new_evaluations_expected, folded_evaluations);
    }

    #[test]
    fn test_A_B_C_eval_zk() {
        let (shape, eval_point_1_x, eval_point_1_y) = matrix_evaluation_setup::<F, E, G1>();
        let cs = ConstraintSystem::<F>::new_ref();

        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let eval_point_2_x: Vec<F> = (0..eval_point_1_x.len()).map(|_| F::rand(&mut rng)).collect();
        let eval_point_2_y: Vec<F> = (0..eval_point_1_y.len()).map(|_| F::rand(&mut rng)).collect();


        let evals_1: (F, F, F) = shape.inst.inst.evaluate(&eval_point_1_x, &eval_point_1_y);
        let evals_2: (F, F, F) = shape.inst.inst.evaluate(&eval_point_2_x, &eval_point_2_y);

        let mut prover_transcript = Transcript::<F>::new(b"new transcript");
        let mut verifier_transcript = prover_transcript.clone();
        let mut verifier_transcript_var = TranscriptVar::from_transcript(
            cs.clone(),
            verifier_transcript.clone(),
        );

        let (_, matrix_evaluation_proof) = fold_matrices_evaluations(
            &shape,
            (eval_point_1_x.clone(), eval_point_1_y.clone()),
            (eval_point_2_x.clone(), eval_point_2_y.clone()),
            &mut prover_transcript,
            evals_1,
            evals_2,
            true,
        );

        let verifier = MatrixEvaluationAccVerifier {
            eval_point_1: (eval_point_1_x.clone(), eval_point_1_y.clone()),
            eval_point_2: (eval_point_2_x.clone(), eval_point_2_y.clone()),
            evals_1: evals_1,
            evals_2: evals_2,
            proof: matrix_evaluation_proof,
        };

        let verifier_var = MatrixEvaluationAccVerifierVar::new_variable(
            cs.clone(),
            || Ok(verifier.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // assert value function works correctly
        assert_eq!(verifier, verifier_var.value().unwrap());

        // call accumulate and see if equal to the non-zk version
        let ((folded_input_x, folded_input_y), folded_evaluations) = verifier.accumulate(&mut verifier_transcript);
        let ((folded_input_x_var, folded_input_y_var), folded_evaluations_var) = verifier_var.accumulate(&mut verifier_transcript_var);

        // check cs status
        assert!(cs.is_satisfied().unwrap());
        println!("constraint num: {}", cs.num_constraints());

        // Assert that the vectors are equal
        assert_eq!(folded_input_x, folded_input_x_var.iter().map(|val| val.value().unwrap()).collect::<Vec<F>>(), "r_x_folded vectors are not equal");
        assert_eq!(folded_input_y, folded_input_y_var.iter().map(|val| val.value().unwrap()).collect::<Vec<F>>(), "r_y_folded vectors are not equal");
        assert_eq!(folded_evaluations, (
            folded_evaluations_var.0.value().unwrap(),
            folded_evaluations_var.1.value().unwrap(),
            folded_evaluations_var.2.value().unwrap(),
        ), "new_evaluations tuples are not equal");
    }
}

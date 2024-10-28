use crate::nexus_spartan::matrix_evaluation_accumulation::prover::to_sponge_vector;
use crate::transcript::transcript::Transcript;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;
use std::ops::Mul;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixEvaluationVerifier<F: PrimeField + Absorb> {
    pub running_input: (Vec<F>, Vec<F>),
    pub current_input: (Vec<F>, Vec<F>),
    pub running_evaluations: (F, F, F),
    pub current_evaluations: (F, F, F),
    pub proof: (F, F, F),
}

impl<F: PrimeField + Absorb> MatrixEvaluationVerifier<F> {
    pub fn accumulate(&self, transcript: &mut Transcript<F>) -> ((Vec<F>, Vec<F>), (F, F, F)) {
        // add the whole struct to transcript
        transcript.append_scalars(
            b"whole struct",
            to_sponge_vector(
                &self.running_input,
                &self.current_input,
                self.running_evaluations,
                self.current_evaluations)
                .as_slice(),
        );
        let beta = transcript.challenge_scalar(b"beta");

        // parse variables
        let (running_input_x, running_input_y) = self.running_input.clone();
        let (current_input_x, current_input_y) = self.current_input.clone();

        // Perform the random combination for r_x_folded and r_y_folded
        let folded_input_x: Vec<F> = running_input_x.iter()
            .zip(current_input_x.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        let folded_input_y: Vec<F> = running_input_y.iter()
            .zip(current_input_y.iter())
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
        let next_evaluation = expected_eval(self.running_evaluations, self.current_evaluations, self.proof);

        ((folded_input_x, folded_input_y), next_evaluation)
    }
}

pub struct MatrixEvaluationVerifierVar<F: PrimeField + Absorb> {
    pub running_input: (Vec<FpVar<F>>, Vec<FpVar<F>>),
    pub current_input: (Vec<FpVar<F>>, Vec<FpVar<F>>),
    pub running_evaluations: (FpVar<F>, FpVar<F>, FpVar<F>),
    pub current_evaluations: (FpVar<F>, FpVar<F>, FpVar<F>),
    pub proof: (FpVar<F>, FpVar<F>, FpVar<F>),
}

// Implement the `AllocVar` and `R1CSVar` traits to integrate with R1CS constraint systems
impl<F: PrimeField + Absorb> AllocVar<MatrixEvaluationVerifier<F>, F> for MatrixEvaluationVerifierVar<F> {
    fn new_variable<T: Borrow<MatrixEvaluationVerifier<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let binding = f()?;
        let instance = binding.borrow();

        let running_input_x = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.running_input.0.clone()), mode)?;
        let running_input_y = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.running_input.1.clone()), mode)?;
        let current_input_x = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.current_input.0.clone()), mode)?;
        let current_input_y = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.current_input.1.clone()), mode)?;

        let running_eval_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.running_evaluations.0), mode)?;
        let running_eval_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.running_evaluations.1), mode)?;
        let running_eval_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.running_evaluations.2), mode)?;

        let current_eval_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.current_evaluations.0), mode)?;
        let current_eval_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.current_evaluations.1), mode)?;
        let current_eval_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.current_evaluations.2), mode)?;

        let proof_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.0), mode)?;
        let proof_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.1), mode)?;
        let proof_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.2), mode)?;

        Ok(MatrixEvaluationVerifierVar {
            running_input: (running_input_x, running_input_y),
            current_input: (current_input_x, current_input_y),
            running_evaluations: (running_eval_0, running_eval_1, running_eval_2),
            current_evaluations: (current_eval_0, current_eval_1, current_eval_2),
            proof: (proof_0, proof_1, proof_2),
        })
    }
}

impl<F: PrimeField + Absorb> R1CSVar<F> for MatrixEvaluationVerifierVar<F> {
    type Value = MatrixEvaluationVerifier<F>;

    /// Returns a reference to the constraint system.
    fn cs(&self) -> ConstraintSystemRef<F> {
        unreachable!()
    }

    /// Retrieves the underlying values from the variable if they exist.
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        Ok(MatrixEvaluationVerifier {
            running_input: (
                self.running_input.0.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?,
                self.running_input.1.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?
            ),
            current_input: (
                self.current_input.0.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?,
                self.current_input.1.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?
            ),
            running_evaluations: (
                self.running_evaluations.0.value()?,
                self.running_evaluations.1.value()?,
                self.running_evaluations.2.value()?,
            ),
            current_evaluations: (
                self.current_evaluations.0.value()?,
                self.current_evaluations.1.value()?,
                self.current_evaluations.2.value()?,
            ),
            proof: (
                self.proof.0.value()?,
                self.proof.1.value()?,
                self.proof.2.value()?,
            ),
        })
    }
}


impl<F: PrimeField + Absorb> MatrixEvaluationVerifierVar<F> {
    /// Constraint-friendly accumulate function.
    pub fn accumulate(&self, transcript: &mut TranscriptVar<F>) -> ((Vec<FpVar<F>>, Vec<FpVar<F>>), (FpVar<F>, FpVar<F>, FpVar<F>)) {
        // add the whole struct to transcript
        transcript.append_scalars(
            b"whole struct",
            to_sponge_vector(
                &self.running_input,
                &self.current_input,
                self.running_evaluations.clone(),
                self.current_evaluations.clone(),
            ).as_slice());
        let beta = transcript.challenge_scalar(b"beta");

        // Parse variables
        let (running_input_x, running_input_y) = &self.running_input;
        let (current_input_x, current_input_y) = &self.current_input;

        // Compute r_x_folded and r_y_folded
        let folded_input_x: Vec<FpVar<F>> = running_input_x.iter()
            .zip(current_input_x.iter())
            .map(|(rx, rx_prime)| {
                rx.mul(&(FpVar::one() - &beta)) + rx_prime.mul(&beta)
            })
            .collect::<Vec<_>>();

        let folded_input_y: Vec<FpVar<F>> = running_input_y.iter()
            .zip(current_input_y.iter())
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
            &self.running_evaluations,
            &self.current_evaluations,
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
    use crate::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationVerifier, MatrixEvaluationVerifierVar};
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
    fn test_non_zk() {
        let (shape, running_input_x, running_input_y) = matrix_evaluation_setup::<F, E, G1>();

        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let current_input_x: Vec<F> = (0..running_input_x.len()).map(|_| F::rand(&mut rng)).collect();
        let current_input_y: Vec<F> = (0..running_input_y.len()).map(|_| F::rand(&mut rng)).collect();

        let running_evaluations: (F, F, F) = shape.inst.inst.evaluate(&running_input_x, &running_input_y);
        let current_evaluations: (F, F, F) = shape.inst.inst.evaluate(&current_input_x, &current_input_y);

        let mut prover_transcript = Transcript::<F>::new(b"new transcript");
        let mut verifier_transcript = prover_transcript.clone();

        let (beta, matrix_evaluation_proof) = fold_matrices_evaluations(
            &shape,
            (running_input_x.clone(), running_input_y.clone()),
            (current_input_x.clone(), current_input_y.clone()),
            &mut prover_transcript,
            running_evaluations, current_evaluations,
            true,
        );

        let verifier = MatrixEvaluationVerifier {
            running_input: (running_input_x.clone(), running_input_y.clone()),
            current_input: (current_input_x.clone(), current_input_y.clone()),
            running_evaluations,
            current_evaluations,
            proof: matrix_evaluation_proof,
        };

        let ((folded_input_x, folded_input_y), folded_evaluations) = verifier.accumulate(&mut verifier_transcript);

        // Perform the random combination for r_x_folded and r_y_folded
        let folded_input_x_expected: Vec<F> = running_input_x.iter()
            .zip(current_input_x.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        assert_eq!(folded_input_x_expected, folded_input_x);

        let folded_input_y_expected: Vec<F> = running_input_y.iter()
            .zip(current_input_y.iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        assert_eq!(folded_input_y_expected, folded_input_y);

        let new_evaluations_expected: (F, F, F) = shape.inst.inst.evaluate(&folded_input_x, &folded_input_y);

        assert_eq!(new_evaluations_expected, folded_evaluations);
    }

    #[test]
    fn test_zk() {
        let (shape, running_input_x, running_input_y) = matrix_evaluation_setup::<F, E, G1>();
        let cs = ConstraintSystem::<F>::new_ref();

        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let current_input_x: Vec<F> = (0..running_input_x.len()).map(|_| F::rand(&mut rng)).collect();
        let current_input_y: Vec<F> = (0..running_input_y.len()).map(|_| F::rand(&mut rng)).collect();

        let running_evaluations: (F, F, F) = shape.inst.inst.evaluate(&running_input_x, &running_input_y);
        let current_evaluations: (F, F, F) = shape.inst.inst.evaluate(&current_input_x, &current_input_y);

        let mut prover_transcript = Transcript::<F>::new(b"new transcript");
        let mut verifier_transcript = prover_transcript.clone();
        let mut verifier_transcript_var = TranscriptVar::from_transcript(
            cs.clone(),
            verifier_transcript.clone(),
        );

        let (beta, matrix_evaluation_proof) = fold_matrices_evaluations(
            &shape,
            (running_input_x.clone(), running_input_y.clone()),
            (current_input_x.clone(), current_input_y.clone()),
            &mut prover_transcript,
            running_evaluations,
            current_evaluations,
            true,
        );

        let verifier = MatrixEvaluationVerifier {
            running_input: (running_input_x.clone(), running_input_y.clone()),
            current_input: (current_input_x.clone(), current_input_y.clone()),
            running_evaluations,
            current_evaluations,
            proof: matrix_evaluation_proof,
        };

        let verifier_var = MatrixEvaluationVerifierVar::new_variable(
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

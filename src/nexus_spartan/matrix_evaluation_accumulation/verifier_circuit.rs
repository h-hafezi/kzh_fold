use std::borrow::Borrow;
use std::ops::Mul;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixEvaluationVerifier<F: PrimeField + Absorb> {
    pub running_input: (Vec<F>, Vec<F>),
    pub new_input: (Vec<F>, Vec<F>),
    pub running_evaluations: (F, F, F),
    pub new_evaluations: (F, F, F),
    pub proof: (F, F, F)
}

impl<F: PrimeField + Absorb> MatrixEvaluationVerifier<F> {
    pub fn accumulate(&self, beta: F) -> ((Vec<F>, Vec<F>), (F, F, F)) {
        // parse variables
        let (r_x, r_y) = self.running_input.clone();
        let (r_x_prime, r_y_prime) = self.new_input.clone();

        // Perform the random combination for r_x_folded and r_y_folded
        let r_x_folded: Vec<F> = r_x.iter()
            .zip(r_x_prime.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        let r_y_folded: Vec<F> = r_y.iter()
            .zip(r_y_prime.iter())
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
        let next_evaluation = expected_eval(self.running_evaluations, self.new_evaluations, self.proof);

        ((r_x_folded, r_y_folded), next_evaluation)
    }
}

pub struct MatrixEvaluationVerifierVar<F: PrimeField + Absorb> {
    pub running_input: (Vec<FpVar<F>>, Vec<FpVar<F>>),
    pub new_input: (Vec<FpVar<F>>, Vec<FpVar<F>>),
    pub running_evaluations: (FpVar<F>, FpVar<F>, FpVar<F>),
    pub new_evaluations: (FpVar<F>, FpVar<F>, FpVar<F>),
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
        let new_input_x = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.new_input.0.clone()), mode)?;
        let new_input_y = Vec::<FpVar<F>>::new_variable(cs.clone(), || Ok(instance.new_input.1.clone()), mode)?;

        let running_eval_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.running_evaluations.0), mode)?;
        let running_eval_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.running_evaluations.1), mode)?;
        let running_eval_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.running_evaluations.2), mode)?;

        let new_eval_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.new_evaluations.0), mode)?;
        let new_eval_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.new_evaluations.1), mode)?;
        let new_eval_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.new_evaluations.2), mode)?;

        let proof_0 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.0), mode)?;
        let proof_1 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.1), mode)?;
        let proof_2 = FpVar::<F>::new_variable(cs.clone(), || Ok(instance.proof.2), mode)?;

        Ok(MatrixEvaluationVerifierVar {
            running_input: (running_input_x, running_input_y),
            new_input: (new_input_x, new_input_y),
            running_evaluations: (running_eval_0, running_eval_1, running_eval_2),
            new_evaluations: (new_eval_0, new_eval_1, new_eval_2),
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
            new_input: (
                self.new_input.0.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?,
                self.new_input.1.iter().map(|var| var.value()).collect::<Result<Vec<_>, _>>()?
            ),
            running_evaluations: (
                self.running_evaluations.0.value()?,
                self.running_evaluations.1.value()?,
                self.running_evaluations.2.value()?,
            ),
            new_evaluations: (
                self.new_evaluations.0.value()?,
                self.new_evaluations.1.value()?,
                self.new_evaluations.2.value()?,
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
    pub fn accumulate(&self, beta: &FpVar<F>) ->((Vec<FpVar<F>>, Vec<FpVar<F>>), (FpVar<F>, FpVar<F>, FpVar<F>)) {
        // Parse variables
        let (r_x, r_y) = &self.running_input;
        let (r_x_prime, r_y_prime) = &self.new_input;

        // Compute r_x_folded and r_y_folded
        let r_x_folded: Vec<FpVar<F>> = r_x.iter()
            .zip(r_x_prime.iter())
            .map(|(rx, rx_prime)| {
                rx.mul(&(FpVar::one() - beta.clone())) + rx_prime.mul(beta)
            })
            .collect::<Vec<_>>();

        let r_y_folded: Vec<FpVar<F>> = r_y.iter()
            .zip(r_y_prime.iter())
            .map(|(ry, ry_prime)| {
                ry.mul(&(FpVar::one() - beta.clone())) + ry_prime.mul(beta)
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
            &self.new_evaluations,
            &self.proof,
        );

        ((r_x_folded, r_y_folded), next_evaluation)
    }
}

#[cfg(test)]
pub mod tests {
    use ark_ff::One;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::{ScalarField, E, G1};
    use crate::nexus_spartan::matrix_evaluation_accumulation::prover::fold_matrices_evaluations;
    use crate::nexus_spartan::matrix_evaluation_accumulation::prover::tests::matrix_evaluation_setup;
    use crate::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationVerifier, MatrixEvaluationVerifierVar};

    type F = ScalarField;

    #[test]
    fn test_non_zk() {
        let (shape, r_x, r_y) = matrix_evaluation_setup::<F, E, G1>();

        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let r_x_prime: Vec<F> = (0..r_x.len()).map(|_| F::rand(&mut rng)).collect();
        let r_y_prime: Vec<F> = (0..r_y.len()).map(|_| F::rand(&mut rng)).collect();

        let z: (F, F, F) = shape.inst.inst.evaluate(&r_x, &r_y);
        let z_prime: (F, F, F) = shape.inst.inst.evaluate(&r_x_prime, &r_y_prime);

        let beta = F::rand(&mut rng);

        let proof = fold_matrices_evaluations(
            &shape,
            (r_x.clone(), r_y.clone()),
            (r_x_prime.clone(), r_y_prime.clone()),
            beta,
            (z, z_prime),
            true,
        );

        let verifier = MatrixEvaluationVerifier{
            running_input: (r_x.clone(), r_y.clone()),
            new_input: (r_x_prime.clone(), r_y_prime.clone()),
            running_evaluations: z,
            new_evaluations: z_prime,
            proof,
        };

        let ((r_x_folded, r_y_folded), new_evaluations) = verifier.accumulate(beta);

        // Perform the random combination for r_x_folded and r_y_folded
        let r_x_folded_expected: Vec<F> = r_x.iter()
            .zip(r_x_prime.iter())
            .map(|(rx, rx_prime)| *rx * (F::one() - beta) + *rx_prime * beta)
            .collect();

        assert_eq!(r_x_folded_expected, r_x_folded);

        let r_y_folded_expected: Vec<F> = r_y.iter()
            .zip(r_y_prime.iter())
            .map(|(ry, ry_prime)| *ry * (F::one() - beta) + *ry_prime * beta)
            .collect();

        assert_eq!(r_y_folded_expected, r_y_folded);

        let new_evaluations_expected: (F, F, F) = shape.inst.inst.evaluate(&r_x_folded, &r_y_folded);

        assert_eq!(new_evaluations_expected, new_evaluations);
    }

    #[test]
    fn test_zk() {
        let (shape, r_x, r_y) = matrix_evaluation_setup::<F, E, G1>();

        let mut rng = thread_rng();

        // Generate random elements in the field for r_x_prime and r_y_prime
        let r_x_prime: Vec<F> = (0..r_x.len()).map(|_| F::rand(&mut rng)).collect();
        let r_y_prime: Vec<F> = (0..r_y.len()).map(|_| F::rand(&mut rng)).collect();

        let z: (F, F, F) = shape.inst.inst.evaluate(&r_x, &r_y);
        let z_prime: (F, F, F) = shape.inst.inst.evaluate(&r_x_prime, &r_y_prime);

        let beta = F::rand(&mut rng);

        let proof = fold_matrices_evaluations(
            &shape,
            (r_x.clone(), r_y.clone()),
            (r_x_prime.clone(), r_y_prime.clone()),
            beta,
            (z, z_prime),
            true,
        );

        let verifier = MatrixEvaluationVerifier{
            running_input: (r_x.clone(), r_y.clone()),
            new_input: (r_x_prime.clone(), r_y_prime.clone()),
            running_evaluations: z,
            new_evaluations: z_prime,
            proof,
        };

        let cs = ConstraintSystem::<F>::new_ref();

        let verifier_var = MatrixEvaluationVerifierVar::new_variable(
            cs.clone(),
            || Ok(verifier.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // assert value function works correctly
        assert_eq!(verifier, verifier_var.value().unwrap());

        let beta = F::rand(&mut thread_rng());
        let beta_var = FpVar::new_variable(
            cs.clone(),
            || Ok(beta.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // call accumulate and see if equal to the non-zk version
        let ((r_x_folded, r_y_folded), new_evaluations) = verifier.accumulate(beta);
        let ((r_x_folded_var, r_y_folded_var), new_evaluations_var) = verifier_var.accumulate(&beta_var);

        // check cs status
        assert!(cs.is_satisfied().unwrap());
        println!("constraint num: {}", cs.num_constraints());

        let r_x_folded_var_value: Vec<F> = r_x_folded_var.iter().map(|val| val.value().unwrap()).collect();
        let r_y_folded_var_value: Vec<F> = r_y_folded_var.iter().map(|val| val.value().unwrap()).collect();

        let new_evaluations_var_value: (F, F, F) = (
            new_evaluations_var.0.value().unwrap(),
            new_evaluations_var.1.value().unwrap(),
            new_evaluations_var.2.value().unwrap(),
        );

        // Assert that the vectors are equal
        assert_eq!(r_x_folded, r_x_folded_var_value, "r_x_folded vectors are not equal");
        assert_eq!(r_y_folded, r_y_folded_var_value, "r_y_folded vectors are not equal");
        assert_eq!(new_evaluations, new_evaluations_var_value, "new_evaluations tuples are not equal");
    }
}
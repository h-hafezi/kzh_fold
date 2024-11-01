use std::borrow::Borrow;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use crate::polynomial::univariate::univariate::PolynomialInterpolator;

/// Struct for evaluating a polynomial in the constraint system with `FpVar<F>` coefficients.
#[derive(Debug)]
pub struct PolynomialInterpolatorVar<F: PrimeField> {
    pub(crate) coefficients: Vec<FpVar<F>>,  // Coefficients in FpVar<F> for the polynomial in the constraint system.
}

impl<F: PrimeField + Absorb> AllocVar<PolynomialInterpolator<F>, F> for PolynomialInterpolatorVar<F> {
    fn new_variable<T: Borrow<PolynomialInterpolator<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        // Fetch the vector of coefficients to allocate as FpVars
        let binding = f()?;
        let values = binding.borrow();

        // Allocate each coefficient as an FpVar
        let coefficients = values
            .coefficients
            .iter()
            .map(|v| FpVar::new_variable(cs.clone(), || Ok(v), mode))
            .collect::<Result<Vec<_>, SynthesisError>>()?;

        Ok(PolynomialInterpolatorVar { coefficients })
    }
}

impl<F: PrimeField + Absorb> R1CSVar<F> for PolynomialInterpolatorVar<F> {
    type Value = PolynomialInterpolator<F>;

    // This method returns the constraint system that the variables are attached to
    fn cs(&self) -> ConstraintSystemRef<F> {
        let mut result = ConstraintSystemRef::None;
        for coeff in &self.coefficients {
            result = coeff.cs().or(result);
        }
        result
    }

    // This method returns the underlying values of the variables, if available
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let mut coeffs = Vec::new();
        for c in &self.coefficients {
            coeffs.push(c.value()?);
        }
        Ok(PolynomialInterpolator { coefficients: coeffs })
    }
}

impl<F: PrimeField> PolynomialInterpolatorVar<F> {

    /// Evaluates the polynomial at a given input `x` in `FpVar<F>`.
    pub fn evaluate(&self, x: &FpVar<F>) -> Result<FpVar<F>, SynthesisError> {
        let mut result = FpVar::<F>::zero();  // Initialize the result to zero.
        let mut x_power = FpVar::<F>::one();  // Start with x^0 = 1.

        for coeff in &self.coefficients {
            // Accumulate the term coefficient * x^i
            result += coeff * &x_power;
            // Update x_power to the next power of x (x^(i+1))
            x_power = x_power * x;
        }

        Ok(result)
    }
}

#[cfg(test)]
mod test {
    use ark_ff::{One, Zero};
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::{ConstraintSystem, SynthesisError};
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::ScalarField;
    use crate::polynomial::univariate::univariate::PolynomialInterpolator;
    use crate::polynomial::univariate::univariate_var::PolynomialInterpolatorVar;

    #[test]
    fn test_polynomial_interpolator_var() -> Result<(), SynthesisError> {
        type F = ScalarField;

        // Create a random polynomial with known coefficients.
        let coefficients = vec![
            F::from(9u8),
            F::from(1u8),
            F::from(100u8),
            F::from(5u8),
            F::zero(),
            F::zero(),
            F::from(12u8),
            F::one(),
        ];
        let original_poly = PolynomialInterpolator {
            coefficients: coefficients.clone(),
        };

        // Set up a constraint system.
        let cs = ConstraintSystem::<F>::new_ref();

        // Allocate the PolynomialInterpolatorVar in the constraint system.
        let poly_var = PolynomialInterpolatorVar::new_variable(
            cs.clone(),
            || Ok(original_poly.clone()),
            ark_r1cs_std::alloc::AllocationMode::Witness,
        )?;

        // 1) Check that calling `.value()` on the allocated variable gives back the original coefficients.
        let reconstructed_poly = poly_var.value().unwrap();
        assert_eq!(reconstructed_poly.coefficients, coefficients);

        // 2) Evaluate both the original and variable polynomials at a random point and compare results.
        let x = F::rand(&mut thread_rng());
        let x_var = FpVar::new_witness(cs.clone(), || Ok(x))?;

        // Evaluate the original polynomial at x.
        let expected_result = original_poly.evaluate(x);

        // Evaluate the variable polynomial at x_var.
        let actual_result_var = poly_var.evaluate(&x_var)?.value()?;

        // Assert that both evaluations give the same result.
        assert_eq!(expected_result, actual_result_var);

        Ok(())
    }

}


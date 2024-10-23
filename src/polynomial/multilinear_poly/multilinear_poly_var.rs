use crate::math::Math;
use crate::polynomial::eq_poly::eq_poly_var::EqPolynomialVar;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use itertools::izip;
use std::borrow::Borrow;

pub struct MultilinearPolynomialVar<F: PrimeField> {
    /// the number of variables in the multilinear polynomial
    pub num_variables: usize,
    /// evaluations of the polynomial in all the 2^num_vars Boolean inputs
    pub evaluation_over_boolean_hypercube: Vec<FpVar<F>>,
    /// length of Z = 2^num_vars
    pub len: usize,
}

impl<F: PrimeField + Absorb> AllocVar<MultilinearPolynomial<F>, F> for MultilinearPolynomialVar<F> {
    fn new_variable<T: Borrow<MultilinearPolynomial<F>>>(
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
        let evaluation_over_boolean_hypercube = values
            .evaluation_over_boolean_hypercube
            .iter()
            .map(|v| FpVar::new_variable(cs.clone(), || Ok(v), mode))
            .collect::<Result<Vec<_>, SynthesisError>>()?;

        Ok(MultilinearPolynomialVar {
            num_variables: values.num_variables,
            evaluation_over_boolean_hypercube,
            len: values.len,
        })
    }
}

impl<F: PrimeField + Absorb> R1CSVar<F> for MultilinearPolynomialVar<F> {
    type Value = MultilinearPolynomial<F>;

    // This method returns the constraint system that the variables are attached to
    fn cs(&self) -> ConstraintSystemRef<F> {
        let mut result = ConstraintSystemRef::None;
        for val in &self.evaluation_over_boolean_hypercube {
            result = val.cs().or(result);
        }
        result
    }

    // This method returns the underlying values of the variables, if available
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let mut evaluation_over_boolean_hypercube = Vec::new();
        for val in &self.evaluation_over_boolean_hypercube {
            evaluation_over_boolean_hypercube.push(val.value()?);
        }
        Ok(MultilinearPolynomial {
            num_variables: self.num_variables,
            evaluation_over_boolean_hypercube,
            len: self.len,
        })
    }
}

impl<F: PrimeField> MultilinearPolynomialVar<F> {
    pub fn new(evaluation_over_boolean_hypercube: Vec<FpVar<F>>) -> Self {
        MultilinearPolynomialVar {
            num_variables: evaluation_over_boolean_hypercube.len().log_2(),
            len: evaluation_over_boolean_hypercube.len(),
            evaluation_over_boolean_hypercube,
        }
    }

    pub fn evaluate(&self, r: &[FpVar<F>]) -> FpVar<F> {
        // r must have a value for each variable
        assert_eq!(r.len(), self.num_variables);

        let chis = EqPolynomialVar::new(r.to_vec()).evals();
        assert_eq!(chis.len(), self.evaluation_over_boolean_hypercube.len());

        let output = {
            let mut res = FpVar::zero();
            for (a, b) in izip!(&self.evaluation_over_boolean_hypercube, chis) {
                res = res + a * b;
            }
            res
        };

        output
    }
}

#[cfg(test)]
mod tests {
    use crate::constant_for_curves::ScalarField;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use crate::polynomial::multilinear_poly::multilinear_poly_var::MultilinearPolynomialVar;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::polynomial::{field_vector_into_fpvar, get_random_vector};

    type F = ScalarField;

    #[test]
    fn test_poly_var() {
        let num_variables = 5;

        let polynomial = MultilinearPolynomial::rand(num_variables, &mut thread_rng());
        let cs = ConstraintSystem::<F>::new_ref();
        let polynomial_var = MultilinearPolynomialVar::new_variable(
            cs.clone(),
            || Ok(polynomial.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // test .value() function to see if works correctly
        assert_eq!(polynomial, polynomial_var.value().unwrap());

        let r :Vec<F> = get_random_vector(num_variables, &mut thread_rng());
        let r_var = field_vector_into_fpvar(cs.clone(), r.as_slice());

        // test .evaluate() works correctly
        assert_eq!(
            polynomial.evaluate(r.as_slice()),
            polynomial_var.evaluate(r_var.as_slice()).value().unwrap()
        )
    }
}
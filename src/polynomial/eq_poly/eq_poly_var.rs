use crate::math::Math;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::{Field, PrimeField};
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

pub struct EqPolynomialVar<F: PrimeField + Copy> {
    r: Vec<FpVar<F>>,
}

impl<F: PrimeField + Absorb> AllocVar<EqPolynomial<F>, F> for EqPolynomialVar<F> {
    fn new_variable<T: Borrow<EqPolynomial<F>>>(
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
        let r = values
            .r
            .iter()
            .map(|v| FpVar::new_variable(cs.clone(), || Ok(v), mode))
            .collect::<Result<Vec<_>, SynthesisError>>()?;

        Ok(EqPolynomialVar { r })
    }
}

impl<F: PrimeField + Absorb> R1CSVar<F> for EqPolynomialVar<F> {
    type Value = EqPolynomial<F>;

    // This method returns the constraint system that the variables are attached to
    fn cs(&self) -> ConstraintSystemRef<F> {
        let mut result = ConstraintSystemRef::None;
        for val in &self.r {
            result = val.cs().or(result);
        }
        result
    }

    // This method returns the underlying values of the variables, if available
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let mut r = Vec::new();
        for val in &self.r {
            r.push(val.value()?);
        }
        Ok(EqPolynomial { r })
    }
}

impl<F: PrimeField + Copy> EqPolynomialVar<F> {
    /// Creates a new EqPolynomial from a vector `w`
    pub fn new(r: Vec<FpVar<F>>) -> Self {
        EqPolynomialVar { r }
    }

    /// Evaluates the polynomial eq_w(r) = prod_{i} (w_i * r_i + (F::ONE - w_i) * (F::ONE - r_i))
    pub fn evaluate(&self, rx: &[FpVar<F>]) -> FpVar<F> {
        assert_eq!(self.r.len(), rx.len());

        let mut result = FpVar::one();

        for i in 0..rx.len() {
            let term = &self.r[i] * &rx[i] + (FpVar::one() - &self.r[i]) * (FpVar::one() - &rx[i]);
            result *= term;
        }

        result
    }

    pub fn evals(&self) -> Vec<FpVar<F>> {
        let ell = self.r.len();

        let mut evals: Vec<FpVar<F>> = vec![FpVar::one(); ell.pow2()];
        let mut size = 1;
        for j in 0..ell {
            // in each iteration, we double the size of chis
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2].clone();
                evals[i] = &scalar * &self.r[j];
                evals[i - 1] = &scalar - &evals[i];
            }
        }
        evals
    }
}

#[cfg(test)]
mod tests {
    use crate::constant_for_curves::ScalarField;
    use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
    use crate::polynomial::eq_poly::eq_poly_var::EqPolynomialVar;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;
    use rand::thread_rng;

    type F = ScalarField;

    #[test]
    fn test_eq_poly_var() {
        let w = vec![
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng())
        ];
        let eq_poly = EqPolynomial::new(w);

        let cs = ConstraintSystem::<F>::new_ref();

        let eq_poly_var = EqPolynomialVar::new_variable(
            cs.clone(),
            || Ok(eq_poly.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // test .value() function works correctly
        assert_eq!(eq_poly, eq_poly_var.value().unwrap());

        let r = vec![
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng())
        ];

        let r_var = vec![
            FpVar::new_variable(cs.clone(), || Ok(r[0].clone()), AllocationMode::Witness).unwrap(),
            FpVar::new_variable(cs.clone(), || Ok(r[1].clone()), AllocationMode::Witness).unwrap(),
            FpVar::new_variable(cs.clone(), || Ok(r[2].clone()), AllocationMode::Witness).unwrap(),
        ];

        // assert that evaluate function works correctly
        assert_eq!(
            eq_poly.evaluate(r.as_slice()),
            eq_poly_var.evaluate(r_var.as_slice()).value().unwrap()
        );

        // assert that evals() function works correctly
        assert_eq!(
            eq_poly.evals(),
            {
                let eval_values: Vec<F> = eq_poly_var.evals().iter()
                    .map(|eval| eval.value().unwrap())
                    .collect();
                eval_values
            }
        )
    }
}
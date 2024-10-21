pub mod sparse_mlpoly;

use std::borrow::Borrow;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::CurveConfig;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use ark_serialize::CanonicalSerialize;
use crate::math::Math;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparsePoly<F: Absorb> {
    num_vars: usize,
    evals: Vec<F>,
}

impl<F: PrimeField + Absorb> SparsePoly<F> {
    pub fn new(num_vars: usize, evals: Vec<F>) -> Self {
        SparsePoly { num_vars, evals }
    }

    fn compute_chi(a: &[bool], r: &[F]) -> F {
        assert_eq!(a.len(), r.len());
        let mut chi_i = F::one();
        for j in 0..r.len() {
            if a[j] {
                chi_i *= r[j];
            } else {
                chi_i *= F::one() - r[j];
            }
        }
        chi_i
    }

    pub fn evaluate(&self, r: &[F]) -> F {
        assert_eq!(self.num_vars, r.len());

        (0..self.evals.len())
            .map(|i| {
                let bits = i.get_bits_canonical_order(r.len());
                SparsePoly::compute_chi(&bits, r) * self.evals[i]
            })
            .sum()
    }
}

pub struct SparsePolyVar<F: Absorb + PrimeField> {
    num_vars: usize,
    evals: Vec<FpVar<F>>,
}

// Implement AllocVar for SparsePoly
impl<F: PrimeField + Absorb> AllocVar<SparsePoly<F>, F> for SparsePolyVar<F> {
    fn new_variable<T: Borrow<SparsePoly<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        // Convert to Namespace<F>
        let ns = cs.into();
        // Get the constraint system reference
        let cs = ns.cs();

        // Retrieve the SparsePoly<F> object
        let binding = f()?;
        let sparse_poly = binding.borrow();

        // Allocate each field element in the `evals` as circuit variables
        let mut evals_var = Vec::new();
        for i in 0..sparse_poly.evals.len() {
            evals_var.push(FpVar::new_variable(
                cs.clone(),
                || Ok(sparse_poly.evals[i]),
                mode
            )?);
        }

        // Return the allocated SparsePolyVar
        Ok(SparsePolyVar {
            num_vars: sparse_poly.num_vars,
            evals: evals_var,
        })
    }
}

impl<F: PrimeField + Absorb> R1CSVar<F> for SparsePolyVar<F> {
    type Value = SparsePoly<F>;

    // Return the constraint system associated with the variable
    fn cs(&self) -> ConstraintSystemRef<F> {
        self.evals
            .iter()
            .fold(ConstraintSystemRef::None, |cs, eval| cs.or(eval.cs()))
    }

    // Extract the underlying values of the evals and return a SparsePoly<F>
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        // Retrieve the values of the evals (they are Option<F> inside FpVar<F>)
        let evals = self
            .evals
            .iter()
            .map(|eval_var| eval_var.value()) // Extract value for each eval_var
            .collect::<Result<Vec<_>, SynthesisError>>()?; // Collect them into a Vec<F>

        // Return a new SparsePoly<F> with the same number of variables
        Ok(SparsePoly {
            num_vars: self.num_vars,
            evals,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_r1cs_std::prelude::*;
    use ark_ff::UniformRand;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef, SynthesisError};
    use rand::thread_rng;

    // Define the ScalarField type (replace with your actual field type)
    type F = ark_bls12_381::Fr; // Assuming you're using BLS12-381 scalar field

    #[test]
    fn test_sparse_poly_var() -> Result<(), SynthesisError> {
        // Create a new constraint system
        let cs = ConstraintSystem::<F>::new_ref();

        // Generate a random sparse polynomial
        let num_vars = 3;
        let evals: Vec<F> = (0..5)
            .map(|_| F::rand(&mut thread_rng()))
            .collect();
        let sparse_poly = SparsePoly::new(num_vars, evals.clone());

        // Generate a SparsePolyVar from the sparse_poly
        let sparse_poly_var = SparsePolyVar::new_variable(
            cs.clone(),
            || Ok(sparse_poly.clone()),
            AllocationMode::Witness,
        )?;

        // Print the number of constraints
        let num_constraints = cs.num_constraints();
        println!("Number of constraints: {}", num_constraints);

        // Check that value() returns the correct sparse polynomial
        let sparse_poly_value = sparse_poly_var.value()?;
        assert_eq!(sparse_poly_value, sparse_poly);

        Ok(())
    }
}

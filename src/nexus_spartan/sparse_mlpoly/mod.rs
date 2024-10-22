pub mod sparse_mlpoly;
mod util;

use crate::math::Math;
use crate::nexus_spartan::sparse_mlpoly::util::get_bits_canonical_order;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::CurveConfig;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::select::CondSelectGadget;
use ark_r1cs_std::R1CSVar;
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use ark_serialize::CanonicalSerialize;
use std::borrow::Borrow;

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

impl<F: PrimeField + Absorb> SparsePolyVar<F> {
    pub fn new(num_vars: usize, evals: Vec<FpVar<F>>) -> SparsePolyVar<F> {
        SparsePolyVar {
            num_vars,
            evals,
        }
    }
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

        let res = f();
        let sparse_poly = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);

        // Allocate each field element in the `evals` as circuit variables
        let mut evals_var = Vec::new();
        for i in 0..sparse_poly?.evals.len() {
            evals_var.push(FpVar::new_variable(
                ns!(cs, "y"),
                || Ok(sparse_poly?.evals[i]),
                mode,
            )?);
        }

        // Return the allocated SparsePolyVar
        Ok(SparsePolyVar {
            num_vars: sparse_poly?.num_vars,
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

impl<F: PrimeField + Absorb> SparsePolyVar<F> {
    pub fn compute_chi(a: &[Boolean<F>], r: &[FpVar<F>]) -> FpVar<F> {
        assert_eq!(a.len(), r.len(), "Unequal length vectors");
        let mut chi_i = FpVar::one();
        for j in 0..r.len() {
            let val = FpVar::conditionally_select(&a[j], &r[j], &(FpVar::one() - &r[j])).unwrap();
            chi_i = chi_i * val;
        }
        chi_i
    }

    pub fn evaluate(&self, r: &[FpVar<F>]) -> FpVar<F> {
        assert_eq!(self.num_vars, r.len());
        let mut res = FpVar::<F>::zero();
        let mut counter = FpVar::<F>::zero();
        for i in 0..self.evals.len() {
            let bits = get_bits_canonical_order(&counter, r.len());
            counter = counter + FpVar::one();
            res += SparsePolyVar::compute_chi(&bits, r) * &self.evals[i];
        }
        res
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::constant_for_curves::ScalarField;
    use ark_ff::UniformRand;
    use ark_r1cs_std::prelude::*;
    use ark_relations::r1cs::ConstraintSystem;
    use rand::thread_rng;

    type F = ScalarField;

    fn get_random_sparse_poly() -> SparsePoly<F> {
        // Generate a random sparse polynomial
        let num_vars = 3;
        let evals: Vec<F> = (0..5)
            .map(|_| F::rand(&mut thread_rng()))
            .collect();
        SparsePoly::new(num_vars, evals.clone())
    }

    #[test]
    fn test_sparse_poly_var() {
        // Create a new constraint system
        let cs = ConstraintSystem::<F>::new_ref();

        let sparse_poly = get_random_sparse_poly();

        // Generate a SparsePolyVar from the sparse_poly
        let sparse_poly_var = SparsePolyVar::new_variable(
            cs.clone(),
            || Ok(sparse_poly.clone()),
            AllocationMode::Witness,
        ).unwrap();

        ///////////////////////////////// test .value()

        // Check that value() returns the correct sparse polynomial
        let sparse_poly_value = sparse_poly_var.value().unwrap();
        assert_eq!(sparse_poly_value, sparse_poly);

        ///////////////////////////////// test .compute_chi()

        // compute chi for both instances
        let (a1, a2, a3) = (
            bool::rand(&mut thread_rng()),
            bool::rand(&mut thread_rng()),
            bool::rand(&mut thread_rng())
        );
        let (a1_var, a2_var, a3_var) = (
            Boolean::new_witness(cs.clone(), || Ok(a1)).unwrap(),
            Boolean::new_witness(cs.clone(), || Ok(a2)).unwrap(),
            Boolean::new_witness(cs.clone(), || Ok(a3)).unwrap(),
        );

        let (r1, r2, r3) = (
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng())
        );
        let (r1_var, r2_var, r3_var) = (
            FpVar::new_witness(cs.clone(), || Ok(r1)).unwrap(),
            FpVar::new_witness(cs.clone(), || Ok(r2)).unwrap(),
            FpVar::new_witness(cs.clone(), || Ok(r3)).unwrap(),
        );

        let t = SparsePoly::compute_chi(&[a1, a2, a3], &[r1, r2, r3]);

        let t_var = SparsePolyVar::compute_chi(
            &[a1_var, a2_var, a3_var],
            &[r1_var.clone(), r2_var.clone(), r3_var.clone()],
        );

        assert_eq!(t, t_var.value().unwrap());

        ///////////////////////////////// test .evaluate()

        let res = sparse_poly.evaluate(&[r1, r2, r3]);
        let res_var = sparse_poly_var.evaluate(&[r1_var, r2_var, r3_var]);

        assert_eq!(res, res_var.value().unwrap());
    }
}

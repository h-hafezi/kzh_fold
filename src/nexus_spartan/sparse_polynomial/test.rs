#[cfg(test)]
mod test {
    use crate::constant_for_curves::ScalarField;
    use ark_ff::UniformRand;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::prelude::*;
    use ark_relations::r1cs::ConstraintSystem;
    use rand::thread_rng;
    use crate::nexus_spartan::sparse_polynomial::sparse_polynomial::SparsePoly;
    use crate::nexus_spartan::sparse_polynomial::sparse_polynomial_var::SparsePolyVar;

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

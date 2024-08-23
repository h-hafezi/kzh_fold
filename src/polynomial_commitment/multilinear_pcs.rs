/// Implementation of the multilinear KZH

#[cfg(test)]
mod tests {
  use super::*;
  use crate::spartan::dense_mlpoly::DensePolynomial as DenseMLPoly;
  use ark_ff::{One, Zero};

    use crate::constant_for_curves::{E, ScalarField};

  #[test]
  fn test_multilinear() {
      // Generate p(X,Y) from evaluations, such that:
      //  p(0,0) = 0
      //  p(0,1) = 1
      //  p(1,0) = 2
      //  p(1,1) = 3
      let p_coeffs = (0..4u8).map(ScalarField::from).collect::<Vec<_>>();
      let mut p = DenseMLPoly::<ScalarField>::new(p_coeffs);
      assert_eq!(p.get_num_vars(), 2);

      // Test p(0,0) == 0
      assert_eq!(
          p.evaluate(&[ScalarField::zero(), ScalarField::zero()]),
          ScalarField::from(0u8)
      );

      // Test p(0,1) == 1
      assert_eq!(
          p.evaluate(&[ScalarField::zero(), ScalarField::one()]),
          ScalarField::from(1u8)
      );

      // Test p(1,0) == 2
      assert_eq!(
          p.evaluate(&[ScalarField::one(), ScalarField::zero()]),
          ScalarField::from(2u8)
      );

      // Test p(1,1) == 3
      assert_eq!(
          p.evaluate(&[ScalarField::one(), ScalarField::one()]),
          ScalarField::from(3u8)
      );

      // Get polynomial from the partial evaluation of p(0,Y)
      let r = ScalarField::zero();
      p.bound_poly_var_top(&r);
      assert_eq!(p.get_num_vars(), 1);

      // Now evaluate p(0,Y) to 1: Hence compute p(0,1) == 1
      assert_eq!(
          p.evaluate(&[ScalarField::one()]),
          ScalarField::from(1u8)
      );
  }
}

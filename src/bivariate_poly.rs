use std::fmt;
use ark_ff::{Field, Zero, PrimeField, FftField};
use itertools::Itertools;
use rand::RngCore;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use crate::lagrange_basis::{LagrangeBasis, LagrangeTraits};
use crate::univariate_poly::UnivariatePolynomial;
use crate::utils::{compute_powers, is_power_of_two};


/// A bivariate polynomial can be represented in two different forms:
///
/// 1. **Coefficient Form**:
///
///    f(X, Y) = sum_{i=0}^{d-1} sum_{j=0}^{d-1} f_{i,j} * X^i * Y^j
///
///    The coefficients of the bivariate polynomial are stored in a flat vector and logically
///    organized in a matrix of dimensions `d x d`. The matrix representation
///    allows us to associate the coefficients with the powers of X and Y.
///
///    For example, given a polynomial with 16 coefficients and degree `d = 4`:
///
///    f(X, Y) = f_{0,0} + f_{0,1}Y + f_{0,2}Y^2 + f_{0,3}Y^3 +
///              f_{1,0}X + f_{1,1}XY + f_{1,2}XY^2 + f_{1,3}XY^3 +
///              f_{2,0}X^2 + f_{2,1}X^2Y + f_{2,2}X^2Y^2 + f_{2,3}X^2Y^3 +
///              f_{3,0}X^3 + f_{3,1}X^3Y + f_{3,2}X^3Y^2 + f_{3,3}X^3Y^3
///
///    The coefficients are stored in a vector as:
///    [f_{0,0}, f_{0,1}, f_{0,2}, f_{0,3}, f_{1,0}, f_{1,1}, f_{1,2}, f_{1,3},
///     f_{2,0}, f_{2,1}, f_{2,2}, f_{2,3}, f_{3,0}, f_{3,1}, f_{3,2}, f_{3,3}]
///
///    Each row in the matrix corresponds to increasing powers of Y, and within each row,
///    the columns correspond to increasing powers of X.
///
/// 2. **Lagrange Basis Form**:
///
///    In the Lagrange basis, the polynomial is represented using the evaluation points and
///    corresponding Lagrange basis polynomials:
///
///    f(X, Y) = L_{0,0}(w_0, w_0)f(w_0, w_0) + L_{0,1}(w_0, w_1)f(w_0, w_1) + L_{0,2}(w_0, w_2)f(w_0, w_2) + f(w_0, w_3) +
///              L_{1,0}(w_1, w_0)f(w_1, w_0) + L_{1,1}(w_1, w_1)f(w_1, w_1) + L_{1,2}(w_1, w_2)f(w_1, w_2) + f(w_1, w_3) +
///              L_{2,0}(w_2, w_0)f(w_2, w_0) + L_{2,1}(w_2, w_1)f(w_2, w_1) + L_{2,2}(w_2, w_2)f(w_2, w_2) + f(w_2, w_3) +
///              L_{3,0}(w_3, w_0)f(w_3, w_0) + L_{3,1}(w_3, w_1)f(w_3, w_1) + L_{3,2}(w_3, w_2)f(w_3, w_2) + f(w_3, w_3)
///
///    Here, L_{i,j}(w_i, w_j) are the Lagrange basis polynomials evaluated at the points w_i and w_j, and f(w_i, w_j)
///    are the evaluations of the polynomial at those points. This form is particularly useful for polynomial interpolation.

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BivariatePolynomial<F: FftField> {
    // evaluation[i][j] corresponds to f(w_i, w_j)
    pub evaluations: Vec<Vec<F>>,
    // the lagrange basis used
    pub lagrange_basis_x: LagrangeBasis<F>,
    pub lagrange_basis_y: LagrangeBasis<F>,
    // Degree of the polynomial in both X and Y
    pub degree_x: usize,
    pub degree_y: usize,
}

pub trait BivariatePolynomialTrait<F: FftField> {
    fn new(evaluations: Vec<Vec<F>>,
           domain_x: GeneralEvaluationDomain<F>,
           domain_y: GeneralEvaluationDomain<F>,
           degree_x: usize,
           degree_y: usize,
    ) -> Self;
    fn random<T: RngCore>(rng: &mut T,
                          domain_x: GeneralEvaluationDomain<F>,
                          domain_y: GeneralEvaluationDomain<F>,
                          degree_x: usize,
                          degree_y: usize,
    ) -> Self;
    fn evaluate(&self, x: &F, y: &F) -> F;
    fn partially_evaluate_at_x(&self, x: &F) -> UnivariatePolynomial<F>;
    fn partially_evaluate_at_y(&self, y: &F) -> UnivariatePolynomial<F>;
}

impl<F: FftField + fmt::Display> fmt::Display for BivariatePolynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "f(X, Y) =")?;
        for i in 0..self.degree_x {
            if i > 0 {
                write!(f, "+ ")?;
            }
            let mut first_term = true;
            for j in 0..self.degree_y {
                if !first_term {
                    write!(f, " + ")?;
                }
                write!(f, "f(w_{}, w_{})", i, j)?;
                write!(f, " ({})", self.evaluations[i][j])?;
                first_term = false;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}


impl<F: FftField> BivariatePolynomialTrait<F> for BivariatePolynomial<F> {
    fn new(
        evaluations: Vec<Vec<F>>,
        domain_x: GeneralEvaluationDomain<F>,
        domain_y: GeneralEvaluationDomain<F>,
        degree_x: usize,
        degree_y: usize,
    ) -> Self {
        assert!(is_power_of_two(degree_x), "degree_x (upper bound) must be a power of two");
        assert!(is_power_of_two(degree_y), "degree_y (upper bound) must be a power of two");
        // return the instance
        Self {
            evaluations,
            lagrange_basis_x: LagrangeBasis { domain: domain_x },
            lagrange_basis_y: LagrangeBasis { domain: domain_y },
            degree_x,
            degree_y,
        }
    }

    /// Generates a random BivariatePolynomial
    fn random<T: RngCore>(
        rng: &mut T,
        domain_x: GeneralEvaluationDomain<F>,
        domain_y: GeneralEvaluationDomain<F>,
        degree_x: usize,
        degree_y: usize,
    ) -> Self {
        assert!(is_power_of_two(degree_x), "degree_x (upper bound) must be a power of two");
        assert!(is_power_of_two(degree_y), "degree_y (upper bound) must be a power of two");

        let mut evaluations = Vec::with_capacity(degree_x);

        for _ in 0..degree_x {
            let mut row = Vec::with_capacity(degree_y);
            for _ in 0..degree_y {
                row.push(F::rand(rng));
            }
            evaluations.push(row);
        }

        BivariatePolynomial {
            evaluations,
            lagrange_basis_x: LagrangeBasis { domain: domain_x },
            lagrange_basis_y: LagrangeBasis { domain: domain_y },
            degree_x,
            degree_y,
        }
    }

    /// evaluation requires O(n^2) additions
    fn evaluate(&self, x: &F, y: &F) -> F {
        let l_x = self.lagrange_basis_x.evaluate(x);
        let l_y = self.lagrange_basis_y.evaluate(y);
        // the final result
        let mut sum = F::ZERO;
        for i in 0..self.degree_x {
            for j in 0..self.degree_y {
                sum += l_x[i] * l_y[j] * self.evaluations[i][j];
            }
        }
        sum
    }

    /// f(x, Y) = sum_{i} L_i(x) * sum_{j} (L_j(Y) * f(w_i, w_j)) ===>
    /// f(x, w_t) = sum_{i} L_i(x) * sum_{j} (L_j(w_t) * f(w_i, w_j))
    ///           = sum_{i} L_i(x) * f(w_i, w_t))
    fn partially_evaluate_at_x(&self, x: &F) -> UnivariatePolynomial<F> {
        let l_x = self.lagrange_basis_x.evaluate(x);
        let mut evaluations = vec![];
        for t in 0..self.degree_y {
            let mut sum = F::ZERO;
            for i in 0..self.degree_x {
                sum += l_x[i] * self.evaluations[i][t];
            }
            evaluations.push(sum);
        }
        UnivariatePolynomial { evaluations, lagrange_basis: self.lagrange_basis_y.clone() }
    }

    /// f(X, y) = sum_{j} L_j(y) * sum_{i} (L_i(X) * f(w_i, w_j)) ===>
    /// f(w_t, y) = sum_{j} L_j(y) * sum_{i} (L_i(w_t) * f(w_i, w_j))
    ///           = sum_{j} L_j(y) * f(w_t, w_j))
    fn partially_evaluate_at_y(&self, y: &F) -> UnivariatePolynomial<F> {
        let l_y = self.lagrange_basis_y.evaluate(y);
        let mut evaluations = vec![];
        for t in 0..self.degree_x {
            let mut sum = F::ZERO;
            for j in 0..self.degree_y {
                sum += l_y[j] * self.evaluations[t][j];
            }
            evaluations.push(sum);
        }
        UnivariatePolynomial { evaluations, lagrange_basis: self.lagrange_basis_x.clone() }
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::bivariate_poly::{BivariatePolynomial, BivariatePolynomialTrait};
    use crate::univariate_poly::UnivariatePolynomialTrait;

    type F = Fr;

    #[test]
    fn test_random_bivariate() {
        let degree_x = 4usize;
        let degree_y = 16usize;
        let domain_x = GeneralEvaluationDomain::<F>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<F>::new(degree_y).unwrap();
        let r: BivariatePolynomial<F> = BivariatePolynomial::random(&mut thread_rng(), domain_x, domain_y, degree_x, degree_y);
        println!("{}", r);
    }

    #[test]
    fn test_partial_evaluation_x() {
        let degree_x = 4usize;
        let degree_y = 16usize;
        let domain_x = GeneralEvaluationDomain::<F>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<F>::new(degree_y).unwrap();
        let r: BivariatePolynomial<F> = BivariatePolynomial::random(&mut thread_rng(), domain_x, domain_y, degree_x, degree_y);
        let x = F::rand(&mut thread_rng());
        let y = F::rand(&mut thread_rng());
        let r_x = r.partially_evaluate_at_x(&x);
        let r_xy_indirect = r_x.evaluate(&y);
        let r_xy_direct = r.evaluate(&x, &y);
        assert_eq!(r_xy_direct, r_xy_indirect);
    }

    #[test]
    fn test_partial_evaluation_y() {
        let degree_x = 16usize;
        let degree_y = 4usize;
        let domain_x = GeneralEvaluationDomain::<F>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<F>::new(degree_y).unwrap();
        println!("{} {}", domain_x.size(), domain_y.size());
        let r: BivariatePolynomial<F> = BivariatePolynomial::random(&mut thread_rng(), domain_x, domain_y, degree_x, degree_y);
        let x = F::rand(&mut thread_rng());
        let y = F::rand(&mut thread_rng());
        let r_y = r.partially_evaluate_at_y(&y);
        let r_xy_indirect = r_y.evaluate(&x);
        let r_xy_direct = r.evaluate(&x, &y);
        assert_eq!(r_xy_direct, r_xy_indirect);
    }
}
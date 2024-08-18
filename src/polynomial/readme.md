
---

# Polynomial Implementations 

This repository provides implementations of univariate and bivariate polynomials in Rust, with a focus on working in the Lagrange basis. The library is designed for cryptographic applications, particularly in contexts where efficient polynomial operations are critical.

## Overview

### Lagrange Basis Evaluation

The `LagrangeBasis` struct is initialized using a `GeneralEvaluationDomain` and supports two primary operations: **evaluation** and **vanishing polynomial evaluation**.

```rust
/// Represents a polynomial in the Lagrange basis, supporting evaluation and vanishing polynomial computation.
pub struct LagrangeBasis {
    domain: GeneralEvaluationDomain,
}

impl LagrangeBasis {
    /// Evaluates the Lagrange basis polynomials at a given point.
    pub fn evaluate(&self, point: &Self::ScalarField) -> Vec<Self::ScalarField> {
        // Implementation goes here...
    }

    /// Computes the vanishing polynomial of the domain at a given point.
    pub fn evaluate_vanishing_polynomial(&self, point: &Self::ScalarField) -> Self::ScalarField {
        // Implementation goes here...
    }
}
```

- **Evaluation**: The `evaluate` method outputs a vector containing the values \( Z_{w_i}(z) \), where \( Z_{w_i}(z) \) is the Lagrange basis polynomial that evaluates to 1 at \( w_i \) and 0 at all other points \( w_j \) in the domain.

- **Vanishing Polynomial Evaluation**: The `evaluate_vanishing_polynomial` method computes the value of the vanishing polynomial of the domain at a given point \( z \). The vanishing polynomial is defined as \( Z(X) = \prod_{i=0}^{n-1} (X - w_i) \), where \( w_i \) are the elements of the domain.

### Univariate Polynomials in Lagrange Basis

The `UnivariatePolynomial` struct is implemented on top of `DensePolynomial` for polynomials represented in the Lagrange basis. It provides the following functionalities: **evaluation** and **addition**.

```rust
/// Represents a univariate polynomial in the Lagrange basis.
pub struct UnivariatePolynomial {
    coefficients: DensePolynomial,
}

impl Polynomial for UnivariatePolynomial {
    /// Evaluates the polynomial at a given point.
    fn evaluate(&self, point: &Self::ScalarField) -> Self::ScalarField {
        // Implementation goes here...
    }

    /// Adds two polynomials in the Lagrange basis.
    fn add(&self, other: &Self) -> Self {
        // Implementation goes here...
    }
}
```

- **Evaluate**: The `evaluate` method computes the value of the univariate polynomial at a given point.

- **Add**: The `add` method performs the addition of two polynomials represented in the Lagrange basis.

### Bivariate Polynomials

The `BivariatePolynomial` struct represents bivariate polynomials, which can be structured in two distinct forms: **coefficient form** and **Lagrange basis form**. It also supports **partial evaluation**.

```rust
/// Represents a bivariate polynomial in both coefficient and Lagrange basis forms.
pub struct BivariatePolynomial {
    coefficients: Vec<Self::ScalarField>,
    degree: usize,
}

impl Polynomial for BivariatePolynomial {
    /// Evaluates the bivariate polynomial at a given pair of points.
    fn evaluate(&self, point: &(Self::ScalarField, Self::ScalarField)) -> Self::ScalarField {
        // Implementation goes here...
    }

    /// Adds two bivariate polynomials.
    fn add(&self, other: &Self) -> Self {
        // Implementation goes here...
    }
}

impl PartialEvaluation for BivariatePolynomial {
    /// Partially evaluates a bivariate polynomial at a given X value, resulting in a univariate polynomial in Y.
    fn partial_evaluate_x(&self, x: &Self::ScalarField) -> UnivariatePolynomial {
        // Implementation goes here...
    }

    /// Partially evaluates a bivariate polynomial at a given Y value, resulting in a univariate polynomial in X.
    fn partial_evaluate_y(&self, y: &Self::ScalarField) -> UnivariatePolynomial {
        // Implementation goes here...
    }
}
```

- **Coefficient Form**: The bivariate polynomial \( f(X, Y) \) is represented as a sum of monomials with coefficients:
  \[
  f(X, Y) = \sum_{i=0}^{d-1} \sum_{j=0}^{d-1} f_{i,j} \cdot X^i \cdot Y^j
  \]
  The coefficients \( f_{i,j} \) are stored in a flat vector, logically organized as a \( d \times d \) matrix.

- **Lagrange Basis Form**: In this form, the polynomial is represented using evaluation points and corresponding Lagrange basis polynomials:
  \[
  f(X, Y) = \sum_{i=0}^{d-1} \sum_{j=0}^{d-1} L_{i,j}(w_i, w_j) \cdot f(w_i, w_j)
  \]
  This form is particularly useful for polynomial interpolation.

- **Partial Evaluation**: The `partial_evaluate_x` and `partial_evaluate_y` methods perform partial evaluation of the bivariate polynomial. For a given \( x \), `partial_evaluate_x` outputs a univariate polynomial \( g(Y) = f(x, Y) \). Similarly, `partial_evaluate_y` evaluates \( g(X) = f(X, y) \) for a given \( y \).

## Traits

The library defines several key traits to ensure flexibility and interoperability:

```rust
/// Defines core polynomial operations such as evaluation, addition, and multiplication.
pub trait Polynomial {
    /// Evaluates the polynomial at a given point.
    fn evaluate(&self, point: &Self::ScalarField) -> Self::ScalarField;

    /// Adds two polynomials.
    fn add(&self, other: &Self) -> Self;
}

/// Provides specialized operations for polynomials in the Lagrange basis.
pub trait LagrangeBasisOperations {
    /// Interpolates the polynomial from its evaluations.
    fn interpolate(&self, evaluations: &[Self::ScalarField]) -> Self;

    /// Transforms the polynomial from coefficient form to Lagrange basis form.
    fn from_coefficients(coefficients: &[Self::ScalarField]) -> Self;
}

/// Defines the partial evaluation function for bivariate polynomials.
pub trait PartialEvaluation {
    /// Partially evaluates a bivariate polynomial at a given point.
    fn partial_evaluate_x(&self, x: &Self::ScalarField) -> UnivariatePolynomial;

    /// Partially evaluates a bivariate polynomial at a given point.
    fn partial_evaluate_y(&self, y: &Self::ScalarField) -> UnivariatePolynomial;
}
```


# Polynomial Commitment Scheme

## Overview

This repository implements a **polynomial commitment scheme** tailored for bivariate polynomials, featuring an opening proof size of `degree_x + degree_y`. A polynomial commitment scheme allows one to commit to a polynomial while keeping the polynomial itself hidden, and later provide a proof that the committed value evaluates to a certain result at a given point. This is useful in cryptographic protocols where integrity and privacy of polynomial evaluations are crucial.

## Structs

```rust
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SRS<E: Pairing> {
    pub degree_x: usize,
    pub degree_y: usize,
    pub matrix_H: Vec<Vec<E::G1Affine>>,
    pub vec_H: Vec<E::G1Affine>,
    pub vec_V: Vec<E::G2>,
    pub V_prime: E::G2,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<E: Pairing> {
    pub C: E::G1Affine,
    pub aux: Vec<E::G1>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpeningProof<E: Pairing> {
    pub vec_D: Vec<E::G1Affine>,
    pub f_star_poly: UnivariatePolynomial<E::ScalarField>,
}

// Define the new struct that encapsulates the functionality of polynomial commitment
pub struct PolyCommit<E: Pairing> {
    pub srs: SRS<E>,
}
```

## Traits

```rust
pub trait PolyCommitTrait<E: Pairing> {
    fn setup<T: RngCore>(n: usize, m: usize, rng: &mut T) -> SRS<E>;

    fn commit(&self, poly: &BivariatePolynomial<E::ScalarField>) -> Commitment<E>;

    fn open(&self,
            poly: &BivariatePolynomial<E::ScalarField>,
            com: Commitment<E>,
            b: &E::ScalarField,
    ) -> OpeningProof<E>;

    fn verify(&self,
              lagrange_x: LagrangeBasis<E::ScalarField>,
              C: &Commitment<E>,
              proof: &OpeningProof<E>,
              b: &E::ScalarField,
              c: &E::ScalarField,
              y: &E::ScalarField,
    ) -> bool;
}
```

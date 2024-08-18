---

# Accumulator

This Rust library implements a cryptographic accumulator, harnessing the power of pairing-based cryptography to streamline verification processes.

## Overview

An **accumulator** is a cryptographic construct that enables the aggregation (or "folding") of multiple NP-instance verifications into a single instance verification. This library provides an accumulation scheme specifically for polynomial commitment schemes (PCS), where multiple relaxed PCS openings are consolidated into a single opening.

Key components include:

- **`AccSRS`**: The Setup Reference String (SRS) used by the accumulator.
- **`AccInstance`**: Represents a single instance of the accumulator, containing a relaxed polynomial opening proof.
- **`AccWitness`**: A relaxed polynomial opening.
- **`Accumulator`**: A structure combining both an `AccInstance` and an `AccWitness`.
- **`AccumulatorTrait`**: A trait defining the core functions for interacting with accumulators.

### Structs

#### 1. AccSRS

```rust
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccSRS<E: Pairing> {
    pub degree_x: usize,
    pub degree_y: usize,
    pub lagrange_basis_x: LagrangeBasis<E::ScalarField>,
    pub lagrange_basis_y: LagrangeBasis<E::ScalarField>,
    pub k_vec_b: Vec<E::G1Affine>,
    pub k_vec_c: Vec<E::G1Affine>,
    pub k_prime: E::G1Affine,
    pub pc_srs: SRS<E>,
}
```

#### 2. AccInstance

```rust
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccInstance<E: Pairing> {
    pub C: E::G1Affine,
    pub T: E::G1Affine,
    pub E: E::G1Affine,
    pub b: E::ScalarField,
    pub c: E::ScalarField,
    pub y: E::ScalarField,
    pub z_b: E::ScalarField,
    pub z_c: E::ScalarField,
}
```

The `AccInstance` struct encapsulates an accumulator instance, holding group elements (`C`, `T`, `E`) and scalar fields (`b`, `c`, `y`, `z_b`, `z_c`).

#### 3. AccWitness

```rust
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AccWitness<E: Pairing> {
    pub vec_D: Vec<E::G1Affine>,
    pub f_star_poly: UnivariatePolynomial<E::ScalarField>,
    pub vec_b: Vec<E::ScalarField>,
    pub vec_c: Vec<E::ScalarField>,
}
```

#### 4. Accumulator

```rust
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Accumulator<E: Pairing> {
    pub witness: AccWitness<E>,
    pub instance: AccInstance<E>,
}
```

### Traits

#### AccumulatorTrait

```rust
pub trait AccumulatorTrait<E: Pairing>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
{
    fn setup<T: RngCore>(degree_x: usize, degree_y: usize, lagrange_basis_x: LagrangeBasis<E::ScalarField>, lagrange_basis_y: LagrangeBasis<E::ScalarField>, pc_srs: SRS<E>, rng: &mut T) -> AccSRS<E>;

    fn new_accumulator(instance: &AccInstance<E>, witness: &AccWitness<E>) -> Accumulator<E>;

    fn compute_randomness(instance_1: &AccInstance<E>, instance_2: &AccInstance<E>, Q: E::G1Affine) -> E::ScalarField;

    fn new_accumulator_instance_from_proof(srs: &AccSRS<E>, C: &E::G1Affine, b: &E::ScalarField, c: &E::ScalarField, y: &E::ScalarField) -> AccInstance<E>;

    fn new_accumulator_witness_from_proof(srs: &AccSRS<E>, proof: OpeningProof<E>, b: &E::ScalarField, c: &E::ScalarField) -> AccWitness<E>;

    fn prove(srs: &AccSRS<E>, acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> (AccInstance<E>, AccWitness<E>, E::G1Affine);

    fn verify(instance_1: &AccInstance<E>, instance_2: &AccInstance<E>, Q: E::G1Affine) -> AccInstance<E>;

    fn decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> bool;

    fn helper_function_decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> E::G1Affine;

    fn helper_function_Q(srs: &AccSRS<E>, acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> E::G1Affine;
}
```

The `AccumulatorTrait` outlines the essential methods for working with accumulators:

- `setup`: Generates the `AccSRS` using the provided parameters.
- `new_accumulator`: Constructs a new `Accumulator` from an `AccInstance` and its corresponding `AccWitness`.
- `compute_randomness`: Calculates a randomness value from two instances and a group element `Q`, facilitating non-interactive Fiat-Shamir challenges.
- `new_accumulator_instance_from_proof`: Creates a new `AccInstance` from a polynomial commitment scheme (PCS) proof.
- `new_accumulator_witness_from_proof`: Constructs a new `AccWitness` from a PCS opening.
- `prove`: Combines two `Accumulator` objects (including their instances and witnesses) into a single accumulator.
- `verify`: Aggregates two `AccInstance` objects into one using the accumulation proof `Q`.
- `decide`: Verifies the correctness of an `Accumulator`, ensuring the witness is valid concerning the instance.
- Helper functions (`helper_function_decide` and `helper_function_Q`) aid in simplifying the code structure by isolating specific operations into standalone functions.

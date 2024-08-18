
# R1CS

This code has been borrowed from [Nexus](https://github.com/nexus-xyz/nexus-zkvm). It implements R1CS (Rank-1 Constraint System) instances, witnesses, and relaxed R1CS based on the Nova paper. The R1CS structure is supported by a sparse matrix implementation found in the `sparse.rs` file, which provides fast matrix-vector multiplication similar to the SciPy library. This code constructs R1CS instances and witnesses, checks satisfiability, and provides functions to fold (relaxed) R1CS instances and witnesses, merging them into new relaxed R1CS instances and witnesses.

## R1CSShape

This struct includes the public input for a (relaxed) R1CS instance.

```rust
#[derive(Debug, Clone, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct R1CSShape<G: CurveGroup> {
    pub num_constraints: usize,
    pub num_vars: usize,
    pub num_io: usize,
    pub A: SparseMatrix<G::ScalarField>,
    pub B: SparseMatrix<G::ScalarField>,
    pub C: SparseMatrix<G::ScalarField>,
}
```

### Important Functions

```rust
impl<G: CurveGroup> R1CSShape<G> {
    fn validate(
        num_constraints: usize,
        num_vars: usize,
        num_io: usize,
        M: MatrixRef<'_, G::ScalarField>,
    ) -> Result<(), Error>;

    pub fn new(
        num_constraints: usize,
        num_vars: usize,
        num_io: usize,
        A: MatrixRef<'_, G::ScalarField>,
        B: MatrixRef<'_, G::ScalarField>,
        C: MatrixRef<'_, G::ScalarField>,
    ) -> Result<R1CSShape<G>, Error>;

    pub fn is_satisfied<C: CommitmentScheme<G>>(
        &self,
        U: &R1CSInstance<G, C>,
        W: &R1CSWitness<G>,
        pp: &C::PP,
    ) -> Result<(), Error>;

    pub fn is_relaxed_satisfied<C: CommitmentScheme<G>>(
        &self,
        U: &RelaxedR1CSInstance<G, C>,
        W: &RelaxedR1CSWitness<G>,
        pp: &C::PP,
    ) -> Result<(), Error>;
}
```

## R1CSWitness

```rust
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct R1CSWitness<G: CurveGroup> {
    pub W: Vec<G::ScalarField>,
}
```

### Important Functions

```rust
impl<G: CurveGroup> R1CSWitness<G> {
    pub fn new(shape: &R1CSShape<G>, W: &[G::ScalarField]) -> Result<Self, Error>;

    pub fn zero(shape: &R1CSShape<G>) -> Self;

    pub fn commit<C: CommitmentScheme<G>>(&self, pp: &C::PP) -> C::Commitment;
}
```

## R1CSInstance

```rust
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct R1CSInstance<G: CurveGroup, C: CommitmentScheme<G>> {
    pub commitment_W: C::Commitment,
    pub X: Vec<G::ScalarField>,
}
```

### Important Functions

```rust
impl<G: CurveGroup, C: CommitmentScheme<G>> R1CSInstance<G, C> {
    pub fn new(
        shape: &R1CSShape<G>,
        commitment_W: &C::Commitment,
        X: &[G::ScalarField],
    ) -> Result<Self, Error>;
}
```

## RelaxedR1CSWitness

```rust
#[derive(Default, Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct RelaxedR1CSWitness<G: CurveGroup> {
    pub W: Vec<G::ScalarField>,
    pub E: Vec<G::ScalarField>,
}
```

### Important Functions

```rust
impl<G: CurveGroup> RelaxedR1CSWitness<G> {
    pub fn zero(shape: &R1CSShape<G>) -> Self;

    pub fn from_r1cs_witness(shape: &R1CSShape<G>, witness: &R1CSWitness<G>) -> Self;

    pub fn fold(
        &self,
        W2: &R1CSWitness<G>,
        T: &[G::ScalarField],
        r: &G::ScalarField,
    ) -> Result<Self, Error>;

    pub fn fold_with_relaxed(
        &self,
        W2: &RelaxedR1CSWitness<G>,
        T: &[G::ScalarField],
        r: &G::ScalarField,
    ) -> Result<Self, Error>;
}
```

## RelaxedR1CSInstance

```rust
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct RelaxedR1CSInstance<G: CurveGroup, C: CommitmentScheme<G>> {
    pub commitment_W: C::Commitment,
    pub commitment_E: C::Commitment,
    pub X: Vec<G::ScalarField>,
}
```

### Important Functions

```rust
impl<G: CurveGroup, C: CommitmentScheme<G>> RelaxedR1CSInstance<G, C> {
    pub fn new(shape: &R1CSShape<G>) -> Self;

    pub fn from_r1cs_instance(
        shape: &R1CSShape<G>,
        instance: &R1CSInstance<G, C>,
    ) -> Result<Self, Error>;

    pub fn fold(
        &self,
        U2: &R1CSInstance<G, C>,
        comm_T: &C::Commitment,
        r: &G::ScalarField,
    ) -> Result<Self, Error>;

    pub fn fold_with_relaxed(
        &self,
        U2: &RelaxedR1CSInstance<G, C>,
        comm_T: &C::Commitment,
        r: &G::ScalarField,
    ) -> Result<Self, Error>;
}
```

## commit_T Function

This function computes the cross-term error `T` and its commitment, which is calculated by the prover and provided to the verifier as folding proof to aid the verifier in computing the folded instance from two (relaxed) R1CS instances.

```rust
pub fn commit_T<G: CurveGroup, C: CommitmentScheme<G>>(
    shape: &R1CSShape<G>,
    pp: &C::PP,
    U1: &RelaxedR1CSInstance<G, C>,
    W1: &RelaxedR1CSWitness<G>,
    U2: &R1CSInstance<G, C>,
    W2: &R1CSWitness<G>,
) -> Result<C::Commitment, Error>;
```

## commit_T_with_relaxed Function

Similar to `commit_T`, this function computes the cross-term error `T` and its commitment, but it is used when both instances are relaxed.

```rust
pub fn commit_T_with_relaxed<G: CurveGroup, C: CommitmentScheme<G>>(
    shape: &R1CSShape<G>,
    pp: &C::PP,
    U1: &RelaxedR1CSInstance<G, C>,
    W1: &RelaxedR1CSWitness<G>,
    U2: &RelaxedR1CSInstance<G, C>,
    W2: &RelaxedR1CSWitness<G>,
) -> Result<C::Commitment, Error>;
```

## Absorb Trait Implementation

The `Absorb` trait is implemented for the following structs in this code, including  `R1CSInstance`  and `RelaxedR1CSInstance`. The `Absorb` trait allows these structs to absorb scalar field elements, which is typically used in the context of Poseidon hashing. 

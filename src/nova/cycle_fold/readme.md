
# Coprocessor Circuit Implementation

## `SecondaryCircuit` Structure

The `SecondaryCircuit` struct is the core component of this implementation. It accepts the following inputs:

- **`g1`**: A point on the elliptic curve \( G \).
- **`g2`**: Another point on the elliptic curve \( G \).
- **`g_out`**: The expected result of the multiscalar multiplication on \( G \).
- **`r`**: A scalar field element used in the multiplication.
- **`flag`**: A boolean value that determines the operation mode.

Depending on the value of the `flag`, the circuit enforces one of two possible operations:

1. **If `flag` is `true`**:
   \[
   g_{\text{out}} = r \cdot g_1 + g_2
   \]

2. **If `flag` is `false`**:
   \[
   g_{\text{out}} = r \cdot (g_1 - g_2) + g_2 = r \cdot g_1 + (1 - r) \cdot g_2
   \]

The circuit is optimized to minimize the number of expensive scalar multiplication operations. It achieves this by conditionally selecting between `g1` and `g1 - g2` before performing a single scalar multiplication and adding `g2`.

The instance size, when translated to Rank-1 Constraint System (R1CS) form, includes 12 elements: 3 elliptic curve points (in projective form, i.e., 9 coordinates), 1 scalar field element, 1 boolean flag, and a leading 1 (a convention for R1CS inputs).

### `SecondaryCircuit` Struct

```rust
pub struct SecondaryCircuit<G1: SWCurveConfig> {
    pub(crate) g1: Projective<G1>,
    pub(crate) g2: Projective<G1>,
    pub(crate) g_out: Projective<G1>,
    pub(crate) r: G1::BaseField,
    pub(crate) flag: bool,
}
```

### Function `generate_constraints`

Generates the constraints needed for the circuit based on the provided inputs.

```rust
fn generate_constraints(
    self,
    cs: ConstraintSystemRef<G1::BaseField>,
) -> Result<(), SynthesisError>
```

### Function `setup_shape`

Sets up the R1CS shape for the secondary circuit, defining it over the base field of the secondary curve.

```rust
pub fn setup_shape<G1, G2>() -> Result<R1CSShape<G2>, SynthesisError>
```

### Function `synthesize`

Synthesizes the public input and witness for the circuit.

```rust
pub fn synthesize<G1, G2, C2>(
    circuit: SecondaryCircuit<G1>,
    pp_secondary: &C2::PP,
) -> Result<(R1CSInstance<G2, C2>, R1CSWitness<G2>), SynthesisError>
```

## `coprocessor_constraints.rs`

This file extends the functionality of the `SecondaryCircuit` by providing mechanisms to fold the secondary circuit's R1CS instances onto the primary curve. The main components in this file are `R1CSInstanceVar` and `RelaxedR1CSInstanceVar`.

### `R1CSInstanceVar` Struct

The `R1CSInstanceVar` struct represents the commitments (group points) and non-native field elements for the secondary circuit's R1CS instance. This struct is used for efficient scalar multiplication on the main curve.

```rust
/// Struct native on G2::BaseField = G1::ScalarField, hence native on primary curve
pub struct R1CSInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
{
    /// Commitment to witness.
    pub commitment_W: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    /// Public input of non-relaxed instance.
    pub X: Vec<NonNativeFieldVar<G2::ScalarField, G2::BaseField>>,

    _commitment_scheme: PhantomData<C2>,
}
```

### `RelaxedR1CSInstanceVar` Struct

The `RelaxedR1CSInstanceVar` struct functions as the running instance on the main curve, and in each cycle, it folds the `R1CSInstanceVar` with itself. Note this folding operation requires scalar multiplication on the native field which is cheap and linear combination of non-native field elements which is expensive. In fact, it requires 12 non-native linear combination since the instance size of the R1CS instance for the secondary curve is 12, as previously pointed out.

```rust
pub struct RelaxedR1CSInstanceVar<G2, C2>
where
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
{
    /// Commitment to witness.
    pub commitment_W: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    /// Commitment to error vector.
    pub commitment_E: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    /// Public input of relaxed instance.
    pub X: Vec<NonNativeFieldVar<G2::ScalarField, G2::BaseField>>,

    _commitment_scheme: PhantomData<C2>,
}
```


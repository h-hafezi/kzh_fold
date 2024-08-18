
# Non Native Group/Field

## `util.rs`

This file contains utility functions for converting between field elements in different prime fields, namely:

- **`convert_field_one_to_field_two<Fr, Fq>`**: Converts an element of the field `Fr` into an element of the field `Fq` by taking the big integer values. The conversion is straightforward but note that if the size of `Fr` is greater than `Fq`, there is a negligible probability of overflow (approximately \(2^{-128}\)) when using curves like Pasta or BN254 with Grumpkin.

- **`non_native_to_fpvar<ScalarField, BaseField>`**: Converts a `NonNativeFieldVar<BaseField, ScalarField>` into an `FpVar<ScalarField>` within a constraint system. This is done by transforming the non-native field variable into a boolean vector and then parsing it into an `FpVar`. The function has a constraint count of **1078**.

#### Function Signatures

```rust
pub fn convert_field_one_to_field_two<Fr, Fq>(first_field: Fr) -> Fq
where
    Fr: PrimeField,
    Fq: PrimeField,
{
    
}

pub fn non_native_to_fpvar<ScalarField, BaseField>(
    non_native_var: &NonNativeFieldVar<BaseField, ScalarField>,
) -> FpVar<ScalarField>
where
    ScalarField: PrimeField,
    BaseField: PrimeField,
{
    
}
```

## `non_native_affine_var.rs`

This code has been partially borrowed from [Nexus](https://github.com/nexus-xyz/nexus-zkvm). This file implements an affine variable for group points where the base is in a non-native field, meaning that the constraint system operates on the scalar field `Fr` while the group points (e.g., projective coordinates) are defined over the base field `Fq`. The `NonNativeAffineVar` struct is not meant for arithmetic operations but for representing group points, hashing, and enforcing equality.

Initialization of a `NonNativeAffineVar` requires two non-native field values and a boolean value. The boolean value takes **543** constraints, and transforming it into a vector of field elements for hashing takes **2760** constraints.

#### Struct Definition

```rust
#[derive(Debug)]
pub struct NonNativeAffineVar<G1>
where
    G1: SWCurveConfig,
    G1::BaseField: PrimeField,
{
    pub x: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,
    pub y: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,
    pub infinity: Boolean<G1::ScalarField>,
}
```

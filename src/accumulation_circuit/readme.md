# AccumulatorVerifier

## Overview

The `AccumulatorVerifier` and its corresponding variable struct `AccumulatorVerifierVar` are designed to perform cryptographic accumulation and verification in zk-SNARK circuits. 

### AccumulatorVerifier
The accumulator circuit in fact checks that one accumulator instance is result of the two given instances. More concretely apart from the two instances, the randomness and the accumulation proof `Q` and the result accumulation instance, it takes some auxiliary inputs which are instance of CycleFold to help it with multi scalar operations and as a result since these instances will be folded with a running cycle fold instances, it will also require to takes some `ProjectiveVar` e.g. `com_C_var`, `com_T_var`, `com_E_1_var` and `com_E_2_var` as folding proof for CycleFold. The struct looks as it follows:
```rust
#[derive(Clone)]
pub struct AccumulatorVerifierVar<G1, G2, C2>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
{
    /// auxiliary input which helps to have C'' = (1-beta) * C + beta * C' without scalar multiplication
    pub auxiliary_input_C_var: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have T'' = (1-beta) * T + beta * T' without scalar multiplication
    pub auxiliary_input_T_var: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E_{temp} = (1-beta) * E + beta * E' without scalar multiplication
    pub auxiliary_input_E_1_var: R1CSInstanceVar<G2, C2>,
    /// auxiliary input which helps to have E'' = E_{temp} + beta * (1-beta) * Q without scalar multiplication
    pub auxiliary_input_E_2_var: R1CSInstanceVar<G2, C2>,

    /// the randomness used for taking linear combination and its non-native counterpart
    pub beta_var: FpVar<G1::ScalarField>,
    pub beta_var_non_native: NonNativeFieldVar<G1::BaseField, G1::ScalarField>,

    /// accumulation proof
    pub Q_var: NonNativeAffineVar<G1>,

    /// accumulation proof for cycle fold (this is also the order of accumulating with cycle_fold_running_instance)
    pub com_C_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub com_T_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub com_E_1_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,
    pub com_E_2_var: ProjectiveVar<G2, FpVar<G2::BaseField>>,

    pub current_accumulator_instance_var: AccumulatorInstanceVar<G1>,
    pub running_accumulator_instance_var: AccumulatorInstanceVar<G1>,
    pub final_accumulator_instance_var: AccumulatorInstanceVar<G1>,

    pub running_cycle_fold_instance_var: RelaxedR1CSInstanceVar<G2, C2>,
    pub final_cycle_fold_instance_var: RelaxedR1CSInstanceVar<G2, C2>,

    // these are constant values
    pub n: u32,
    pub m: u32,
}

```
### `accumulate` Function

The `accumulate` function checks one verifier_circuit is indeed valid, by checking the following conditions:

#### 1. Checking Consistency Between `beta` and `beta_non_native`

```rust
let beta_bits = self.beta_var_non_native.to_bits_le().unwrap();
let beta_ = Boolean::le_bits_to_fp_var(beta_bits.as_slice()).unwrap();
self.beta_var.enforce_equal(&beta_).expect("error while enforcing equality");
```
- **Operation**: Verifies that the native and non-native representations of `beta` are consistent.

#### 2. Poseidon Hash and Consistency Check with `beta`

```rust
// compute Poseidon hash and make sure it's consistent with input beta
let mut hash_object = PoseidonHashVar::new(self.current_accumulator_instance_var.cs());
let mut sponge = Vec::new();
sponge.extend(self.current_accumulator_instance_var.to_sponge_field_elements().unwrap());
sponge.extend(self.running_accumulator_instance_var.to_sponge_field_elements().unwrap());
sponge.extend(self.Q_var.to_sponge_field_elements().unwrap());
hash_object.update_sponge(sponge);
hash_object.output().enforce_equal(&self.beta_var).expect("error while enforcing equality");
```
- **Operation**: Computes the Poseidon hash of the accumulator inputs and `Q` and ensures it matches the `beta` value.

#### 3. Linear Combination of Commitment `C`

```rust
let (flag, r, g1, g2, C_var) = self.auxiliary_input_C_var.parse_secondary_io::<G1>().unwrap();
// g1 == acc.C
self.running_accumulator_instance_var.C_var.enforce_equal(&g1).expect("error while enforcing equality");
// g2 == instance.C
self.current_accumulator_instance_var.C_var.enforce_equal(&g2).expect("error while enforcing equality");
// enforce flag to be false
flag.enforce_equal(&NonNativeFieldVar::zero()).expect("error while enforcing equality");
// check r to be equal to beta
r.enforce_equal(&self.beta_var_non_native).expect("error while enforcing equality");
// check out the result C_var is consistent with result_acc
C_var.enforce_equal(&self.final_accumulator_instance_var.C_var).expect("error while enforcing equality");
```
- **Operation**:  Since `C'' = (1-beta) * C + beta * C'`, then to compute this with cycle fold, we use a secondary circuit with the following parameters:
```rust
auxiliary_input_C = SecondaryCircuit {
    g1: C',
    g2: C,
    g_out: C'',
    r: beta,
    flag: false,
}
```
This operation makes sure that the secondary circuit instance is indeed consistent with the input accumulator commitment `C`, `C'`, the output accumulator commitment `C''`, the randomness `beta` and its flag parameter is `false`.
#### 4. Linear Combination of Commitment `T`

- **Operation**: Similar to step 3 but applied to the commitment `T`.

#### 5. Compute New Error Term
Despite `C''` and `T''`, the error vector require two multi scalar operation since `E''= (1-beta) * E + beta * E' + beta * (1-beta) * Q`, we compute this by diving it into steps:
 * `E_temp = (1-beta) * E + beta * E'`
 * `E'' = E_temp + beta * (1-beta) * Q'`
Which means we need to secondary circuits with the following parameters:
```rust
auxiliary_input_E_1 = SecondaryCircuit {
    g1: E',
    g2: E,
    g_out: E_temp,
    r: beta,
    flag: false,
}

auxiliary_input_C = SecondaryCircuit {
    g1: Q,
    g2: E_temp,
    g_out: E,
    r: beta * (1-beta),
    flag: true,
}
```

- **Operation**: Similar to step 3, this steps also checks the consistencies w.r.t circuits above.


#### 6. Native Field Operations: Linear Combinations

- **Operation**: Performs linear combinations for the fields `b`, `c`, `y`, `z_b`, and `z_c`, ensuring they match the final accumulator's values.

#### 7. Equality Assertions for `z_b` and `z_c`
If the current instance isn't a relaxed once (`E=0`) then we need to check that `z_b=b^n-1` and `z_c=c^m-1`.


#### 8. Folding of Cycle Instances


- **Operation**: Fold the cycle fold instances with its running instance and check the result to be equal to `final_cycle_fold_instance_var`.


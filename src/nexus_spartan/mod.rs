extern crate core;
extern crate digest;
extern crate rand;
extern crate sha3;

mod commitments;
pub mod committed_relaxed_snark;
pub mod crr1cs;
pub mod crr1csproof;
pub mod errors;
pub mod commitment_traits;
mod product_tree;
pub mod r1csinstance;
pub mod sparse_polynomial;
pub mod sumcheck;
mod timer;
mod unipoly;
mod conversion;
pub mod sumcheck_circuit;
pub mod partial_verifier;
mod sparse_mlpoly;
pub mod matrix_evaluation_accumulation;

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::{BigInteger, PrimeField};
use ark_serialize::*;
use core::cmp::max;
use errors::R1CSError;
use r1csinstance::{
    R1CSCommitment, R1CSDecommitment, R1CSInstance,
};
use crate::kzh::KZH;

/// `ComputationCommitment` holds a public preprocessed NP statement (e.g., R1CS)
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct ComputationCommitment<E: Pairing, PC: KZH<E>>
where
    <E as Pairing>::ScalarField: Absorb,
{
    comm: R1CSCommitment<E, PC>,
}

/// `ComputationDecommitment` holds information to decommit `ComputationCommitment`
#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct ComputationDecommitment<F>
where
    F: Sync + CanonicalDeserialize + CanonicalSerialize + PrimeField + Absorb,
{
    decomm: R1CSDecommitment<F>,
}

/// `Assignment` holds an assignment of values to either the inputs or variables in an `Instance`
#[derive(Clone)]
pub struct Assignment<F> {
    pub assignment: Vec<F>,
}

impl<F: PrimeField + Absorb> Assignment<F> {
    /// Constructs a new `Assignment` from a vector
    pub fn new(assignment: &[F]) -> Result<Self, R1CSError> {
        let bytes_to_scalar = |vec: &[F]| -> Result<Vec<F>, R1CSError> {
            let mut vec_scalar: Vec<F> = Vec::new();
            for v in vec {
                vec_scalar.push(*v);
            }
            Ok(vec_scalar)
        };

        let assignment_scalar = bytes_to_scalar(assignment)?;

        Ok(Assignment {
            assignment: assignment_scalar,
        })
    }

    /// pads Assignment to the specified length
    fn pad(&self, len: usize) -> VarsAssignment<F> {
        // check that the new length is higher than current length
        assert!(len > self.assignment.len());

        let padded_assignment = {
            let mut padded_assignment = self.assignment.clone();
            padded_assignment.extend(vec![F::zero(); len - self.assignment.len()]);
            padded_assignment
        };

        VarsAssignment {
            assignment: padded_assignment,
        }
    }
}

/// `VarsAssignment` holds an assignment of values to variables in an `Instance`
pub type VarsAssignment<F> = Assignment<F>;

/// `InputsAssignment` holds an assignment of values to variables in an `Instance`
pub type InputsAssignment<F> = Assignment<F>;

/// `Instance` holds the description of R1CS matrices
#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct Instance<F: PrimeField + Absorb> {
    pub inst: R1CSInstance<F>,
}

impl<F: PrimeField + Absorb> Instance<F> {
    /// Constructs a new `Instance` and an associated satisfying assignment
    pub fn new(
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
        A: &[(usize, usize, F)],
        B: &[(usize, usize, F)],
        C: &[(usize, usize, F)],
    ) -> Result<Self, R1CSError> {
        // padding
        let (num_vars_padded, num_cons_padded) = {
            let num_vars_padded = {
                let mut num_vars_padded = num_vars;

                // ensure that num_inputs + 1 <= num_vars
                num_vars_padded = max(num_vars_padded, num_inputs + 1);

                // ensure that num_vars_padded a power of two
                if num_vars_padded.next_power_of_two() != num_vars_padded {
                    num_vars_padded = num_vars_padded.next_power_of_two();
                }
                num_vars_padded
            };

            let num_cons_padded = {
                let mut num_cons_padded = num_cons;

                // ensure that num_cons_padded is at least 2
                if num_cons_padded == 0 || num_cons_padded == 1 {
                    num_cons_padded = 2;
                }

                // ensure that num_cons_padded is power of 2
                if num_cons.next_power_of_two() != num_cons {
                    num_cons_padded = num_cons.next_power_of_two();
                }
                num_cons_padded
            };

            (num_vars_padded, num_cons_padded)
        };

        let bytes_to_scalar =
            |tups: &[(usize, usize, F)]| -> Result<Vec<(usize, usize, F)>, R1CSError> {
                let mut mat: Vec<(usize, usize, F)> = Vec::new();
                for &(row, col, val) in tups {
                    // row must be smaller than num_cons
                    if row >= num_cons {
                        return Err(R1CSError::InvalidIndex);
                    }

                    // col must be smaller than num_vars + 1 + num_inputs
                    if col >= num_vars + 1 + num_inputs {
                        return Err(R1CSError::InvalidIndex);
                    }

                    if col >= num_vars {
                        // Column refers to inputs (not witness). Move it past padding
                        mat.push((row, col + num_vars_padded - num_vars, val));
                    } else {
                        // Column refers to witness
                        mat.push((row, col, val));
                    }
                }

                // pad with additional constraints up until num_cons_padded if the original constraints were 0 or 1
                // we do not need to pad otherwise because the dummy constraints are implicit in the sum-check protocol
                if num_cons == 0 || num_cons == 1 {
                    for i in tups.len()..num_cons_padded {
                        mat.push((i, num_vars, F::zero()));
                    }
                }

                Ok(mat)
            };

        let A_scalar = bytes_to_scalar(A);
        if A_scalar.is_err() {
            return Err(A_scalar.err().unwrap());
        }

        let B_scalar = bytes_to_scalar(B);
        if B_scalar.is_err() {
            return Err(B_scalar.err().unwrap());
        }

        let C_scalar = bytes_to_scalar(C);
        if C_scalar.is_err() {
            return Err(C_scalar.err().unwrap());
        }

        let inst = R1CSInstance::<F>::new(
            num_cons_padded,
            num_vars_padded,
            num_inputs,
            &A_scalar.unwrap(),
            &B_scalar.unwrap(),
            &C_scalar.unwrap(),
        );

        Ok(Instance { inst })
    }

    /// Checks if a given R1CSInstance is satisfiable with a given variables and inputs assignments
    pub fn is_sat(
        &self,
        vars: &VarsAssignment<F>,
        inputs: &InputsAssignment<F>,
    ) -> Result<bool, R1CSError> {
        if vars.assignment.len() > self.inst.get_num_vars() {
            return Err(R1CSError::InvalidNumberOfInputs);
        }

        if inputs.assignment.len() != self.inst.get_num_inputs() {
            return Err(R1CSError::InvalidNumberOfInputs);
        }

        // we might need to pad variables
        let padded_vars = {
            let num_padded_vars = self.inst.get_num_vars();
            let num_vars = vars.assignment.len();
            if num_padded_vars > num_vars {
                vars.pad(num_padded_vars)
            } else {
                vars.clone()
            }
        };

        Ok(
            self
                .inst
                .is_sat(&padded_vars.assignment, &inputs.assignment),
        )
    }

    /// Constructs a new synthetic R1CS `Instance` and an associated satisfying assignment
    pub fn produce_synthetic_r1cs(
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
    ) -> (Instance<F>, VarsAssignment<F>, InputsAssignment<F>) {
        let (inst, vars, inputs) = R1CSInstance::produce_synthetic_r1cs(num_cons, num_vars, num_inputs);
        (
            Instance { inst },
            VarsAssignment { assignment: vars },
            InputsAssignment { assignment: inputs },
        )
    }
}

pub fn analyze_vector_sparseness<F: PrimeField>(name: &str, z: &Vec<F>) {
    // Count the number of zeroes and ones in the witness
    let total = z.len();
    let zero_count = z.iter().filter(|&x| x.is_zero()).count();
    let one_count = z.iter().filter(|&x| x.is_one()).count();

    // Calculate relative counts
    let zero_relative = zero_count as f64 / total as f64;
    let one_relative = one_count as f64 / total as f64;

    // Compute the number of bits for each element
    let num_bits: Vec<u32> = z.iter().map(|x| x.into_bigint().num_bits()).collect();

    // Compute average, min, and max num_bits
    let sum_bits: u32 = num_bits.iter().sum();
    let avg_bits = sum_bits as f64 / total as f64;
    let min_bits = num_bits.iter().min().unwrap();
    let max_bits = num_bits.iter().max().unwrap();

    // Get the modulus bit length (maximum bit length per element)
    let modulus_bit_length = F::MODULUS_BIT_SIZE as u32;

    // Compute the maximum possible total bit length
    let max_bit_length_possible = total as u32 * modulus_bit_length;

    // Compute the effective bit length percentage
    let effective_bit_length_percentage = (sum_bits as f64 / max_bit_length_possible as f64) * 100.0;

    // Output the results
    println!("[*] Analyzing {}", name);
    println!("\tsize: {}", total);
    println!("\tzero count: {}", zero_count);
    println!("\tzero percentage: {:.2}%", zero_relative * 100.0);
    println!("\tone count: {}", one_count);
    println!("\tone percentage: {:.2}%", one_relative * 100.0);
    println!("\taverage bit length: {:.2}", avg_bits);
    println!("\tmin bit length: {}", min_bits);
    println!("\tmax bit length: {}", max_bits);
    println!("\teffective bit length (total bits): {}", sum_bits);
    println!("\tmax possible bit length (total bits): {}", max_bit_length_possible);
    println!(
        "\teffective bit length percentage: {:.2}%",
        effective_bit_length_percentage
    );
}

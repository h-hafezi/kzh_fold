#![allow(dead_code)]
/// This is mostly a copy of the SNARK implementation in lib.rs, with minor modifications to work with committed relaxed R1CS.
use core::cmp::max;

use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use super::{crr1csproof::{CRR1CSKey}, polycommitments::PolyCommitmentScheme, r1csinstance::{R1CSCommitmentGens}};

/// `SNARKGens` holds public parameters for producing and verifying proofs with the Spartan SNARK
#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct CRSNARKKey<E: Pairing, PC: PolyCommitmentScheme<E>> {
    pub gens_r1cs_sat: CRR1CSKey<E, PC>,
    pub gens_r1cs_eval: R1CSCommitmentGens<E, PC>,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> CRSNARKKey<E, PC> {
    /// Constructs a new `SNARKGens` given the size of the R1CS statement
    /// `num_nz_entries` specifies the maximum number of non-zero entries in any of the three R1CS matrices
    pub fn new(
        SRS: &PC::SRS,
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
        num_nz_entries: usize,
    ) -> Self {
        let num_vars_padded = Self::get_num_vars_padded(num_vars, num_inputs);
        let gens_r1cs_sat = CRR1CSKey::<E, PC>::new(SRS, num_cons, num_vars_padded);
        let gens_r1cs_eval =
            R1CSCommitmentGens::new(SRS, num_cons, num_vars_padded, num_inputs, num_nz_entries);
        CRSNARKKey {
            gens_r1cs_sat,
            gens_r1cs_eval,
        }
    }
    fn get_num_vars_padded(num_vars: usize, num_inputs: usize) -> usize {
        let mut num_vars_padded = max(num_vars, num_inputs + 1);
        if num_vars_padded != num_vars_padded.next_power_of_two() {
            num_vars_padded = num_vars_padded.next_power_of_two();
        }
        num_vars_padded
    }
    pub fn get_min_num_vars(
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
        num_nz_entries: usize,
    ) -> usize {
        let num_vars_padded = Self::get_num_vars_padded(num_vars, num_inputs);
        let min_num_vars_sat = CRR1CSKey::<E, PC>::get_min_num_vars(num_cons, num_vars_padded);
        let min_num_vars_eval =
            R1CSCommitmentGens::<E, PC>::get_min_num_vars(num_cons, num_vars_padded, num_nz_entries);
        max(min_num_vars_sat, min_num_vars_eval)
    }
}


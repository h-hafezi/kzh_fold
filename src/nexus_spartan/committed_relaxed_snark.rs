#![allow(dead_code)]
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
/// This is mostly a copy of the SNARK implementation in lib.rs, with minor modifications to work with committed relaxed R1CS.
use core::cmp::max;
use rand::thread_rng;
use crate::kzh::KZH;
use crate::math::Math;

/// `SNARKGens` holds public parameters for producing and verifying proofs with the Spartan SNARK
#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct CRSNARKKey<E: Pairing, PC: KZH<E>>
where
    <E as Pairing>::ScalarField: Absorb,
{
    pub gens_r1cs_sat: PC::SRS,
}

impl<E: Pairing, PC: KZH<E>> CRSNARKKey<E, PC>
where
    <E as Pairing>::ScalarField: Absorb,
{
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
        let gens_r1cs_sat = PC::setup(
            max(num_cons, num_vars).log_2(),
            &mut thread_rng(),
        );

        CRSNARKKey {
            gens_r1cs_sat,
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
    ) -> usize {
        let num_vars_padded = Self::get_num_vars_padded(num_vars, num_inputs);
        let n = max(num_cons, num_vars_padded);
        n.log_2()
    }
}


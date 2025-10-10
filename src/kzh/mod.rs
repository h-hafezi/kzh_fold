use std::fmt::Debug;
use ark_crypto_primitives::sponge::Absorb;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use rand::prelude::IteratorRandom;
use rand::Rng;
use crate::nexus_spartan::commitment_traits::ToAffine;
use crate::transcript::transcript::AppendToTranscript;

pub mod kzh2;

pub mod kzh3;

pub mod kzh4;
pub mod zk_kzh2;
mod zk_kzh3;

pub trait KZH<E: Pairing> where <E as Pairing>::ScalarField: Absorb {
    type Degree;
    type SRS: CanonicalSerialize + CanonicalDeserialize + Clone;
    type Commitment: AppendToTranscript<E::ScalarField>
    + Debug
    + CanonicalSerialize
    + CanonicalDeserialize
    + PartialEq
    + Eq
    + Clone
    + AppendToTranscript<E::ScalarField>
    + ToAffine<E>;
    type Aux: AppendToTranscript<E::ScalarField> + Sync + CanonicalSerialize + CanonicalDeserialize + Debug;

    type Opening: Sync + CanonicalSerialize + CanonicalDeserialize + Debug;

    fn split_input<T: Clone>(srs: &Self::SRS, input: &[T], default: T) -> Vec<Vec<T>>;

    fn get_degree_from_maximum_supported_degree(n: usize) -> Self::Degree;

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS;

    fn commit<R: Rng>(
        srs: &Self::SRS,
        poly: &MultilinearPolynomial<E::ScalarField>,
        rng: &mut R
    ) -> (Self::Commitment, Self::Aux);

    fn open<R: Rng>(
        srs: &Self::SRS,
        input: &[E::ScalarField],
        com: &Self::Commitment,
        aux: &Self::Aux,
        poly: &MultilinearPolynomial<E::ScalarField>,
        rng: &mut R
    ) -> Self::Opening;

    fn verify(
        srs: &Self::SRS,
        input: &[E::ScalarField],
        output: &E::ScalarField,
        com: &Self::Commitment,
        open: &Self::Opening,
    );
}

#[derive(Debug, Clone, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
pub struct SparseMultilinearPolynomial<F: PrimeField> {
    /// the number of variables in the multilinear polynomial
    pub num_variables: usize,
    /// non-zero evaluations stored as (index, value) pairs
    /// index ranges from 0..2^num_variables
    pub non_zero_entries: Vec<(usize, F)>,
}

impl<F: PrimeField> SparseMultilinearPolynomial<F> {
    pub fn new(num_variables: usize, non_zero_entries: Vec<(usize, F)>) -> Self {
        // optional: sanity check indices
        let max_index = 1 << num_variables;
        assert!(
            non_zero_entries.iter().all(|(i, _)| *i < max_index),
            "Index out of range for given num_variables"
        );
        Self { num_variables, non_zero_entries }
    }

    /// Convert sparse polynomial into dense multilinear polynomial
    pub fn to_dense(&self) -> MultilinearPolynomial<F> {
        let len = 1 << self.num_variables;
        let mut evals = vec![F::zero(); len];
        for (i, val) in &self.non_zero_entries {
            evals[*i] = *val;
        }
        MultilinearPolynomial {
            num_variables: self.num_variables,
            evaluation_over_boolean_hypercube: evals.clone(),
            len,
        }
    }

    /// Generate a random sparse multilinear polynomial with `num_non_zero` entries
    pub fn rand<R: Rng>(
        rng: &mut R,
        num_variables: usize,
        num_non_zero: usize,
    ) -> Self {
        let len = 1 << num_variables;
        assert!(
            num_non_zero < len,
            "Number of non-zero entries must be less than 2^num_variables"
        );

        // sample `num_non_zero` distinct indices
        let indices: Vec<usize> = (0..len)
            .choose_multiple(rng, num_non_zero);

        // assign random values to these indices
        let mut non_zero_entries = Vec::with_capacity(num_non_zero);
        for idx in indices {
            let val = F::rand(rng);
            non_zero_entries.push((idx, val));
        }

        Self {
            num_variables,
            non_zero_entries,
        }
    }
}


/// Prepads a slice with the field's zero element to a specified length.
///
/// If the input slice's length is already greater than or equal to `len`,
/// it is returned as a new `Vec` without modification. Otherwise, `F::zero()`
/// elements are prepended until the vector's length is `len`.
///
/// # Arguments
/// * `input`: The slice to pad.
/// * `len`: The target minimum length.
///
/// # Returns
/// A new `Vec<F>` with the padded elements.
pub fn pad_at_start<F: PrimeField>(input: &[F], len: usize) -> Vec<F> {
    let current_len = input.len();

    // If the input is already long enough, just return it as a Vec.
    if current_len >= len {
        return input.to_vec();
    }

    // Calculate how many zeros are needed for padding.
    let num_zeros = len - current_len;

    // 1. Create a new vector with the required number of zeros.
    let mut padded_vec = vec![F::zero(); num_zeros];

    // 2. Efficiently append the original input slice to the end of the zeros.
    padded_vec.extend_from_slice(input);

    padded_vec
}
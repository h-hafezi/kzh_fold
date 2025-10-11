#![allow(clippy::type_complexity)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

use crate::math::Math;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::*;
use ark_std::cmp::max;
use crate::kzh::KZH;

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct SparseMatEntry<F: PrimeField> {
    row: usize,
    col: usize,
    val: F,
}

impl<F: PrimeField> SparseMatEntry<F> {
    pub fn new(row: usize, col: usize, val: F) -> Self {
        SparseMatEntry { row, col, val }
    }
}

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct SparseMatPolynomial<F: PrimeField> {
    num_vars_x: usize,
    num_vars_y: usize,
    M: Vec<SparseMatEntry<F>>,
}

impl<F: PrimeField + Absorb> SparseMatPolynomial<F> {
    pub fn new(num_vars_x: usize, num_vars_y: usize, M: Vec<SparseMatEntry<F>>) -> Self {
        SparseMatPolynomial {
            num_vars_x,
            num_vars_y,
            M,
        }
    }

    fn evaluate_with_tables(&self, eval_table_rx: &[F], eval_table_ry: &[F]) -> F {
        assert_eq!(self.num_vars_x.pow2(), eval_table_rx.len());
        assert_eq!(self.num_vars_y.pow2(), eval_table_ry.len());

        (0..self.M.len())
            .map(|i| {
                let row = self.M[i].row;
                let col = self.M[i].col;
                let val = &self.M[i].val;
                eval_table_rx[row] * eval_table_ry[col] * val
            })
            .sum()
    }

    pub fn multi_evaluate(polys: &[&SparseMatPolynomial<F>], rx: &[F], ry: &[F]) -> Vec<F> {
        let eval_table_rx = EqPolynomial::new(rx.to_vec()).evals();
        let eval_table_ry = EqPolynomial::new(ry.to_vec()).evals();

        (0..polys.len())
            .map(|i| polys[i].evaluate_with_tables(&eval_table_rx, &eval_table_ry))
            .collect::<Vec<F>>()
    }

    pub fn multiply_vec(&self, num_rows: usize, num_cols: usize, z: &[F]) -> Vec<F> {
        assert_eq!(z.len(), num_cols);

        (0..self.M.len())
            .map(|i| {
                let row = self.M[i].row;
                let col = self.M[i].col;
                let val = self.M[i].val;
                (row, val * z[col])
            })
            .fold(vec![F::zero(); num_rows], |mut Mz, (r, v)| {
                Mz[r] += v;
                Mz
            })
    }

    pub fn compute_eval_table_sparse(&self, rx: &[F], num_rows: usize, num_cols: usize) -> Vec<F> {
        assert_eq!(rx.len(), num_rows);

        let mut M_evals: Vec<F> = vec![F::zero(); num_cols];

        for i in 0..self.M.len() {
            let entry = &self.M[i];
            M_evals[entry.col] += rx[entry.row] * entry.val;
        }
        M_evals
    }
}

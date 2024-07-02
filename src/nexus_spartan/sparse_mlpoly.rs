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

#[derive(CanonicalDeserialize, CanonicalSerialize)]
struct AddrTimestamps<F>
where
    F: Sync + CanonicalSerialize + CanonicalDeserialize + PrimeField + Absorb,
{
    ops_addr_usize: Vec<Vec<usize>>,
    ops_addr: Vec<MultilinearPolynomial<F>>,
    read_ts: Vec<MultilinearPolynomial<F>>,
    audit_ts: MultilinearPolynomial<F>,
}

#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct MultiSparseMatPolynomialAsDense<F>
where
    F: Sync + CanonicalSerialize + CanonicalDeserialize + PrimeField + Absorb,
{
    batch_size: usize,
    val: Vec<MultilinearPolynomial<F>>,
    row: AddrTimestamps<F>,
    col: AddrTimestamps<F>,
    comb_ops: MultilinearPolynomial<F>,
    comb_mem: MultilinearPolynomial<F>,
}

#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct SparseMatPolyCommitmentKey<E, PC>
where
    E: Pairing,
    PC: KZH<E>,
    E::ScalarField: Absorb,
{
    gens_ops: PC::SRS,
    gens_mem: PC::SRS,
    gens_derefs: PC::SRS,
}

impl<E: Pairing, PC: KZH<E>> SparseMatPolyCommitmentKey<E, PC>
where
    <E as Pairing>::ScalarField: Absorb,
{
    pub fn new(
        SRS: &PC::SRS,
        num_vars_x: usize,
        num_vars_y: usize,
        num_nz_entries: usize,
        batch_size: usize,
    ) -> SparseMatPolyCommitmentKey<E, PC> {
        let (_num_vars_ops, _num_vars_mem, _num_vars_derefs) =
            Self::get_gens_sizes(num_vars_x, num_vars_y, num_nz_entries, batch_size);

        let gens_ops = SRS.clone();
        let gens_mem = SRS.clone();
        let gens_derefs = SRS.clone();
        SparseMatPolyCommitmentKey {
            gens_ops,
            gens_mem,
            gens_derefs,
        }
    }
    fn get_gens_sizes(
        num_vars_x: usize,
        num_vars_y: usize,
        num_nz_entries: usize,
        batch_size: usize,
    ) -> (usize, usize, usize) {
        let num_vars_ops =
            num_nz_entries.next_power_of_two().log_2() + (batch_size * 5).next_power_of_two().log_2();
        let num_vars_mem = if num_vars_x > num_vars_y {
            num_vars_x
        } else {
            num_vars_y
        } + 1;
        let num_vars_derefs =
            num_nz_entries.next_power_of_two().log_2() + (batch_size * 2).next_power_of_two().log_2();
        (num_vars_ops, num_vars_mem, num_vars_derefs)
    }
    pub fn get_min_num_vars(
        num_vars_x: usize,
        num_vars_y: usize,
        num_nz_entries: usize,
        batch_size: usize,
    ) -> usize {
        let (num_vars_ops, num_vars_mem, num_vars_derefs) =
            Self::get_gens_sizes(num_vars_x, num_vars_y, num_nz_entries, batch_size);
        max(num_vars_ops, max(num_vars_mem, num_vars_derefs))
    }
}

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct SparseMatPolyCommitment<E: Pairing, PC: KZH<E>>
where
    <E as Pairing>::ScalarField: Absorb,
{
    batch_size: usize,
    num_ops: usize,
    num_mem_cells: usize,
    comm_comb_ops: PC::Commitment,
    comm_comb_mem: PC::Commitment,
}

impl<E: Pairing, PC: KZH<E>> AppendToTranscript<E::ScalarField>
for SparseMatPolyCommitment<E, PC>
where
    <E as Pairing>::ScalarField: Absorb,
{
    fn append_to_transcript(&self, _label: &'static [u8], transcript: &mut Transcript<E::ScalarField>)
    where
        <E as Pairing>::ScalarField: Absorb,
    {
        transcript.append_u64(b"batch_size", self.batch_size as u64);
        transcript.append_u64(b"num_ops", self.num_ops as u64);
        transcript.append_u64(b"num_mem_cells", self.num_mem_cells as u64);
        self
            .comm_comb_ops
            .append_to_transcript(b"comm_comb_ops", transcript);
        self
            .comm_comb_mem
            .append_to_transcript(b"comm_comb_mem", transcript);
    }
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

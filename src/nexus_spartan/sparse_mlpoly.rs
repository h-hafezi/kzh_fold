#![allow(clippy::type_complexity)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

use crate::math::Math;
use crate::nexus_spartan::polycommitments::{PCSKeys, PolyCommitmentScheme};
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::*;
use ark_std::cmp::max;

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

pub struct Derefs<F>
where
    F: Sync + CanonicalDeserialize + CanonicalSerialize + PrimeField,
{
    row_ops_val: Vec<MultilinearPolynomial<F>>,
    col_ops_val: Vec<MultilinearPolynomial<F>>,
    comb: MultilinearPolynomial<F>,
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

impl<F: PrimeField + Absorb> AddrTimestamps<F> {
    pub fn new(num_cells: usize, num_ops: usize, ops_addr: Vec<Vec<usize>>) -> Self {
        for item in ops_addr.iter() {
            assert_eq!(item.len(), num_ops);
        }

        let mut audit_ts = vec![0usize; num_cells];
        let mut ops_addr_vec: Vec<MultilinearPolynomial<F>> = Vec::new();
        let mut read_ts_vec: Vec<MultilinearPolynomial<F>> = Vec::new();
        for ops_addr_inst in ops_addr.iter() {
            let mut read_ts = vec![0usize; num_ops];

            // since read timestamps are trustworthy, we can simply increment the r-ts to obtain a w-ts
            // this is sufficient to ensure that the write-set, consisting of (addr, val, ts) tuples, is a set
            for i in 0..num_ops {
                let addr = ops_addr_inst[i];
                assert!(addr < num_cells);
                let r_ts = audit_ts[addr];
                read_ts[i] = r_ts;

                let w_ts = r_ts + 1;
                audit_ts[addr] = w_ts;
            }

            ops_addr_vec.push(MultilinearPolynomial::from_usize(ops_addr_inst));
            read_ts_vec.push(MultilinearPolynomial::from_usize(&read_ts));
        }

        AddrTimestamps {
            ops_addr: ops_addr_vec,
            ops_addr_usize: ops_addr,
            read_ts: read_ts_vec,
            audit_ts: MultilinearPolynomial::from_usize(&audit_ts),
        }
    }

    fn deref_mem(addr: &[usize], mem_val: &[F]) -> MultilinearPolynomial<F> {
        MultilinearPolynomial::new(
            (0..addr.len())
                .map(|i| {
                    let a = addr[i];
                    mem_val[a]
                })
                .collect::<Vec<F>>(),
        )
    }

    pub fn deref(&self, mem_val: &[F]) -> Vec<MultilinearPolynomial<F>> {
        (0..self.ops_addr.len())
            .map(|i| AddrTimestamps::deref_mem(&self.ops_addr_usize[i], mem_val))
            .collect::<Vec<MultilinearPolynomial<F>>>()
    }
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
    PC: PolyCommitmentScheme<E>,
    E::ScalarField: Absorb,
{
    gens_ops: PCSKeys<E, PC>,
    gens_mem: PCSKeys<E, PC>,
    gens_derefs: PCSKeys<E, PC>,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> SparseMatPolyCommitmentKey<E, PC>
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

        let gens_ops = PC::trim(SRS);
        let gens_mem = PC::trim(SRS);
        let gens_derefs = PC::trim(SRS);
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
pub struct SparseMatPolyCommitment<E: Pairing, PC: PolyCommitmentScheme<E>>
where
    <E as Pairing>::ScalarField: Absorb,
{
    batch_size: usize,
    num_ops: usize,
    num_mem_cells: usize,
    comm_comb_ops: PC::Commitment,
    comm_comb_mem: PC::Commitment,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> AppendToTranscript<E::ScalarField>
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

    pub fn get_num_nz_entries(&self) -> usize {
        self.M.len().next_power_of_two()
    }

    fn sparse_to_dense_vecs(&self, N: usize) -> (Vec<usize>, Vec<usize>, Vec<F>) {
        assert!(N >= self.get_num_nz_entries());
        let mut ops_row: Vec<usize> = vec![0; N];
        let mut ops_col: Vec<usize> = vec![0; N];
        let mut val: Vec<F> = vec![F::zero(); N];

        for i in 0..self.M.len() {
            ops_row[i] = self.M[i].row;
            ops_col[i] = self.M[i].col;
            val[i] = self.M[i].val;
        }
        (ops_row, ops_col, val)
    }

    fn multi_sparse_to_dense_rep(
        sparse_polys: &[&SparseMatPolynomial<F>],
    ) -> MultiSparseMatPolynomialAsDense<F> {
        assert!(!sparse_polys.is_empty());
        for i in 1..sparse_polys.len() {
            assert_eq!(sparse_polys[i].num_vars_x, sparse_polys[0].num_vars_x);
            assert_eq!(sparse_polys[i].num_vars_y, sparse_polys[0].num_vars_y);
        }

        let N = (0..sparse_polys.len())
            .map(|i| sparse_polys[i].get_num_nz_entries())
            .max()
            .unwrap()
            .next_power_of_two();

        let mut ops_row_vec: Vec<Vec<usize>> = Vec::new();
        let mut ops_col_vec: Vec<Vec<usize>> = Vec::new();
        let mut val_vec: Vec<MultilinearPolynomial<F>> = Vec::new();
        for poly in sparse_polys {
            let (ops_row, ops_col, val) = poly.sparse_to_dense_vecs(N);
            ops_row_vec.push(ops_row);
            ops_col_vec.push(ops_col);
            val_vec.push(MultilinearPolynomial::new(val));
        }

        let any_poly = &sparse_polys[0];

        let num_mem_cells = if any_poly.num_vars_x > any_poly.num_vars_y {
            any_poly.num_vars_x.pow2()
        } else {
            any_poly.num_vars_y.pow2()
        };

        let row = AddrTimestamps::new(num_mem_cells, N, ops_row_vec);
        let col = AddrTimestamps::new(num_mem_cells, N, ops_col_vec);

        // combine polynomials into a single polynomial for commitment purposes
        let comb_ops = MultilinearPolynomial::merge(
            [
                row.ops_addr.as_slice(),
                row.read_ts.as_slice(),
                col.ops_addr.as_slice(),
                col.read_ts.as_slice(),
                val_vec.as_slice(),
            ]
                .concat()
                .as_slice(),
        );
        let mut comb_mem = row.audit_ts.clone();
        comb_mem.extend(&col.audit_ts);

        MultiSparseMatPolynomialAsDense {
            batch_size: sparse_polys.len(),
            row,
            col,
            val: val_vec,
            comb_ops,
            comb_mem,
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

    pub fn multi_commit<E: Pairing<ScalarField=F>, PC: PolyCommitmentScheme<E>>(
        sparse_polys: &[&SparseMatPolynomial<F>],
        gens: &SparseMatPolyCommitmentKey<E, PC>,
    ) -> (
        SparseMatPolyCommitment<E, PC>,
        MultiSparseMatPolynomialAsDense<F>,
    ) {
        let batch_size = sparse_polys.len();
        let dense = SparseMatPolynomial::multi_sparse_to_dense_rep(sparse_polys);

        let comm_comb_ops = PC::commit(&dense.comb_ops, &gens.gens_ops.ck);
        let comm_comb_mem = PC::commit(&dense.comb_mem, &gens.gens_mem.ck);

        (
            SparseMatPolyCommitment {
                batch_size,
                num_mem_cells: dense.row.audit_ts.len(),
                num_ops: dense.row.read_ts[0].len(),
                comm_comb_ops,
                comm_comb_mem,
            },
            dense,
        )
    }
}

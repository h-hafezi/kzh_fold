#![allow(clippy::type_complexity)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]
use super::errors::ProofVerifyError;
use super::math::Math;
use super::polycommitments::{PCSKeys, PolyCommitmentScheme};
use super::product_tree::{DotProductCircuit, ProductCircuit, ProductCircuitEvalProofBatched};
use super::timer::Timer;
use super::transcript::{AppendToTranscript, ProofTranscript};
use ark_ff::{Field, PrimeField};
use ark_serialize::*;
use ark_std::{cmp::max, One, Zero};
use ark_ec::pairing::Pairing;
use merlin::Transcript;
use crate::polynomial::eq_poly::EqPolynomial;
use crate::polynomial::identity::IdentityPolynomial;
use crate::polynomial::multilinear_poly::MultilinearPolynomial;

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

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct DerefsCommitment<E: Pairing, PC: PolyCommitmentScheme<E>> {
    comm_ops_val: PC::Commitment,
}

impl<F: PrimeField> Derefs<F> {
    pub fn new(row_ops_val: Vec<MultilinearPolynomial<F>>, col_ops_val: Vec<MultilinearPolynomial<F>>) -> Self {
        assert_eq!(row_ops_val.len(), col_ops_val.len());

        let derefs = {
            // combine all polynomials into a single polynomial (used below to produce a single commitment)
            let comb = MultilinearPolynomial::merge(
                [row_ops_val.as_slice(), col_ops_val.as_slice()]
                    .concat()
                    .as_slice(),
            );

            Derefs {
                row_ops_val,
                col_ops_val,
                comb,
            }
        };

        derefs
    }

    pub fn commit<E: Pairing<ScalarField = F>, PC: PolyCommitmentScheme<E>>(
        &self,
        gens: &PC::PolyCommitmentKey,
    ) -> DerefsCommitment<E, PC> {
        let comm_ops_val = PC::commit(&self.comb, gens);
        DerefsCommitment { comm_ops_val }
    }
}

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct DerefsEvalProof<E: Pairing, PC: PolyCommitmentScheme<E>> {
    proof_derefs: PC::PolyCommitmentProof,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> DerefsEvalProof<E, PC> {
    fn protocol_name() -> &'static [u8] {
        b"Derefs evaluation proof"
    }

    fn prove_single(
        joint_poly: &MultilinearPolynomial<E::ScalarField>,
        r: &[E::ScalarField],
        evals: Vec<E::ScalarField>,
        ck: &PC::PolyCommitmentKey,
        transcript: &mut Transcript,
    ) -> PC::PolyCommitmentProof {
        assert_eq!(joint_poly.get_num_vars(), r.len() + evals.len().log_2());

        // append the claimed evaluations to transcript
        <Transcript as ProofTranscript<E>>::append_scalars(transcript, b"evals_ops_val", &evals);

        // n-to-1 reduction
        let (r_joint, eval_joint) = {
            let challenges = <Transcript as ProofTranscript<E>>::challenge_vector(
                transcript,
                b"challenge_combine_n_to_one",
                evals.len().log_2(),
            );

            let mut poly_evals = MultilinearPolynomial::new(evals);
            for i in (0..challenges.len()).rev() {
                poly_evals.bound_poly_var_bot(&challenges[i]);
            }
            assert_eq!(poly_evals.len(), 1);
            let joint_claim_eval = poly_evals[0];
            let mut r_joint = challenges;
            r_joint.extend(r);

            debug_assert_eq!(joint_poly.evaluate(&r_joint), joint_claim_eval);
            (r_joint, joint_claim_eval)
        };
        // decommit the joint polynomial at r_joint
        <Transcript as ProofTranscript<E>>::append_scalar(transcript, b"joint_claim_eval", &eval_joint);

        PC::prove(None, joint_poly, &r_joint, &eval_joint, ck, transcript)
    }

    // evalues both polynomials at r and produces a joint proof of opening
    pub fn prove(
        derefs: &Derefs<E::ScalarField>,
        eval_row_ops_val_vec: &[E::ScalarField],
        eval_col_ops_val_vec: &[E::ScalarField],
        r: &[E::ScalarField],
        ck: &PC::PolyCommitmentKey,
        transcript: &mut Transcript,
    ) -> Self {
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            DerefsEvalProof::<E, PC>::protocol_name(),
        );

        let evals = {
            let mut evals = eval_row_ops_val_vec.to_owned();
            evals.extend(eval_col_ops_val_vec);
            evals.resize(evals.len().next_power_of_two(), E::ScalarField::zero());
            evals
        };
        let proof_derefs =
            DerefsEvalProof::<E, PC>::prove_single(&derefs.comb, r, evals, ck, transcript);

        DerefsEvalProof { proof_derefs }
    }

    fn verify_single(
        proof: &PC::PolyCommitmentProof,
        comm: &PC::Commitment,
        r: &[E::ScalarField],
        evals: Vec<E::ScalarField>,
        vk: &PC::EvalVerifierKey,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        // append the claimed evaluations to transcript
        <Transcript as ProofTranscript<E>>::append_scalars(transcript, b"evals_ops_val", &evals);

        // n-to-1 reduction
        let challenges = <Transcript as ProofTranscript<E>>::challenge_vector(
            transcript,
            b"challenge_combine_n_to_one",
            evals.len().log_2(),
        );
        let mut poly_evals = MultilinearPolynomial::new(evals);
        for i in (0..challenges.len()).rev() {
            poly_evals.bound_poly_var_bot(&challenges[i]);
        }
        assert_eq!(poly_evals.len(), 1);
        let joint_claim_eval = poly_evals[0];
        let mut r_joint = challenges;
        r_joint.extend(r);

        // decommit the joint polynomial at r_joint
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"joint_claim_eval",
            &joint_claim_eval,
        );

        PC::verify(comm, proof, vk, transcript, &r_joint, &joint_claim_eval).map_err(|e| e.into())
    }

    // verify evaluations of both polynomials at r
    pub fn verify(
        &self,
        r: &[E::ScalarField],
        eval_row_ops_val_vec: &[E::ScalarField],
        eval_col_ops_val_vec: &[E::ScalarField],
        vk: &PC::EvalVerifierKey,
        comm: &DerefsCommitment<E, PC>,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            DerefsEvalProof::<E, PC>::protocol_name(),
        );
        let mut evals = eval_row_ops_val_vec.to_owned();
        evals.extend(eval_col_ops_val_vec);
        evals.resize(evals.len().next_power_of_two(), E::ScalarField::zero());

        DerefsEvalProof::<E, PC>::verify_single(
            &self.proof_derefs,
            &comm.comm_ops_val,
            r,
            evals,
            vk,
            transcript,
        )
    }
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> AppendToTranscript<E> for DerefsCommitment<E, PC> {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_message(b"derefs_commitment", b"begin_derefs_commitment");
        self.comm_ops_val.append_to_transcript(label, transcript);
        transcript.append_message(b"derefs_commitment", b"end_derefs_commitment");
    }
}

#[derive(CanonicalDeserialize, CanonicalSerialize)]
struct AddrTimestamps<F>
where
    F: Sync + CanonicalSerialize + CanonicalDeserialize + PrimeField,
{
    ops_addr_usize: Vec<Vec<usize>>,
    ops_addr: Vec<MultilinearPolynomial<F>>,
    read_ts: Vec<MultilinearPolynomial<F>>,
    audit_ts: MultilinearPolynomial<F>,
}

impl<F: PrimeField> AddrTimestamps<F> {
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
    F: Sync + CanonicalSerialize + CanonicalDeserialize + PrimeField,
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
{
    gens_ops: PCSKeys<E, PC>,
    gens_mem: PCSKeys<E, PC>,
    gens_derefs: PCSKeys<E, PC>,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> SparseMatPolyCommitmentKey<E, PC> {
    pub fn new(
        SRS: &PC::SRS,
        num_vars_x: usize,
        num_vars_y: usize,
        num_nz_entries: usize,
        batch_size: usize,
    ) -> SparseMatPolyCommitmentKey<E, PC> {
        let (num_vars_ops, num_vars_mem, num_vars_derefs) =
            Self::get_gens_sizes(num_vars_x, num_vars_y, num_nz_entries, batch_size);

        let gens_ops = PC::trim(SRS, num_vars_ops);
        let gens_mem = PC::trim(SRS, num_vars_mem);
        let gens_derefs = PC::trim(SRS, num_vars_derefs);
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
pub struct SparseMatPolyCommitment<E: Pairing, PC: PolyCommitmentScheme<E>> {
    batch_size: usize,
    num_ops: usize,
    num_mem_cells: usize,
    comm_comb_ops: PC::Commitment,
    comm_comb_mem: PC::Commitment,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> AppendToTranscript<E>
for SparseMatPolyCommitment<E, PC>
{
    fn append_to_transcript(&self, _label: &'static [u8], transcript: &mut Transcript) {
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

impl<F: PrimeField> SparseMatPolynomial<F> {
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

    pub fn multi_commit<E: Pairing<ScalarField = F>, PC: PolyCommitmentScheme<E>>(
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

impl<F: PrimeField> MultiSparseMatPolynomialAsDense<F> {
    pub fn deref(&self, row_mem_val: &[F], col_mem_val: &[F]) -> Derefs<F> {
        let row_ops_val = self.row.deref(row_mem_val);
        let col_ops_val = self.col.deref(col_mem_val);

        Derefs::new(row_ops_val, col_ops_val)
    }
}

#[derive(Debug)]
struct ProductLayer<F>
where
    F: Sync + CanonicalSerialize + CanonicalDeserialize + PrimeField,
{
    init: ProductCircuit<F>,
    read_vec: Vec<ProductCircuit<F>>,
    write_vec: Vec<ProductCircuit<F>>,
    audit: ProductCircuit<F>,
}

#[derive(Debug)]
struct Layers<F>
where
    F: Sync + CanonicalSerialize + CanonicalDeserialize + PrimeField,
{
    prod_layer: ProductLayer<F>,
}

impl<F: PrimeField> Layers<F> {
    fn build_hash_layer(
        eval_table: &[F],
        addrs_vec: &[MultilinearPolynomial<F>],
        derefs_vec: &[MultilinearPolynomial<F>],
        read_ts_vec: &[MultilinearPolynomial<F>],
        audit_ts: &MultilinearPolynomial<F>,
        r_mem_check: &(F, F),
    ) -> (
        MultilinearPolynomial<F>,
        Vec<MultilinearPolynomial<F>>,
        Vec<MultilinearPolynomial<F>>,
        MultilinearPolynomial<F>,
    ) {
        let (r_hash, r_multiset_check) = r_mem_check;

        //hash(addr, val, ts) = ts * r_hash_sqr + val * r_hash + addr
        let r_hash_sqr = r_hash.square();
        let hash_func = |addr: &F, val: &F, ts: &F| -> F { *ts * r_hash_sqr + *val * *r_hash + *addr };

        // hash init and audit that does not depend on #instances
        let num_mem_cells = eval_table.len();
        let poly_init_hashed = MultilinearPolynomial::new(
            (0..num_mem_cells)
                .map(|i| {
                    // at init time, addr is given by i, init value is given by eval_table, and ts = 0
                    hash_func(&F::from(i as u64), &eval_table[i], &F::zero()) - r_multiset_check
                })
                .collect::<Vec<F>>(),
        );
        let poly_audit_hashed = MultilinearPolynomial::new(
            (0..num_mem_cells)
                .map(|i| {
                    // at audit time, addr is given by i, value is given by eval_table, and ts is given by audit_ts
                    hash_func(&F::from(i as u64), &eval_table[i], &audit_ts[i]) - r_multiset_check
                })
                .collect::<Vec<F>>(),
        );

        // hash read and write that depends on #instances
        let mut poly_read_hashed_vec: Vec<MultilinearPolynomial<F>> = Vec::new();
        let mut poly_write_hashed_vec: Vec<MultilinearPolynomial<F>> = Vec::new();
        for i in 0..addrs_vec.len() {
            let (addrs, derefs, read_ts) = (&addrs_vec[i], &derefs_vec[i], &read_ts_vec[i]);
            assert_eq!(addrs.len(), derefs.len());
            assert_eq!(addrs.len(), read_ts.len());
            let num_ops = addrs.len();
            let poly_read_hashed = MultilinearPolynomial::new(
                (0..num_ops)
                    .map(|i| {
                        // at read time, addr is given by addrs, value is given by derefs, and ts is given by read_ts
                        hash_func(&addrs[i], &derefs[i], &read_ts[i]) - r_multiset_check
                    })
                    .collect::<Vec<F>>(),
            );
            poly_read_hashed_vec.push(poly_read_hashed);

            let poly_write_hashed = MultilinearPolynomial::new(
                (0..num_ops)
                    .map(|i| {
                        // at write time, addr is given by addrs, value is given by derefs, and ts is given by write_ts = read_ts + 1
                        hash_func(&addrs[i], &derefs[i], &(read_ts[i] + F::one())) - r_multiset_check
                    })
                    .collect::<Vec<F>>(),
            );
            poly_write_hashed_vec.push(poly_write_hashed);
        }

        (
            poly_init_hashed,
            poly_read_hashed_vec,
            poly_write_hashed_vec,
            poly_audit_hashed,
        )
    }

    pub fn new(
        eval_table: &[F],
        addr_timestamps: &AddrTimestamps<F>,
        poly_ops_val: &[MultilinearPolynomial<F>],
        r_mem_check: &(F, F),
    ) -> Self {
        let (poly_init_hashed, poly_read_hashed_vec, poly_write_hashed_vec, poly_audit_hashed) =
            Layers::build_hash_layer(
                eval_table,
                &addr_timestamps.ops_addr,
                poly_ops_val,
                &addr_timestamps.read_ts,
                &addr_timestamps.audit_ts,
                r_mem_check,
            );

        let prod_init = ProductCircuit::new(&poly_init_hashed);
        let prod_read_vec = (0..poly_read_hashed_vec.len())
            .map(|i| ProductCircuit::new(&poly_read_hashed_vec[i]))
            .collect::<Vec<ProductCircuit<F>>>();
        let prod_write_vec = (0..poly_write_hashed_vec.len())
            .map(|i| ProductCircuit::new(&poly_write_hashed_vec[i]))
            .collect::<Vec<ProductCircuit<F>>>();
        let prod_audit = ProductCircuit::new(&poly_audit_hashed);

        // subset audit check
        let hashed_writes: F = (0..prod_write_vec.len())
            .map(|i| prod_write_vec[i].evaluate())
            .product();
        let hashed_write_set: F = prod_init.evaluate() * hashed_writes;

        let hashed_reads: F = (0..prod_read_vec.len())
            .map(|i| prod_read_vec[i].evaluate())
            .product();
        let hashed_read_set: F = hashed_reads * prod_audit.evaluate();

        //assert_eq!(hashed_read_set, hashed_write_set);
        debug_assert_eq!(hashed_read_set, hashed_write_set);

        Layers {
            prod_layer: ProductLayer {
                init: prod_init,
                read_vec: prod_read_vec,
                write_vec: prod_write_vec,
                audit: prod_audit,
            },
        }
    }
}

#[derive(Debug)]
struct PolyEvalNetwork<F>
where
    F: Sync + CanonicalDeserialize + CanonicalSerialize + PrimeField,
{
    row_layers: Layers<F>,
    col_layers: Layers<F>,
}

impl<F: PrimeField> PolyEvalNetwork<F> {
    pub fn new(
        dense: &MultiSparseMatPolynomialAsDense<F>,
        derefs: &Derefs<F>,
        mem_rx: &[F],
        mem_ry: &[F],
        r_mem_check: &(F, F),
    ) -> Self {
        let row_layers = Layers::new(mem_rx, &dense.row, &derefs.row_ops_val, r_mem_check);
        let col_layers = Layers::new(mem_ry, &dense.col, &derefs.col_ops_val, r_mem_check);

        PolyEvalNetwork {
            row_layers,
            col_layers,
        }
    }
}

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
struct HashLayerProof<E: Pairing, PC: PolyCommitmentScheme<E>> {
    eval_row: (Vec<E::ScalarField>, Vec<E::ScalarField>, E::ScalarField),
    eval_col: (Vec<E::ScalarField>, Vec<E::ScalarField>, E::ScalarField),
    eval_val: Vec<E::ScalarField>,
    eval_derefs: (Vec<E::ScalarField>, Vec<E::ScalarField>),
    proof_ops: PC::PolyCommitmentProof,
    proof_mem: PC::PolyCommitmentProof,
    proof_derefs: DerefsEvalProof<E, PC>,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> HashLayerProof<E, PC> {
    fn protocol_name() -> &'static [u8] {
        b"Sparse polynomial hash layer proof"
    }

    fn prove_helper(
        rand: (&Vec<E::ScalarField>, &Vec<E::ScalarField>),
        addr_timestamps: &AddrTimestamps<E::ScalarField>,
    ) -> (Vec<E::ScalarField>, Vec<E::ScalarField>, E::ScalarField) {
        let (rand_mem, rand_ops) = rand;

        // decommit ops-addr at rand_ops
        let mut eval_ops_addr_vec: Vec<E::ScalarField> = Vec::new();
        for i in 0..addr_timestamps.ops_addr.len() {
            let eval_ops_addr = addr_timestamps.ops_addr[i].evaluate(rand_ops);
            eval_ops_addr_vec.push(eval_ops_addr);
        }

        // decommit read_ts at rand_ops
        let mut eval_read_ts_vec: Vec<E::ScalarField> = Vec::new();
        for i in 0..addr_timestamps.read_ts.len() {
            let eval_read_ts = addr_timestamps.read_ts[i].evaluate(rand_ops);
            eval_read_ts_vec.push(eval_read_ts);
        }

        // decommit audit-ts at rand_mem
        let eval_audit_ts = addr_timestamps.audit_ts.evaluate(rand_mem);

        (eval_ops_addr_vec, eval_read_ts_vec, eval_audit_ts)
    }

    fn prove(
        rand: (&Vec<E::ScalarField>, &Vec<E::ScalarField>),
        dense: &MultiSparseMatPolynomialAsDense<E::ScalarField>,
        derefs: &Derefs<E::ScalarField>,
        gens: &SparseMatPolyCommitmentKey<E, PC>,
        transcript: &mut Transcript,
    ) -> Self {
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            HashLayerProof::<E, PC>::protocol_name(),
        );

        let (rand_mem, rand_ops) = rand;

        // decommit derefs at rand_ops
        let eval_row_ops_val = (0..derefs.row_ops_val.len())
            .map(|i| derefs.row_ops_val[i].evaluate(rand_ops))
            .collect::<Vec<E::ScalarField>>();
        let eval_col_ops_val = (0..derefs.col_ops_val.len())
            .map(|i| derefs.col_ops_val[i].evaluate(rand_ops))
            .collect::<Vec<E::ScalarField>>();
        let proof_derefs = DerefsEvalProof::prove(
            derefs,
            &eval_row_ops_val,
            &eval_col_ops_val,
            rand_ops,
            &gens.gens_derefs.ck,
            transcript,
        );
        let eval_derefs = (eval_row_ops_val, eval_col_ops_val);

        // evaluate row_addr, row_read-ts, col_addr, col_read-ts, val at rand_ops
        // evaluate row_audit_ts and col_audit_ts at rand_mem
        let (eval_row_addr_vec, eval_row_read_ts_vec, eval_row_audit_ts) =
            HashLayerProof::<E, PC>::prove_helper((rand_mem, rand_ops), &dense.row);
        let (eval_col_addr_vec, eval_col_read_ts_vec, eval_col_audit_ts) =
            HashLayerProof::<E, PC>::prove_helper((rand_mem, rand_ops), &dense.col);
        let eval_val_vec = (0..dense.val.len())
            .map(|i| dense.val[i].evaluate(rand_ops))
            .collect::<Vec<E::ScalarField>>();

        // form a single decommitment using comm_comb_ops
        let mut evals_ops: Vec<E::ScalarField> = Vec::new();
        evals_ops.extend(&eval_row_addr_vec);
        evals_ops.extend(&eval_row_read_ts_vec);
        evals_ops.extend(&eval_col_addr_vec);
        evals_ops.extend(&eval_col_read_ts_vec);
        evals_ops.extend(&eval_val_vec);
        evals_ops.resize(evals_ops.len().next_power_of_two(), E::ScalarField::zero());

        <Transcript as ProofTranscript<E>>::append_scalars(transcript, b"claim_evals_ops", &evals_ops);

        let challenges_ops = <Transcript as ProofTranscript<E>>::challenge_vector(
            transcript,
            b"challenge_combine_n_to_one",
            evals_ops.len().log_2(),
        );

        let mut poly_evals_ops = MultilinearPolynomial::new(evals_ops);
        for i in (0..challenges_ops.len()).rev() {
            poly_evals_ops.bound_poly_var_bot(&challenges_ops[i]);
        }
        assert_eq!(poly_evals_ops.len(), 1);
        let joint_claim_eval_ops = poly_evals_ops[0];
        let mut r_joint_ops = challenges_ops;
        r_joint_ops.extend(rand_ops);
        debug_assert_eq!(
            dense.comb_ops.evaluate(&r_joint_ops),
            joint_claim_eval_ops
        );
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"joint_claim_eval_ops",
            &joint_claim_eval_ops,
        );

        let proof_ops = PC::prove(
            None,
            &dense.comb_ops,
            &r_joint_ops,
            &joint_claim_eval_ops,
            &gens.gens_ops.ck,
            transcript,
        );

        // form a single decommitment using comb_comb_mem at rand_mem
        let evals_mem: Vec<E::ScalarField> = vec![eval_row_audit_ts, eval_col_audit_ts];

        <Transcript as ProofTranscript<E>>::append_scalars(transcript, b"claim_evals_mem", &evals_mem);
        let challenges_mem = <Transcript as ProofTranscript<E>>::challenge_vector(
            transcript,
            b"challenge_combine_two_to_one",
            evals_mem.len().log_2(),
        );

        let mut poly_evals_mem = MultilinearPolynomial::new(evals_mem);
        for i in (0..challenges_mem.len()).rev() {
            poly_evals_mem.bound_poly_var_bot(&challenges_mem[i]);
        }
        assert_eq!(poly_evals_mem.len(), 1);
        let joint_claim_eval_mem = poly_evals_mem[0];
        let mut r_joint_mem = challenges_mem;
        r_joint_mem.extend(rand_mem);
        debug_assert_eq!(
            dense.comb_mem.evaluate(&r_joint_mem),
            joint_claim_eval_mem
        );
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"joint_claim_eval_mem",
            &joint_claim_eval_mem,
        );

        let proof_mem = PC::prove(
            None,
            &dense.comb_mem,
            &r_joint_mem,
            &joint_claim_eval_mem,
            &gens.gens_mem.ck,
            transcript,
        );

        HashLayerProof {
            eval_row: (eval_row_addr_vec, eval_row_read_ts_vec, eval_row_audit_ts),
            eval_col: (eval_col_addr_vec, eval_col_read_ts_vec, eval_col_audit_ts),
            eval_val: eval_val_vec,
            eval_derefs,
            proof_ops,
            proof_mem,
            proof_derefs,
        }
    }

    fn verify_helper(
        rand: &(&Vec<E::ScalarField>, &Vec<E::ScalarField>),
        claims: &(
            E::ScalarField,
            Vec<E::ScalarField>,
            Vec<E::ScalarField>,
            E::ScalarField,
        ),
        eval_ops_val: &[E::ScalarField],
        eval_ops_addr: &[E::ScalarField],
        eval_read_ts: &[E::ScalarField],
        eval_audit_ts: &E::ScalarField,
        r: &[E::ScalarField],
        r_hash: &E::ScalarField,
        r_multiset_check: &E::ScalarField,
    ) -> Result<(), ProofVerifyError> {
        let r_hash_sqr = r_hash.square();
        let hash_func = |addr: &E::ScalarField,
                         val: &E::ScalarField,
                         ts: &E::ScalarField|
                         -> E::ScalarField { *ts * r_hash_sqr + *val * *r_hash + *addr };

        let (rand_mem, _rand_ops) = rand;
        let (claim_init, claim_read, claim_write, claim_audit) = claims;

        // init
        let eval_init_addr = IdentityPolynomial::new(rand_mem.len()).evaluate(rand_mem);
        let eval_init_val = EqPolynomial::new(r.to_vec()).evaluate(rand_mem);
        let hash_init_at_rand_mem =
            hash_func(&eval_init_addr, &eval_init_val, &E::ScalarField::zero()) - r_multiset_check; // verify the claim_last of init chunk
        assert_eq!(&hash_init_at_rand_mem, claim_init);

        // read
        for i in 0..eval_ops_addr.len() {
            let hash_read_at_rand_ops =
                hash_func(&eval_ops_addr[i], &eval_ops_val[i], &eval_read_ts[i]) - r_multiset_check; // verify the claim_last of init chunk
            assert_eq!(&hash_read_at_rand_ops, &claim_read[i]);
        }

        // write: shares addr, val component; only decommit write_ts
        for i in 0..eval_ops_addr.len() {
            let eval_write_ts = eval_read_ts[i] + E::ScalarField::one();
            let hash_write_at_rand_ops =
                hash_func(&eval_ops_addr[i], &eval_ops_val[i], &eval_write_ts) - r_multiset_check; // verify the claim_last of init chunk
            assert_eq!(&hash_write_at_rand_ops, &claim_write[i]);
        }

        // audit: shares addr and val with init
        let eval_audit_addr = eval_init_addr;
        let eval_audit_val = eval_init_val;
        let hash_audit_at_rand_mem =
            hash_func(&eval_audit_addr, &eval_audit_val, eval_audit_ts) - r_multiset_check;
        assert_eq!(&hash_audit_at_rand_mem, claim_audit); // verify the last step of the sum-check for audit

        Ok(())
    }

    fn verify(
        &self,
        rand: (&Vec<E::ScalarField>, &Vec<E::ScalarField>),
        claims_row: &(
            E::ScalarField,
            Vec<E::ScalarField>,
            Vec<E::ScalarField>,
            E::ScalarField,
        ),
        claims_col: &(
            E::ScalarField,
            Vec<E::ScalarField>,
            Vec<E::ScalarField>,
            E::ScalarField,
        ),
        claims_dotp: &[E::ScalarField],
        comm: &SparseMatPolyCommitment<E, PC>,
        gens: &SparseMatPolyCommitmentKey<E, PC>,
        comm_derefs: &DerefsCommitment<E, PC>,
        rx: &[E::ScalarField],
        ry: &[E::ScalarField],
        r_hash: &E::ScalarField,
        r_multiset_check: &E::ScalarField,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        let timer = Timer::new("verify_hash_proof");
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            HashLayerProof::<E, PC>::protocol_name(),
        );

        let (rand_mem, rand_ops) = rand;

        // verify derefs at rand_ops
        let (eval_row_ops_val, eval_col_ops_val) = &self.eval_derefs;
        assert_eq!(eval_row_ops_val.len(), eval_col_ops_val.len());
        self.proof_derefs.verify(
            rand_ops,
            eval_row_ops_val,
            eval_col_ops_val,
            &gens.gens_derefs.vk,
            comm_derefs,
            transcript,
        )?;

        // verify the decommitments used in evaluation sum-check
        let eval_val_vec = &self.eval_val;
        assert_eq!(claims_dotp.len(), 3 * eval_row_ops_val.len());
        for i in 0..claims_dotp.len() / 3 {
            let claim_row_ops_val = claims_dotp[3 * i];
            let claim_col_ops_val = claims_dotp[3 * i + 1];
            let claim_val = claims_dotp[3 * i + 2];

            assert_eq!(claim_row_ops_val, eval_row_ops_val[i]);
            assert_eq!(claim_col_ops_val, eval_col_ops_val[i]);
            assert_eq!(claim_val, eval_val_vec[i]);
        }

        // verify addr-timestamps using comm_comb_ops at rand_ops
        let (eval_row_addr_vec, eval_row_read_ts_vec, eval_row_audit_ts) = &self.eval_row;
        let (eval_col_addr_vec, eval_col_read_ts_vec, eval_col_audit_ts) = &self.eval_col;

        let mut evals_ops: Vec<E::ScalarField> = Vec::new();
        evals_ops.extend(eval_row_addr_vec);
        evals_ops.extend(eval_row_read_ts_vec);
        evals_ops.extend(eval_col_addr_vec);
        evals_ops.extend(eval_col_read_ts_vec);
        evals_ops.extend(eval_val_vec);
        evals_ops.resize(evals_ops.len().next_power_of_two(), E::ScalarField::zero());

        <Transcript as ProofTranscript<E>>::append_scalars(transcript, b"claim_evals_ops", &evals_ops);

        let challenges_ops = <Transcript as ProofTranscript<E>>::challenge_vector(
            transcript,
            b"challenge_combine_n_to_one",
            evals_ops.len().log_2(),
        );

        let mut poly_evals_ops = MultilinearPolynomial::new(evals_ops);
        for i in (0..challenges_ops.len()).rev() {
            poly_evals_ops.bound_poly_var_bot(&challenges_ops[i]);
        }
        assert_eq!(poly_evals_ops.len(), 1);
        let joint_claim_eval_ops = poly_evals_ops[0];
        let mut r_joint_ops = challenges_ops;
        r_joint_ops.extend(rand_ops);
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"joint_claim_eval_ops",
            &joint_claim_eval_ops,
        );
        PC::verify(
            &comm.comm_comb_ops,
            &self.proof_ops,
            &gens.gens_ops.vk,
            transcript,
            &r_joint_ops,
            &joint_claim_eval_ops,
        )?;

        // verify proof-mem using comm_comb_mem at rand_mem
        // form a single decommitment using comb_comb_mem at rand_mem
        let evals_mem: Vec<E::ScalarField> = vec![*eval_row_audit_ts, *eval_col_audit_ts];
        <Transcript as ProofTranscript<E>>::append_scalars(transcript, b"claim_evals_mem", &evals_mem);
        let challenges_mem = <Transcript as ProofTranscript<E>>::challenge_vector(
            transcript,
            b"challenge_combine_two_to_one",
            evals_mem.len().log_2(),
        );

        let mut poly_evals_mem = MultilinearPolynomial::new(evals_mem);
        for i in (0..challenges_mem.len()).rev() {
            poly_evals_mem.bound_poly_var_bot(&challenges_mem[i]);
        }
        assert_eq!(poly_evals_mem.len(), 1);
        let joint_claim_eval_mem = poly_evals_mem[0];
        let mut r_joint_mem = challenges_mem;
        r_joint_mem.extend(rand_mem);
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"joint_claim_eval_mem",
            &joint_claim_eval_mem,
        );
        PC::verify(
            &comm.comm_comb_mem,
            &self.proof_mem,
            &gens.gens_mem.vk,
            transcript,
            &r_joint_mem,
            &joint_claim_eval_mem,
        )?;

        // verify the claims from the product layer
        let (eval_ops_addr, eval_read_ts, eval_audit_ts) = &self.eval_row;
        HashLayerProof::<E, PC>::verify_helper(
            &(rand_mem, rand_ops),
            claims_row,
            eval_row_ops_val,
            eval_ops_addr,
            eval_read_ts,
            eval_audit_ts,
            rx,
            r_hash,
            r_multiset_check,
        )?;

        let (eval_ops_addr, eval_read_ts, eval_audit_ts) = &self.eval_col;
        HashLayerProof::<E, PC>::verify_helper(
            &(rand_mem, rand_ops),
            claims_col,
            eval_col_ops_val,
            eval_ops_addr,
            eval_read_ts,
            eval_audit_ts,
            ry,
            r_hash,
            r_multiset_check,
        )?;

        timer.stop();
        Ok(())
    }
}

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
struct ProductLayerProof<F: PrimeField> {
    eval_row: (F, Vec<F>, Vec<F>, F),
    eval_col: (F, Vec<F>, Vec<F>, F),
    eval_val: (Vec<F>, Vec<F>),
    proof_mem: ProductCircuitEvalProofBatched<F>,
    proof_ops: ProductCircuitEvalProofBatched<F>,
}

impl<F: PrimeField> ProductLayerProof<F> {
    fn protocol_name() -> &'static [u8] {
        b"Sparse polynomial product layer proof"
    }

    pub fn prove<E>(
        row_prod_layer: &mut ProductLayer<F>,
        col_prod_layer: &mut ProductLayer<F>,
        dense: &MultiSparseMatPolynomialAsDense<F>,
        derefs: &Derefs<F>,
        eval: &[F],
        transcript: &mut Transcript,
    ) -> (Self, Vec<F>, Vec<F>)
    where
        E: Pairing<ScalarField = F>,
    {
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            ProductLayerProof::<F>::protocol_name(),
        );

        let row_eval_init = row_prod_layer.init.evaluate();
        let row_eval_audit = row_prod_layer.audit.evaluate();
        let row_eval_read = (0..row_prod_layer.read_vec.len())
            .map(|i| row_prod_layer.read_vec[i].evaluate())
            .collect::<Vec<F>>();
        let row_eval_write = (0..row_prod_layer.write_vec.len())
            .map(|i| row_prod_layer.write_vec[i].evaluate())
            .collect::<Vec<F>>();

        // subset check
        let ws: F = (0..row_eval_write.len())
            .map(|i| row_eval_write[i])
            .product();
        let rs: F = (0..row_eval_read.len()).map(|i| row_eval_read[i]).product();
        assert_eq!(row_eval_init * ws, rs * row_eval_audit);

        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_row_eval_init",
            &row_eval_init,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_row_eval_read",
            &row_eval_read,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_row_eval_write",
            &row_eval_write,
        );
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_row_eval_audit",
            &row_eval_audit,
        );

        let col_eval_init = col_prod_layer.init.evaluate();
        let col_eval_audit = col_prod_layer.audit.evaluate();
        let col_eval_read: Vec<F> = (0..col_prod_layer.read_vec.len())
            .map(|i| col_prod_layer.read_vec[i].evaluate())
            .collect();
        let col_eval_write: Vec<F> = (0..col_prod_layer.write_vec.len())
            .map(|i| col_prod_layer.write_vec[i].evaluate())
            .collect();

        // subset check
        let ws: F = (0..col_eval_write.len())
            .map(|i| col_eval_write[i])
            .product();
        let rs: F = (0..col_eval_read.len()).map(|i| col_eval_read[i]).product();
        assert_eq!(col_eval_init * ws, rs * col_eval_audit);

        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_col_eval_init",
            &col_eval_init,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_col_eval_read",
            &col_eval_read,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_col_eval_write",
            &col_eval_write,
        );
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_col_eval_audit",
            &col_eval_audit,
        );

        // prepare dotproduct circuit for batching then with ops-related product circuits
        assert_eq!(eval.len(), derefs.row_ops_val.len());
        assert_eq!(eval.len(), derefs.col_ops_val.len());
        assert_eq!(eval.len(), dense.val.len());
        let mut dotp_circuit_left_vec: Vec<DotProductCircuit<F>> = Vec::new();
        let mut dotp_circuit_right_vec: Vec<DotProductCircuit<F>> = Vec::new();
        let mut eval_dotp_left_vec: Vec<F> = Vec::new();
        let mut eval_dotp_right_vec: Vec<F> = Vec::new();
        for i in 0..derefs.row_ops_val.len() {
            // evaluate sparse polynomial evaluation using two dotp checks
            let left = derefs.row_ops_val[i].clone();
            let right = derefs.col_ops_val[i].clone();
            let weights = dense.val[i].clone();

            // build two dot product circuits to prove evaluation of sparse polynomial
            let mut dotp_circuit = DotProductCircuit::new(left, right, weights);
            let (dotp_circuit_left, dotp_circuit_right) = dotp_circuit.split();

            let (eval_dotp_left, eval_dotp_right) =
                (dotp_circuit_left.evaluate(), dotp_circuit_right.evaluate());

            <Transcript as ProofTranscript<E>>::append_scalar(
                transcript,
                b"claim_eval_dotp_left",
                &eval_dotp_left,
            );

            <Transcript as ProofTranscript<E>>::append_scalar(
                transcript,
                b"claim_eval_dotp_right",
                &eval_dotp_right,
            );

            assert_eq!(eval_dotp_left + eval_dotp_right, eval[i]);
            eval_dotp_left_vec.push(eval_dotp_left);
            eval_dotp_right_vec.push(eval_dotp_right);

            dotp_circuit_left_vec.push(dotp_circuit_left);
            dotp_circuit_right_vec.push(dotp_circuit_right);
        }

        // The number of operations into the memory encoded by rx and ry are always the same (by design)
        // So we can produce a batched product proof for all of them at the same time.
        // prove the correctness of claim_row_eval_read, claim_row_eval_write, claim_col_eval_read, and claim_col_eval_write
        // TODO: we currently only produce proofs for 3 batched sparse polynomial evaluations
        assert_eq!(row_prod_layer.read_vec.len(), 3);
        let (row_read_A, row_read_B, row_read_C) = {
            let (vec_A, vec_BC) = row_prod_layer.read_vec.split_at_mut(1);
            let (vec_B, vec_C) = vec_BC.split_at_mut(1);
            (vec_A, vec_B, vec_C)
        };

        let (row_write_A, row_write_B, row_write_C) = {
            let (vec_A, vec_BC) = row_prod_layer.write_vec.split_at_mut(1);
            let (vec_B, vec_C) = vec_BC.split_at_mut(1);
            (vec_A, vec_B, vec_C)
        };

        let (col_read_A, col_read_B, col_read_C) = {
            let (vec_A, vec_BC) = col_prod_layer.read_vec.split_at_mut(1);
            let (vec_B, vec_C) = vec_BC.split_at_mut(1);
            (vec_A, vec_B, vec_C)
        };

        let (col_write_A, col_write_B, col_write_C) = {
            let (vec_A, vec_BC) = col_prod_layer.write_vec.split_at_mut(1);
            let (vec_B, vec_C) = vec_BC.split_at_mut(1);
            (vec_A, vec_B, vec_C)
        };

        let (dotp_left_A, dotp_left_B, dotp_left_C) = {
            let (vec_A, vec_BC) = dotp_circuit_left_vec.split_at_mut(1);
            let (vec_B, vec_C) = vec_BC.split_at_mut(1);
            (vec_A, vec_B, vec_C)
        };

        let (dotp_right_A, dotp_right_B, dotp_right_C) = {
            let (vec_A, vec_BC) = dotp_circuit_right_vec.split_at_mut(1);
            let (vec_B, vec_C) = vec_BC.split_at_mut(1);
            (vec_A, vec_B, vec_C)
        };

        let (proof_ops, rand_ops) = ProductCircuitEvalProofBatched::<F>::prove::<E>(
            &mut [
                &mut row_read_A[0],
                &mut row_read_B[0],
                &mut row_read_C[0],
                &mut row_write_A[0],
                &mut row_write_B[0],
                &mut row_write_C[0],
                &mut col_read_A[0],
                &mut col_read_B[0],
                &mut col_read_C[0],
                &mut col_write_A[0],
                &mut col_write_B[0],
                &mut col_write_C[0],
            ],
            &mut [
                &mut dotp_left_A[0],
                &mut dotp_right_A[0],
                &mut dotp_left_B[0],
                &mut dotp_right_B[0],
                &mut dotp_left_C[0],
                &mut dotp_right_C[0],
            ],
            transcript,
        );

        // produce a batched proof of memory-related product circuits
        let (proof_mem, rand_mem) = ProductCircuitEvalProofBatched::<F>::prove::<E>(
            &mut [
                &mut row_prod_layer.init,
                &mut row_prod_layer.audit,
                &mut col_prod_layer.init,
                &mut col_prod_layer.audit,
            ],
            &mut Vec::new(),
            transcript,
        );

        let product_layer_proof = ProductLayerProof {
            eval_row: (row_eval_init, row_eval_read, row_eval_write, row_eval_audit),
            eval_col: (col_eval_init, col_eval_read, col_eval_write, col_eval_audit),
            eval_val: (eval_dotp_left_vec, eval_dotp_right_vec),
            proof_mem,
            proof_ops,
        };

        let mut product_layer_proof_encoded = vec![];
        product_layer_proof
            .serialize_compressed(&mut product_layer_proof_encoded)
            .unwrap();

        let msg = format!(
            "len_product_layer_proof {:?}",
            product_layer_proof_encoded.len()
        );
        Timer::print(&msg);

        (product_layer_proof, rand_mem, rand_ops)
    }

    pub fn verify<E>(
        &self,
        num_ops: usize,
        num_cells: usize,
        eval: &[F],
        transcript: &mut Transcript,
    ) -> Result<(Vec<F>, Vec<F>, Vec<F>, Vec<F>, Vec<F>), ProofVerifyError>
    where
        E: Pairing<ScalarField = F>,
    {
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            ProductLayerProof::<F>::protocol_name(),
        );

        let timer = Timer::new("verify_prod_proof");
        let num_instances = eval.len();

        // subset check
        let (row_eval_init, row_eval_read, row_eval_write, row_eval_audit) = &self.eval_row;
        assert_eq!(row_eval_write.len(), num_instances);
        assert_eq!(row_eval_read.len(), num_instances);
        let ws: F = (0..row_eval_write.len())
            .map(|i| row_eval_write[i])
            .product();
        let rs: F = (0..row_eval_read.len()).map(|i| row_eval_read[i]).product();
        assert_eq!(*row_eval_init * ws, rs * row_eval_audit);

        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_row_eval_init",
            row_eval_init,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_row_eval_read",
            row_eval_read,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_row_eval_write",
            row_eval_write,
        );
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_row_eval_audit",
            row_eval_audit,
        );

        // subset check
        let (col_eval_init, col_eval_read, col_eval_write, col_eval_audit) = &self.eval_col;
        assert_eq!(col_eval_write.len(), num_instances);
        assert_eq!(col_eval_read.len(), num_instances);
        let ws: F = (0..col_eval_write.len())
            .map(|i| col_eval_write[i])
            .product();
        let rs: F = (0..col_eval_read.len()).map(|i| col_eval_read[i]).product();
        assert_eq!(*col_eval_init * ws, rs * col_eval_audit);

        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_col_eval_init",
            col_eval_init,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_col_eval_read",
            col_eval_read,
        );
        <Transcript as ProofTranscript<E>>::append_scalars(
            transcript,
            b"claim_col_eval_write",
            col_eval_write,
        );
        <Transcript as ProofTranscript<E>>::append_scalar(
            transcript,
            b"claim_col_eval_audit",
            col_eval_audit,
        );

        // verify the evaluation of the sparse polynomial
        let (eval_dotp_left, eval_dotp_right) = &self.eval_val;
        assert_eq!(eval_dotp_left.len(), eval_dotp_left.len());
        assert_eq!(eval_dotp_left.len(), num_instances);
        let mut claims_dotp_circuit: Vec<F> = Vec::new();
        for i in 0..num_instances {
            assert_eq!(eval_dotp_left[i] + eval_dotp_right[i], eval[i]);

            <Transcript as ProofTranscript<E>>::append_scalar(
                transcript,
                b"claim_eval_dotp_left",
                &eval_dotp_left[i],
            );

            <Transcript as ProofTranscript<E>>::append_scalar(
                transcript,
                b"claim_eval_dotp_right",
                &eval_dotp_right[i],
            );

            claims_dotp_circuit.push(eval_dotp_left[i]);
            claims_dotp_circuit.push(eval_dotp_right[i]);
        }

        // verify the correctness of claim_row_eval_read, claim_row_eval_write, claim_col_eval_read, and claim_col_eval_write
        let mut claims_prod_circuit: Vec<F> = Vec::new();
        claims_prod_circuit.extend(row_eval_read);
        claims_prod_circuit.extend(row_eval_write);
        claims_prod_circuit.extend(col_eval_read);
        claims_prod_circuit.extend(col_eval_write);

        let (claims_ops, claims_dotp, rand_ops) = self.proof_ops.verify::<E>(
            &claims_prod_circuit,
            &claims_dotp_circuit,
            num_ops,
            transcript,
        );
        // verify the correctness of claim_row_eval_init and claim_row_eval_audit
        let (claims_mem, _claims_mem_dotp, rand_mem) = self.proof_mem.verify::<E>(
            &[
                *row_eval_init,
                *row_eval_audit,
                *col_eval_init,
                *col_eval_audit,
            ],
            &Vec::new(),
            num_cells,
            transcript,
        );
        timer.stop();

        Ok((claims_mem, rand_mem, claims_ops, claims_dotp, rand_ops))
    }
}

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
struct PolyEvalNetworkProof<E: Pairing, PC: PolyCommitmentScheme<E>> {
    proof_prod_layer: ProductLayerProof<E::ScalarField>,
    proof_hash_layer: HashLayerProof<E, PC>,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> PolyEvalNetworkProof<E, PC> {
    fn protocol_name() -> &'static [u8] {
        b"Sparse polynomial evaluation proof"
    }

    pub fn prove(
        network: &mut PolyEvalNetwork<E::ScalarField>,
        dense: &MultiSparseMatPolynomialAsDense<E::ScalarField>,
        derefs: &Derefs<E::ScalarField>,
        evals: &[E::ScalarField],
        gens: &SparseMatPolyCommitmentKey<E, PC>,
        transcript: &mut Transcript,
    ) -> Self {
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            PolyEvalNetworkProof::<E, PC>::protocol_name(),
        );

        let (proof_prod_layer, rand_mem, rand_ops) = ProductLayerProof::<E::ScalarField>::prove::<E>(
            &mut network.row_layers.prod_layer,
            &mut network.col_layers.prod_layer,
            dense,
            derefs,
            evals,
            transcript,
        );

        // proof of hash layer for row and col
        let proof_hash_layer =
            HashLayerProof::prove((&rand_mem, &rand_ops), dense, derefs, gens, transcript);

        PolyEvalNetworkProof {
            proof_prod_layer,
            proof_hash_layer,
        }
    }

    pub fn verify(
        &self,
        comm: &SparseMatPolyCommitment<E, PC>,
        comm_derefs: &DerefsCommitment<E, PC>,
        evals: &[E::ScalarField],
        gens: &SparseMatPolyCommitmentKey<E, PC>,
        rx: &[E::ScalarField],
        ry: &[E::ScalarField],
        r_mem_check: &(E::ScalarField, E::ScalarField),
        nz: usize,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        let timer = Timer::new("verify_polyeval_proof");
        <Transcript as ProofTranscript<E>>::append_protocol_name(
            transcript,
            PolyEvalNetworkProof::<E, PC>::protocol_name(),
        );

        let num_instances = evals.len();
        let (r_hash, r_multiset_check) = r_mem_check;

        let num_ops = nz.next_power_of_two();
        let num_cells = rx.len().pow2();
        assert_eq!(rx.len(), ry.len());

        let (claims_mem, rand_mem, mut claims_ops, claims_dotp, rand_ops) = self
            .proof_prod_layer
            .verify::<E>(num_ops, num_cells, evals, transcript)?;
        assert_eq!(claims_mem.len(), 4);
        assert_eq!(claims_ops.len(), 4 * num_instances);
        assert_eq!(claims_dotp.len(), 3 * num_instances);

        let (claims_ops_row, claims_ops_col) = claims_ops.split_at_mut(2 * num_instances);
        let (claims_ops_row_read, claims_ops_row_write) = claims_ops_row.split_at_mut(num_instances);
        let (claims_ops_col_read, claims_ops_col_write) = claims_ops_col.split_at_mut(num_instances);

        // verify the proof of hash layer
        self.proof_hash_layer.verify(
            (&rand_mem, &rand_ops),
            &(
                claims_mem[0],
                claims_ops_row_read.to_vec(),
                claims_ops_row_write.to_vec(),
                claims_mem[1],
            ),
            &(
                claims_mem[2],
                claims_ops_col_read.to_vec(),
                claims_ops_col_write.to_vec(),
                claims_mem[3],
            ),
            &claims_dotp,
            comm,
            gens,
            comm_derefs,
            rx,
            ry,
            r_hash,
            r_multiset_check,
            transcript,
        )?;
        timer.stop();

        Ok(())
    }
}


pub struct SparsePolyEntry<F> {
    idx: usize,
    val: F,
}

impl<F> SparsePolyEntry<F> {
    pub fn new(idx: usize, val: F) -> Self {
        SparsePolyEntry { idx, val }
    }
}

pub struct SparsePolynomial<F> {
    num_vars: usize,
    Z: Vec<SparsePolyEntry<F>>,
}

impl<F: PrimeField> SparsePolynomial<F> {
    pub fn new(num_vars: usize, Z: Vec<SparsePolyEntry<F>>) -> Self {
        SparsePolynomial { num_vars, Z }
    }

    fn compute_chi(a: &[bool], r: &[F]) -> F {
        assert_eq!(a.len(), r.len());
        let mut chi_i = F::one();
        for j in 0..r.len() {
            if a[j] {
                chi_i *= r[j];
            } else {
                chi_i *= F::one() - r[j];
            }
        }
        chi_i
    }

    // Takes O(n log n). TODO: do this in O(n) where n is the number of entries in Z
    pub fn evaluate(&self, r: &[F]) -> F {
        assert_eq!(self.num_vars, r.len());

        (0..self.Z.len())
            .map(|i| {
                let bits = self.Z[i].idx.get_bits(r.len());
                SparsePolynomial::compute_chi(&bits, r) * self.Z[i].val
            })
            .sum()
    }
}

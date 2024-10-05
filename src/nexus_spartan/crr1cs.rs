use ark_ec::CurveGroup;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{cmp::max, test_rng, One, UniformRand, Zero};
use rayon::iter::ParallelIterator;
use crate::polynomial::multilinear_poly::MultilinearPolynomial;
use super::polycommitments::{PCSKeys, PolyCommitmentScheme, VectorCommitmentScheme};
use super::{committed_relaxed_snark::CRSNARKKey, errors::R1CSError, InputsAssignment, Instance, math::Math, VarsAssignment};
#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct CRR1CSKey<E: Pairing, PC: PolyCommitmentScheme<E>> {
    pub keys: PCSKeys<E, PC>,
}

impl<E: Pairing, PC: PolyCommitmentScheme<E>> CRR1CSKey<E, PC> {
    pub fn new(SRS: &PC::SRS, num_cons: usize, num_vars: usize) -> Self {
        // Since we have commitments both to the witness and the error vectors
        // we need the commitment key to hold the larger of the two
        let n = max(num_cons, num_vars);
        CRR1CSKey {
            keys: PC::trim(SRS, n.log_2()),
        }
    }
    pub fn get_min_num_vars(num_cons: usize, num_vars: usize) -> usize {
        let n = max(num_cons, num_vars);
        n.log_2()
    }
}

#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct CRR1CSShape<F: PrimeField> {
    pub inst: Instance<F>,
}

impl<F: PrimeField> CRR1CSShape<F> {
    pub fn get_num_cons(&self) -> usize {
        self.inst.inst.get_num_cons()
    }
    pub fn get_num_vars(&self) -> usize {
        self.inst.inst.get_num_vars()
    }
    pub fn get_num_inputs(&self) -> usize {
        self.inst.inst.get_num_inputs()
    }
}

pub struct CRR1CSInstance<E: Pairing, PC: PolyCommitmentScheme<E>> {
    pub input: InputsAssignment<E::ScalarField>,
    pub u: E::ScalarField,
    pub comm_W: PC::Commitment,
    pub comm_E: PC::Commitment,
}

#[derive(Clone)]
pub struct CRR1CSWitness<F: PrimeField> {
    pub W: VarsAssignment<F>,
    pub E: Vec<F>,
}

pub fn relaxed_r1cs_is_sat<E: Pairing, PC: PolyCommitmentScheme<E>>(
    shape: &CRR1CSShape<E::ScalarField>,
    instance: &CRR1CSInstance<E, PC>,
    witness: &CRR1CSWitness<E::ScalarField>,
) -> Result<bool, R1CSError> {
    let CRR1CSWitness { W, E } = witness;
    let CRR1CSInstance { input, u, .. } = instance;
    let CRR1CSShape { inst } = shape;

    if W.assignment.len() > inst.inst.get_num_vars() {
        return Err(R1CSError::InvalidNumberOfInputs);
    }

    if input.assignment.len() != inst.inst.get_num_inputs() {
        return Err(R1CSError::InvalidNumberOfInputs);
    }

    // we might need to pad variables
    let padded_vars = {
        let num_padded_vars = inst.inst.get_num_vars();
        let num_vars = W.assignment.len();
        if num_padded_vars > num_vars {
            W.pad(num_padded_vars)
        } else {
            W.clone()
        }
    };

    // similarly we might need to pad the error vector
    let padded_E = {
        let num_padded_cons = inst.inst.get_num_cons();
        let num_cons = E.len();
        if num_padded_cons > num_cons {
            let mut padded_E = E.clone();
            padded_E.resize(num_padded_cons, E::ScalarField::zero());
            padded_E
        } else {
            E.clone()
        }
    };

    let (num_cons, num_vars, num_inputs) = (
        inst.inst.get_num_cons(),
        inst.inst.get_num_vars(),
        inst.inst.get_num_inputs(),
    );

    let z = {
        let mut z = padded_vars.assignment.to_vec();
        z.extend(&vec![*u]);
        z.extend(input.assignment.clone());
        z
    };

    // verify if Az * Bz - u * Cz = E
    let Az = inst
        .inst
        .A
        .multiply_vec(num_cons, num_vars + num_inputs + 1, &z);
    let Bz = inst
        .inst
        .B
        .multiply_vec(num_cons, num_vars + num_inputs + 1, &z);
    let Cz = inst
        .inst
        .C
        .multiply_vec(num_cons, num_vars + num_inputs + 1, &z);

    assert_eq!(Az.len(), num_cons);
    assert_eq!(Bz.len(), num_cons);
    assert_eq!(Cz.len(), num_cons);
    assert_eq!(padded_E.len(), num_cons);
    let res: usize = (0..num_cons)
        .map(|i| {
            if Az[i] * Bz[i] == *u * Cz[i] + padded_E[i] {
                0
            } else {
                1
            }
        })
        .sum();

    Ok(res == 0)
}

pub fn check_commitments<E: Pairing, PC: PolyCommitmentScheme<E>>(
    instance: &CRR1CSInstance<E, PC>,
    witness: &CRR1CSWitness<E::ScalarField>,
    key: &CRR1CSKey<E, PC>,
) -> bool {
    let CRR1CSWitness { W, E } = witness;
    let CRR1CSInstance { comm_W, comm_E, .. } = instance;

    let W = W.assignment.clone();
    let E = E.clone();

    let poly_W = MultilinearPolynomial::new(W);
    let poly_E = MultilinearPolynomial::new(E);

    let expected_comm_W = PC::commit(&poly_W, &key.keys.ck);
    let expected_comm_E = PC::commit(&poly_E, &key.keys.ck);

    expected_comm_W == *comm_W && expected_comm_E == *comm_E
}

pub fn is_sat<E: Pairing, PC: PolyCommitmentScheme<E>>(
    shape: &CRR1CSShape<E::ScalarField>,
    instance: &CRR1CSInstance<E, PC>,
    witness: &CRR1CSWitness<E::ScalarField>,
    key: &CRR1CSKey<E, PC>,
) -> Result<bool, R1CSError> {
    if !check_commitments(instance, witness, key) {
        return Ok(false);
    }
    relaxed_r1cs_is_sat(shape, instance, witness)
}

#[allow(clippy::type_complexity)]
// This produces a random satisfying structure, instance, witness, and public parameters for testing and benchmarking purposes.
pub fn produce_synthetic_crr1cs<E: Pairing, PC: PolyCommitmentScheme<E>>(
    num_cons: usize,
    num_vars: usize,
    num_inputs: usize,
) -> (
    CRR1CSShape<E::ScalarField>,
    CRR1CSInstance<E, PC>,
    CRR1CSWitness<E::ScalarField>,
    CRSNARKKey<E, PC>,
) {
    // compute random satisfying assignment for r1cs
    let (inst, vars, inputs) = Instance::produce_synthetic_r1cs(num_cons, num_vars, num_inputs);
    // the `Instance` initializer may have padded the variable lengths
    let (num_cons, num_vars, num_inputs) = (
        inst.inst.get_num_cons(),
        inst.inst.get_num_vars(),
        inst.inst.get_num_inputs(),
    );
    assert_eq!(num_vars, vars.assignment.len());
    assert_eq!(num_inputs, inputs.assignment.len());
    let shape = CRR1CSShape { inst };

    // Note that `produce_synthetic_r1cs` produces a satisfying assignment for Z = [vars, 1, inputs].
    let mut Z = vars.assignment.clone();
    Z.extend(&vec![E::ScalarField::one()]);
    Z.extend(inputs.assignment.clone());

    // Choose a random u and set Z[num_vars] = u.
    let u = E::ScalarField::rand(&mut test_rng());
    Z[num_vars] = u;

    let (poly_A, poly_B, poly_C) =
        shape
            .inst
            .inst
            .multiply_vec(num_cons, num_vars + num_inputs + 1, Z.as_slice());

    // Compute the error vector E = (AZ * BZ) - (u * CZ)
    let mut E = vec![E::ScalarField::zero(); num_cons];
    for i in 0..num_cons {
        let AB_val = poly_A[i] * poly_B[i];
        let C_val = poly_C[i];
        E[i] = AB_val - u * C_val;
    }

    // produce public parameters
    let min_num_vars =
        CRSNARKKey::<E, PC>::get_min_num_vars(num_cons, num_vars, num_inputs, num_cons);
    let SRS = PC::setup(min_num_vars, b"CRSNARK_profiler_SRS", &mut test_rng()).unwrap();
    let gens = CRSNARKKey::<E, PC>::new(&SRS, num_cons, num_vars, num_inputs, num_cons);

    // compute commitments to the vectors `vars` and `E`.
    let comm_W = <PC as VectorCommitmentScheme<E>>::commit(
        vars.assignment.as_slice(),
        &gens.gens_r1cs_sat.keys.ck,
    );
    let comm_E = <PC as VectorCommitmentScheme<E>>::commit(E.as_slice(), &gens.gens_r1cs_sat.keys.ck);
    (
        shape,
        CRR1CSInstance::<E, PC> {
            input: inputs,
            u,
            comm_W,
            comm_E,
        },
        CRR1CSWitness::<E::ScalarField> {
            W: vars.clone(),
            E: E.clone(),
        },
        gens,
    )
}
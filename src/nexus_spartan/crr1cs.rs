use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{test_rng, One};
use rand::thread_rng;
use crate::constant_for_curves::E;
use crate::kzh::KZH;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

use super::commitment_traits::{VectorCommitmentScheme};
use super::{committed_relaxed_snark::CRSNARKKey, errors::R1CSError, InputsAssignment, Instance, VarsAssignment};


#[derive(CanonicalDeserialize, CanonicalSerialize)]
pub struct CRR1CSShape<F: PrimeField + Absorb> {
    pub inst: Instance<F>,
}

impl<F: PrimeField + Absorb> CRR1CSShape<F> {
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

pub struct CRR1CSInstance<E: Pairing, PC: KZH<E>>
where
    <E as Pairing>::ScalarField: Absorb,
{
    pub input: InputsAssignment<E::ScalarField>,
    pub comm_W: PC::Commitment,
}

#[derive(Clone)]
pub struct CRR1CSWitness<E: Pairing, PC: KZH<E>>
where
    <E as Pairing>::ScalarField: Absorb,
{
    pub W: VarsAssignment<E::ScalarField>,
    pub aux_W: PC::Aux,
}

pub fn relaxed_r1cs_is_sat<E: Pairing, PC: KZH<E>>(
    shape: &CRR1CSShape<E::ScalarField>,
    instance: &CRR1CSInstance<E, PC>,
    witness: &CRR1CSWitness<E, PC>,
) -> Result<bool, R1CSError> where
    <E as Pairing>::ScalarField: Absorb,
{
    let CRR1CSWitness { W, .. } = witness;
    let CRR1CSInstance { input, .. } = instance;
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

    let (num_cons, num_vars, num_inputs) = (
        inst.inst.get_num_cons(),
        inst.inst.get_num_vars(),
        inst.inst.get_num_inputs(),
    );

    let z = {
        let mut z = padded_vars.assignment.to_vec();
        z.extend(&vec![E::ScalarField::ONE]);
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
    let res: usize = (0..num_cons)
        .map(|i| {
            if Az[i] * Bz[i] == Cz[i] {
                0
            } else {
                1
            }
        })
        .sum();

    Ok(res == 0)
}

pub fn check_commitments<E: Pairing, PC: KZH<E>>(
    instance: &CRR1CSInstance<E, PC>,
    witness: &CRR1CSWitness<E, PC>,
    key: &PC::SRS,
) -> bool where
    <E as Pairing>::ScalarField: Absorb,
{
    let CRR1CSWitness { W, .. } = witness;
    let CRR1CSInstance { comm_W, .. } = instance;

    let W = W.assignment.clone();

    let poly_W = MultilinearPolynomial::new(W);

    let (expected_comm_W, _) = PC::commit(&key, &poly_W, &mut thread_rng());

    expected_comm_W == *comm_W
}

pub fn is_sat<E: Pairing, PC: KZH<E>>(
    shape: &CRR1CSShape<E::ScalarField>,
    instance: &CRR1CSInstance<E, PC>,
    witness: &CRR1CSWitness<E, PC>,
    key: &PC::SRS,
) -> Result<bool, R1CSError> where
    <E as Pairing>::ScalarField: Absorb,
{
    if !check_commitments(instance, witness, key) {
        return Ok(false);
    }
    relaxed_r1cs_is_sat(shape, instance, witness)
}

#[allow(clippy::type_complexity)]
// This produces a random satisfying structure, instance, witness, and public parameters for testing and benchmarking purposes.
pub fn produce_synthetic_crr1cs<E: Pairing, PC: KZH<E>>(
    num_cons: usize,
    num_vars: usize,
    num_inputs: usize,
) -> (
    CRR1CSShape<E::ScalarField>,
    CRR1CSInstance<E, PC>,
    CRR1CSWitness<E, PC>,
    CRSNARKKey<E, PC>,
) where
    <E as Pairing>::ScalarField: Absorb,
{
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

    Z[num_vars] = E::ScalarField::ONE;

    // produce public parameters
    let min_num_vars = CRSNARKKey::<E, PC>::get_min_num_vars(num_cons, num_vars, num_inputs);
    let SRS = PC::setup(min_num_vars, &mut test_rng());
    let gens = CRSNARKKey::<E, PC>::new(&SRS, num_cons, num_vars, num_inputs, num_cons);

    // compute commitments to the vectors `vars` and `E`.
    let comm_W = <PC as VectorCommitmentScheme<E>>::commit(
        vars.assignment.as_slice(),
        &gens.gens_r1cs_sat,
    );

    (
        shape,
        CRR1CSInstance::<E, PC> {
            input: inputs,
            comm_W: comm_W.0,
        },
        CRR1CSWitness::<E, PC> {
            W: vars.clone(),
            aux_W: comm_W.1,
        },
        gens,
    )
}
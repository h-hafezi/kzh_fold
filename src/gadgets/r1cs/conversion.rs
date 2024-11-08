// we do conversion from a constraint system into the R1CS abstract used by Nexus

use crate::commitment::{CommitmentScheme, Len};
use crate::gadgets::r1cs::r1cs::{commit_T, R1CSInstance, R1CSWitness, RelaxedR1CSInstance, RelaxedR1CSWitness};
use crate::gadgets::r1cs::R1CSShape;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef, SynthesisMode};
use rand::thread_rng;

/// given a ConstraintSystem object, returns corresponding R1CSShape, R1CSInstance and R1CSWitness
pub fn convert_constraint_system_into_instance_witness<F, C1, G1>(
    cs: ConstraintSystemRef<F>,
    pp: &C1::PP,
) -> (R1CSShape<G1>, R1CSInstance<Projective<G1>, C1>, R1CSWitness<Projective<G1>>)
where
    F: PrimeField + Absorb,
    C1: CommitmentScheme<Projective<G1>>,
    G1: SWCurveConfig<ScalarField=F>,
{
    // the cs has to be finalised before the matrices A, B and C are generated
    cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
    cs.finalize();

    // make sure cs is already satisfied
    assert!(cs.is_satisfied().unwrap(), "the passed constraint system isn't satisfied");

    let shape = R1CSShape::from(cs.clone());
    let cs_borrow = cs.borrow().unwrap();
    let W = cs_borrow.witness_assignment.clone();
    let X = cs_borrow.instance_assignment.clone();

    // make sure pp has the right length
    assert_eq!(W.len(), pp.len(), "the length of the witness and commitment public params are inconsistent");

    let w = R1CSWitness { W };

    let commitment_W = w.commit::<C1>(pp);
    let u = R1CSInstance { commitment_W, X };

    assert!(shape.is_satisfied(&u, &w, pp).is_ok());

    (shape, u, w)
}

/// Generates a random constraint system with specified parameters.
///
/// # Arguments
/// * `num_constraints` - The number of constraints to generate.
/// * `num_vars` - The number of witness variables.
/// * `num_io` - The number of public input/output variables.
pub fn generate_random_constraint_system<F: PrimeField>(
    num_constraints: usize,
    num_vars: usize,
    num_io: usize,
) -> ConstraintSystemRef<F> {
    assert!(num_vars >= num_constraints);
    let cs = ConstraintSystem::<F>::new_ref();
    let rng = &mut thread_rng();

    // Create random variables for the public inputs, the first one is always one
    for _ in 0..num_io - 1 {
        let _ = FpVar::new_input(cs.clone(), || Ok(F::rand(rng))).unwrap();
    }

    // Create random variables for the witness
    let mut witness_vars = vec![];
    for _ in 0..num_vars - num_constraints {
        let var = FpVar::new_witness(cs.clone(), || Ok(F::rand(rng))).unwrap();
        witness_vars.push(var);
    }

    // This one is important because it affects the shape of A, B, C which has to be the same for all instance/witnesses we create
    for i in 0..num_constraints {
        // define j in case num_constraint > num_vars, otherwise it wouldn't make a difference
        let j = i % (num_vars - num_constraints);
        let _ = &witness_vars[j] * &witness_vars[j];
    }

    assert_eq!(cs.num_constraints(), num_constraints);
    assert_eq!(cs.num_witness_variables(), num_vars);
    assert_eq!(cs.num_instance_variables(), num_io);

    cs
}

/// outputs a new instance/witness following the shape of generate_random_constraint_system
pub fn get_random_r1cs_instance_witness<F, C1, G1>(num_constraints: usize,
                                                   num_vars: usize,
                                                   num_io: usize,
                                                   pp: &C1::PP,
) -> (R1CSShape<G1>, R1CSInstance<Projective<G1>, C1>, R1CSWitness<Projective<G1>>) where
    F: PrimeField + Absorb,
    C1: CommitmentScheme<Projective<G1>>,
    G1: SWCurveConfig<ScalarField=F>,
{
    assert!(num_vars >= num_constraints);

    // generate a constraint system corresponding the shape
    let cs = generate_random_constraint_system(num_constraints, num_vars, num_io);

    // generate corresponding instance/witness
    let (shape, instance, witness) = convert_constraint_system_into_instance_witness(cs.clone(), pp);

    assert_eq!(cs.num_constraints(), shape.num_constraints);
    assert_eq!(cs.num_witness_variables(), shape.num_vars);
    assert_eq!(cs.num_instance_variables(), shape.num_io);
    assert_eq!(cs.num_constraints(), num_constraints);
    assert_eq!(cs.num_witness_variables(), num_vars);
    assert_eq!(cs.num_instance_variables(), num_io);

    (shape, instance, witness)
}

pub fn get_random_relaxed_r1cs_instance_witness<F, C1, G1>(num_constraints: usize,
                                                           num_vars: usize,
                                                           num_io: usize,
                                                           pp: &C1::PP,
) -> (R1CSShape<G1>, RelaxedR1CSInstance<Projective<G1>, C1>, RelaxedR1CSWitness<Projective<G1>>) where
    F: PrimeField + Absorb,
    C1: CommitmentScheme<Projective<G1>>,
    G1: SWCurveConfig<ScalarField=F>,
{
    assert!(num_vars >= num_constraints);

    // generate a constraint system corresponding the shape
    let cs1 = generate_random_constraint_system(num_constraints, num_vars, num_io);
    let cs2 = generate_random_constraint_system(num_constraints, num_vars, num_io);

    // generate corresponding instance/witness
    let (shape1, instance1, witness1) = convert_constraint_system_into_instance_witness(cs1.clone(), pp);
    let (shape2, instance2, witness2) = convert_constraint_system_into_instance_witness(cs2.clone(), pp);

    assert_eq!(shape1, shape2);

    // make instance1/witness1 into relaxed
    let relaxed_instance = RelaxedR1CSInstance::from(&instance1);
    let relaxed_witness = RelaxedR1CSWitness::from_r1cs_witness(&shape1, &witness1);

    shape1.is_relaxed_satisfied(&relaxed_instance, &relaxed_witness, pp).expect("not satisfied r1cs instance");

    // fold the two instance/witness

    let (t, com_t) = commit_T(&shape1, pp, &relaxed_instance, &relaxed_witness, &instance2, &witness2).unwrap();
    let folding_randomness = F::rand(&mut thread_rng());

    let folded_instance = relaxed_instance.fold(&instance2, &com_t, &folding_randomness).unwrap();
    let folded_witness = relaxed_witness.fold(&witness2, &t, &folding_randomness).unwrap();

    shape1.is_relaxed_satisfied(&folded_instance, &folded_witness, pp).expect("not satisfied r1cs instance");

    (shape1, folded_instance, folded_witness)
}


#[cfg(test)]
mod tests {
    use crate::commitment::CommitmentScheme;
    use crate::constant_for_curves::{G1Projective, ScalarField, G1};
    use crate::gadgets::r1cs::conversion::{convert_constraint_system_into_instance_witness, generate_random_constraint_system, get_random_r1cs_instance_witness, get_random_relaxed_r1cs_instance_witness};
    use crate::hash::pederson::PedersenCommitment;
    use ark_ec::short_weierstrass::Affine;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::eq::EqGadget;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::FieldVar;
    use ark_relations::r1cs::ConstraintSystem;

    type F = ScalarField;
    type C1 = PedersenCommitment<G1Projective>;

    #[test]
    fn test_convert_constraint_system_into_instance_witness() {
        let cs = ConstraintSystem::<F>::new_ref();

        // add some random constraints
        for i in 0u8..10 {
            let x = FpVar::new_witness(cs.clone(), || Ok(F::from(i))).unwrap();
            let y = FpVar::new_witness(cs.clone(), || Ok(F::from(i + 1))).unwrap();

            y.enforce_equal(&(x + FpVar::one())).unwrap()
        }

        assert!(cs.is_satisfied().unwrap());

        let pp: Vec<Affine<G1>> = C1::setup(cs.num_witness_variables(), b"test", &());

        let (shape, u, w) = convert_constraint_system_into_instance_witness::<F, C1, G1>(cs.clone(), &pp);

        // assert that it's in fact satisfied
        shape.is_satisfied(&u, &w, &pp).expect("cs is not satisfied");
    }

    #[test]
    fn test_generate_random_constraint_system() {
        let (num_constraints, num_io, num_vars) = (10, 3, 17);
        let cs = generate_random_constraint_system::<F>(num_constraints, num_vars, num_io);
        cs.finalize();

        let pp: Vec<Affine<G1>> = C1::setup(num_vars, b"test", &());

        let (shape, u, w) = convert_constraint_system_into_instance_witness::<F, C1, G1>(cs.clone(), &pp);

        // assert that it's in fact satisfied
        shape.is_satisfied(&u, &w, &pp).expect("cs is not satisfied");

        // assert the shape is well formatted
        assert_eq!((shape.num_io, shape.num_vars, shape.num_constraints), (num_io, num_vars, num_constraints));
    }

    #[test]
    fn test_random_r1cs() {
        let (num_constraints, num_io, num_vars) = (10, 3, 17);
        let pp: Vec<Affine<G1>> = C1::setup(num_vars, b"test", &());

        let (shape, instance, witness) = get_random_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &pp);

        // assert it's satisfied
        shape.is_satisfied(&instance, &witness, &pp).expect("unsatisfied r1cs");

        // generate a relaxed instance/witness this time
        let (relaxed_shape, relaxed_instance, relaxed_witness) = get_random_relaxed_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &pp);

        // assert the shape is equal to the previous shape
        assert_eq!(shape, relaxed_shape);

        // make sure the instance is satisfied
        shape.is_relaxed_satisfied(&relaxed_instance, &relaxed_witness, &pp).expect("unsatisfied r1cs");
    }
}
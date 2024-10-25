use crate::accumulation_circuit::affine_to_projective;
use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::ova::commit_T;
use crate::gadgets::r1cs::{OvaInstance, OvaWitness, R1CSShape, RelaxedOvaInstance, RelaxedOvaWitness};
use crate::hash::pederson::PedersenCommitment;
use crate::nova::cycle_fold::coprocessor::{setup_shape, synthesize, SecondaryCircuit};
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::{Field, PrimeField};
use ark_std::UniformRand;
use rand::Rng;
use std::marker::PhantomData;

pub struct SignatureVerifierProver<G1, G2, C2, E>
where
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=G1::ScalarField>,
{
    /// the randomness used for taking linear combination, it should be input fiat-shamir
    pub beta: G2::ScalarField,

    /// pk_t = pk_1 + pk_2
    pub pk_1: E::G1Affine,
    pub pk_2: E::G1Affine,
    pub pk_t: E::G1Affine,

    /// running cycle fold instance
    pub shape: R1CSShape<G2>,
    pub commitment_pp: <C2 as CommitmentScheme<Projective<G2>>>::PP,
    pub cycle_fold_running_instance: RelaxedOvaInstance<G2, C2>,
    pub cycle_fold_running_witness: RelaxedOvaWitness<G2>,

    pub phantom: PhantomData<G1>,
}

impl<G1, G2, C2, E> SignatureVerifierProver<G1, G2, C2, E>
where
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=G1::ScalarField>,
{
    pub fn rand<R: Rng>(rng: &mut R, beta: G2::ScalarField) -> Self {
        let (pk_1, pk_2, pk_t) = SignatureVerifierProver::<G1, G2, C2, E>::get_satisfying_public_keys(rng);
        let shape = setup_shape::<G1, G2>().unwrap();
        let commitment_pp: Vec<Affine<G2>> = PedersenCommitment::<Projective<G2>>::setup(shape.num_vars + shape.num_constraints, b"test", &());
        let (cycle_fold_running_instance, cycle_fold_running_witness) = SignatureVerifierProver::<G1, G2, C2, E>::get_satisfying_running_instance_witness(&shape, &commitment_pp, rng);
        SignatureVerifierProver {
            beta,
            pk_1,
            pk_2,
            pk_t,
            shape,
            commitment_pp,
            cycle_fold_running_instance,
            cycle_fold_running_witness,
            phantom: Default::default(),
        }
    }

    pub fn get_satisfying_public_keys<R: Rng>(rng: &mut R) -> (E::G1Affine, E::G1Affine, E::G1Affine) {
        let (pk_1, pk_2) = (E::G1Affine::rand(rng), E::G1Affine::rand(rng));

        (pk_1, pk_2, (pk_1 + pk_2).into())
    }

    pub fn get_auxiliary_input_for_public_keys(shape: &R1CSShape<G2>, commitment_pp: &C2::PP, pk_1: E::G1Affine, pk_2: E::G1Affine) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        assert_eq!(shape.num_constraints + shape.num_vars, commitment_pp.len());
        let pk1 = affine_to_projective(pk_1.clone());
        let pk2 = affine_to_projective(pk_2.clone());

        // pk_t =  pk_1 + pk_2
        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1: pk1,
            g2: pk2,
            g_out: pk1 + pk2,
            r: G2::ScalarField::ONE,
            flag: true,
        }, &commitment_pp[0..shape.num_vars].to_vec(),
        ).unwrap()
    }

    pub fn get_satisfying_running_instance_witness<R: Rng>(shape: &R1CSShape<G2>, commitment_pp: &C2::PP, rng: &mut R) -> (
        RelaxedOvaInstance<G2, C2>,
        RelaxedOvaWitness<G2>
    ) {
        // check the length of the commitment keys
        assert_eq!(shape.num_constraints + shape.num_vars, commitment_pp.len());

        // get a random satisfying instance/witness
        let (instance, witness) = SignatureVerifierProver::<G1, G2, C2, E>::get_auxiliary_input_for_public_keys(
            shape,
            commitment_pp,
            E::G1Affine::rand(rng),
            E::G1Affine::rand(rng),
        );

        // convert to relaxed ova instance/witness
        let relaxed_ova_instance = RelaxedOvaInstance::from(&instance);
        let relaxed_ova_witness = RelaxedOvaWitness::from(&shape, &witness);

        (relaxed_ova_instance, relaxed_ova_witness)
    }

    pub fn compute_cycle_fold_proofs_and_final_instance(&self, instance: &OvaInstance<G2, C2>, witness: &OvaWitness<G2>) -> (
        C2::Commitment,
        RelaxedOvaInstance<G2, C2>,
        RelaxedOvaWitness<G2>
    ) {
        let (T, com_T) = commit_T(
            &self.shape,
            &self.commitment_pp[self.shape.num_vars..].to_vec(),
            &self.cycle_fold_running_instance,
            &self.cycle_fold_running_witness,
            instance,
            witness,
        ).unwrap();

        // Fold the running instance and witness with the first proof
        let new_running_instance = self.cycle_fold_running_instance.fold(instance, &com_T, &self.beta).unwrap();
        let new_running_witness = self.cycle_fold_running_witness.fold(witness, &T, &self.beta).unwrap();

        (com_T, new_running_instance, new_running_witness)
    }
}

#[cfg(test)]
mod test {
    use crate::commitment::CommitmentScheme;
    use crate::constant_for_curves::{BaseField, G1Affine, E, G1, G2};
    use crate::hash::pederson::PedersenCommitment;
    use crate::nova::cycle_fold::coprocessor::setup_shape;
    use crate::signature_aggregation::verifier_circuit::prover::SignatureVerifierProver;
    use ark_ec::short_weierstrass::{Affine, Projective};
    use ark_std::UniformRand;
    use rand::thread_rng;

    type GrumpkinCurveGroup = ark_grumpkin::Projective;
    type C2 = PedersenCommitment<GrumpkinCurveGroup>;
    type Q = BaseField;

    #[test]
    fn test_satisfying_public_key() {
        let rng = &mut thread_rng();
        let (pk_1, pk_2, pk_t): (G1Affine, G1Affine, G1Affine) = SignatureVerifierProver::<G1, G2, C2, E>::get_satisfying_public_keys(rng);

        assert_eq!(pk_1 + pk_2, pk_t)
    }

    #[test]
    fn test_get_auxiliary_input_for_public_keys() {
        let shape = setup_shape::<G1, G2>().unwrap();
        let commitment_pp: Vec<Affine<G2>> = PedersenCommitment::<Projective<G2>>::setup(shape.num_vars + shape.num_constraints, b"test", &());
        let rng = &mut thread_rng();
        let (pk_1, pk_2, pk_t): (G1Affine, G1Affine, G1Affine) = SignatureVerifierProver::<G1, G2, C2, E>::get_satisfying_public_keys(rng);
        let (instance, witness) = SignatureVerifierProver::<G1, G2, C2, E>::get_auxiliary_input_for_public_keys(&shape, &commitment_pp, pk_1, pk_2);
        shape.is_ova_satisfied(&instance, &witness, &commitment_pp).unwrap()
    }

    #[test]
    fn test_get_satisfying_running_instance_witness() {
        let shape = setup_shape::<G1, G2>().unwrap();
        let commitment_pp: Vec<Affine<G2>> = PedersenCommitment::<Projective<G2>>::setup(shape.num_vars + shape.num_constraints, b"test", &());
        let (relaxed_instance, relaxed_witness) = SignatureVerifierProver::<G1, G2, C2, E>::get_satisfying_running_instance_witness(&shape, &commitment_pp, &mut thread_rng());
        shape.is_relaxed_ova_satisfied(&relaxed_instance, &relaxed_witness, &commitment_pp).unwrap()
    }

    #[test]
    fn test_compute_cycle_fold_proofs_and_final_instance() {
        let rng = &mut thread_rng();

        // get a random prover
        let prover = SignatureVerifierProver::<G1, G2, C2, E>::rand(rng, Q::from(2u8));

        // get a set of random public keys
        let (pk_1, pk_2, pk_t): (G1Affine, G1Affine, G1Affine) = SignatureVerifierProver::<G1, G2, C2, E>::get_satisfying_public_keys(rng);

        // get corresponding instance/witness of (pk_1, pk_2, pk_t)
        let (instance, witness) = SignatureVerifierProver::<G1, G2, C2, E>::get_auxiliary_input_for_public_keys(&prover.shape, &prover.commitment_pp, pk_1, pk_2);
        prover.shape.is_ova_satisfied(&instance, &witness, &prover.commitment_pp).unwrap();

        // fold it with the prover's running instance/witness
        let (com_T, new_running_instance, new_running_witness) = prover.compute_cycle_fold_proofs_and_final_instance(&instance, &witness);

        // check validity
        prover.shape.is_relaxed_ova_satisfied(&new_running_instance, &new_running_witness, &prover.commitment_pp).unwrap();
    }
}
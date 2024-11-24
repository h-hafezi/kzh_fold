use crate::accumulation_circuit::affine_to_projective;
use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::ova::commit_T;
use crate::gadgets::r1cs::{OvaInstance, OvaWitness, R1CSShape, RelaxedOvaInstance, RelaxedOvaWitness};
use crate::hash::pederson::PedersenCommitment;
use crate::nova::cycle_fold::coprocessor::{setup_shape, synthesize, SecondaryCircuit};
use crate::signature_aggregation::signature_aggregation::SignatureAggrData;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::CurveConfig;
use ark_ff::{Field, PrimeField};
use ark_std::UniformRand;
use rand::Rng;
use std::marker::PhantomData;
use ark_crypto_primitives::sponge::Absorb;
use crate::gadgets::non_native::util::cast_field;

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
    pub running_pk: E::G1Affine,
    pub current_pk: E::G1Affine,

    /// commitments to bitfields
    pub com_bitfield_C: E::G1Affine,
    pub com_bitfield_B_1: E::G1Affine,
    pub com_bitfield_B_2: E::G1Affine,

    /// running cycle fold instance
    pub ova_shape: R1CSShape<G2>,
    pub ova_commitment_pp: <C2 as CommitmentScheme<Projective<G2>>>::PP,
    pub ova_running_instance: RelaxedOvaInstance<G2, C2>,
    pub ova_running_witness: RelaxedOvaWitness<G2>,

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
    // given a signature_aggregate_data, it returns a random prover
    pub fn rand<R: Rng>(rng: &mut R, signature_aggregate_data: SignatureAggrData<E, G1::ScalarField>, beta: G2::ScalarField) -> Self
    where
        <G2 as CurveConfig>::ScalarField: Absorb,
        <G2 as CurveConfig>::BaseField: Absorb
    {
        // get ova shape
        let ova_shape = setup_shape::<G1, G2>().unwrap();

        // get ova commitment public parameters according to the shape
        let ova_commitment_pp: Vec<Affine<G2>> = PedersenCommitment::<Projective<G2>>::setup(
            ova_shape.num_vars + ova_shape.num_constraints,
            b"test",
            &(),
        );

        SignatureVerifierProver {
            beta,
            running_pk: E::G1Affine::rand(rng),
            current_pk: signature_aggregate_data.pk,
            com_bitfield_C: signature_aggregate_data.bitfield_commitment.C,
            com_bitfield_B_1: signature_aggregate_data.B_1_commitment.C,
            com_bitfield_B_2: signature_aggregate_data.B_2_commitment.C,
            ova_commitment_pp,
            ova_running_instance: RelaxedOvaInstance::new(&ova_shape),
            ova_running_witness: RelaxedOvaWitness::zero(&ova_shape),
            phantom: Default::default(),
            ova_shape,
        }
    }

    // get ova auxiliary input for pk
    pub fn get_ova_auxiliary_input_pk(
        &self,
    ) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        assert_eq!(self.ova_shape.num_constraints + self.ova_shape.num_vars, self.ova_commitment_pp.len());
        let running_pk = affine_to_projective(self.running_pk.clone());
        let current_pk = affine_to_projective(self.current_pk.clone());

        // final_pk =  running_pk + current_pk
        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1: running_pk,
            g2: current_pk,
            g_out: running_pk + current_pk,
            r: G2::ScalarField::ONE,
            flag: true,
        }, &self.ova_commitment_pp[0..self.ova_shape.num_vars].to_vec(),
        ).unwrap()
    }

    // get ova auxiliary input for B1 + B2 * c_0
    pub fn get_ova_auxiliary_input_bitfield_1(
        &self,
        r: G1::ScalarField,
    ) -> (Projective<G1>, (OvaInstance<G2, C2>, OvaWitness<G2>)) {
        assert_eq!(self.ova_shape.num_constraints + self.ova_shape.num_vars, self.ova_commitment_pp.len());
        let B1 = affine_to_projective(self.com_bitfield_B_1.clone());
        let B2 = affine_to_projective(self.com_bitfield_B_2.clone());

        // temp = B1 + B2 * c_0
        (B1 + (B2 * r), synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1: B2,
            g2: B1,
            g_out: B1 + (B2 * r),
            r: cast_field::<G1::ScalarField, G1::BaseField>(r),
            flag: true,
        }, &self.ova_commitment_pp[0..self.ova_shape.num_vars].to_vec(),
        ).unwrap())
    }

    // get ova auxiliary input for res = temp + C * c_1
    pub fn get_ova_auxiliary_input_bitfield_2(
        &self,
        temp: Projective<G1>,
        r: G1::ScalarField,
    ) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        assert_eq!(self.ova_shape.num_constraints + self.ova_shape.num_vars, self.ova_commitment_pp.len());
        let C = affine_to_projective(self.com_bitfield_C.clone());

        // pk_t =  C * r + temp
        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1: C,
            g2: temp,
            g_out: temp + (C * r),
            r: cast_field::<G1::ScalarField, G1::BaseField>(r),
            flag: true,
        }, &self.ova_commitment_pp[0..self.ova_shape.num_vars].to_vec(),
        ).unwrap()
    }

    pub fn compute_ova_final_instance(&self, (c_0, c_1): (G1::ScalarField, G1::ScalarField)) -> (
        RelaxedOvaInstance<G2, C2>,
        RelaxedOvaWitness<G2>,
        (
            C2::Commitment,
            C2::Commitment,
            C2::Commitment,
        )
    ) {
        // *********************************** fold pk ***********************************
        let (instance_pk, witness_pk) = self.get_ova_auxiliary_input_pk();

        let (cross_term_error_pk, cross_term_error_commitment_pk) = commit_T(
            &self.ova_shape,
            &self.ova_commitment_pp[self.ova_shape.num_vars..].to_vec(),
            &self.ova_running_instance,
            &self.ova_running_witness,
            &instance_pk,
            &witness_pk,
        ).unwrap();

        // Fold the running instance and witness with the first proof
        let folded_instance_1 = self.ova_running_instance.fold(
            &instance_pk,
            &cross_term_error_commitment_pk,
            &self.beta,
        ).unwrap();

        let folded_witness_1 = self.ova_running_witness.fold(
            &witness_pk,
            &cross_term_error_pk,
            &self.beta,
        ).unwrap();

        // *********************************** fold bitfield_1 ***********************************
        let (temp, (instance_bitfield_1, witness_bitfield_1)) = self.get_ova_auxiliary_input_bitfield_1(c_0);

        let (cross_term_error_bitfield_1, cross_term_error_commitment_bitfield_1) = commit_T(
            &self.ova_shape,
            &self.ova_commitment_pp[self.ova_shape.num_vars..].to_vec(),
            &folded_instance_1,
            &folded_witness_1,
            &instance_bitfield_1,
            &witness_bitfield_1,
        ).unwrap();

        // Fold the running instance and witness with the first proof
        let folded_instance_2 = folded_instance_1.fold(
            &instance_bitfield_1,
            &cross_term_error_commitment_bitfield_1,
            &self.beta,
        ).unwrap();

        let folded_witness_2 = folded_witness_1.fold(
            &witness_bitfield_1,
            &cross_term_error_bitfield_1,
            &self.beta,
        ).unwrap();

        // *********************************** fold bitfield_2 ***********************************
        let (instance_bitfield_2, witness_bitfield_2) = self.get_ova_auxiliary_input_bitfield_2(temp, c_1);

        let (cross_term_error_bitfield_2, cross_term_error_commitment_bitfield_2) = commit_T(
            &self.ova_shape,
            &self.ova_commitment_pp[self.ova_shape.num_vars..].to_vec(),
            &folded_instance_2,
            &folded_witness_2,
            &instance_bitfield_2,
            &witness_bitfield_2,
        ).unwrap();

        // Fold the running instance and witness with the first proof
        let final_instance = folded_instance_2.fold(
            &instance_bitfield_2,
            &cross_term_error_commitment_bitfield_2,
            &self.beta,
        ).unwrap();

        let final_witness = folded_witness_2.fold(
            &witness_bitfield_2,
            &cross_term_error_bitfield_2,
            &self.beta,
        ).unwrap();

        (
            final_instance, final_witness,
            (
                cross_term_error_commitment_pk,
                cross_term_error_commitment_bitfield_1,
                cross_term_error_commitment_bitfield_2
            )
        )
    }
}


#[cfg(test)]
mod test {
    use std::ops::Mul;
    use ark_crypto_primitives::sponge::Absorb;
    use ark_ec::AffineRepr;
    use ark_ec::pairing::Pairing;
    use crate::constant_for_curves::{BaseField, ScalarField, C2, E, G1, G2};
    use crate::signature_aggregation::verifier_circuit::prover::SignatureVerifierProver;
    use ark_ff::{Field, PrimeField};
    use ark_std::UniformRand;
    use rand::{thread_rng, RngCore};
    use crate::kzh_fold::kzh2_fold::Accumulator2;
    use crate::kzh::KZH;
    use crate::kzh::kzh2::KZH2;
    use crate::signature_aggregation::signature_aggregation::{SignatureAggrData, SignatureAggrSRS};
    use crate::math::Math;

    type Q = BaseField;
    type F = ScalarField;

    impl<E: Pairing> SignatureAggrSRS<E>
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
        <E as Pairing>::ScalarField: Absorb,
    {
        pub(crate) fn new<R: RngCore>(degree_x: usize, degree_y: usize, rng: &mut R) -> Self {
            let pcs_srs = KZH2::setup((degree_x * degree_y).log_2(), rng);

            SignatureAggrSRS {
                acc_srs: Accumulator2::setup(pcs_srs, rng),
            }
        }
    }

    fn get_random_prover() -> SignatureVerifierProver<G1, G2, C2, E> {
        let rng = &mut thread_rng();
        let signature_aggregation_data = {
            let degree_x = 64usize;
            let degree_y = 64usize;
            let num_vars = 12usize;
            let srs = SignatureAggrSRS::<E>::new(degree_x, degree_y, rng);
            SignatureAggrData::rand(num_vars, &srs.acc_srs, rng)
        };
        let beta = Q::rand(rng);
        let prover = SignatureVerifierProver::rand(rng, signature_aggregation_data, beta);

        prover
    }

    #[test]
    fn test_get_auxiliary_input_for_public_keys() {
        let prover = get_random_prover();
        let (instance, witness) = prover.get_ova_auxiliary_input_pk();

        prover.ova_shape.is_ova_satisfied(
            &instance,
            &witness,
            &prover.ova_commitment_pp
        ).unwrap();

        let secondary_circuit_pk = instance.parse_secondary_io::<G1>().unwrap();
        assert_eq!(secondary_circuit_pk.g_out, prover.running_pk + prover.current_pk);
        assert_eq!(secondary_circuit_pk.r, Q::ONE);
    }

    #[test]
    fn test_get_auxiliary_input_for_bitfield() {
        let rng = &mut thread_rng();
        let prover = get_random_prover();
        let (c_0, c_1) = (F::rand(rng), F::rand(rng));
        let (temp, (instance_1, witness_1)) = prover.get_ova_auxiliary_input_bitfield_1(c_0);

        prover.ova_shape.is_ova_satisfied(
            &instance_1,
            &witness_1,
            &prover.ova_commitment_pp
        ).unwrap();

        let secondary_circuit_bitfield_1 = instance_1.parse_secondary_io::<G1>().unwrap();
        assert_eq!(secondary_circuit_bitfield_1.g1, prover.com_bitfield_B_2);
        assert_eq!(secondary_circuit_bitfield_1.g2, prover.com_bitfield_B_1);
        assert_eq!(secondary_circuit_bitfield_1.g_out, temp);
        assert_eq!(secondary_circuit_bitfield_1.flag, true);

        let (instance_2, witness_2) = prover.get_ova_auxiliary_input_bitfield_2(temp, c_1);

        prover.ova_shape.is_ova_satisfied(
            &instance_2,
            &witness_2,
            &prover.ova_commitment_pp
        ).unwrap();

        let secondary_circuit_bitfield_2 = instance_2.parse_secondary_io::<G1>().unwrap();
        assert_eq!(secondary_circuit_bitfield_2.g1, prover.com_bitfield_C);
        assert_eq!(secondary_circuit_bitfield_2.g2, temp);
        assert_eq!(secondary_circuit_bitfield_2.g_out, prover.com_bitfield_B_1 + prover.com_bitfield_B_2.mul(c_0) + prover.com_bitfield_C.mul(c_1));
        assert_eq!(secondary_circuit_bitfield_2.flag, true);
    }
}

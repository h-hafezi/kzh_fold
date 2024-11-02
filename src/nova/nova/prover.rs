use crate::accumulation_circuit::affine_to_projective;
use crate::commitment::{Commitment, CommitmentScheme};
use crate::gadgets::r1cs::{OvaInstance, OvaWitness, R1CSInstance, R1CSShape, R1CSWitness, RelaxedOvaInstance, RelaxedOvaWitness, RelaxedR1CSInstance, RelaxedR1CSWitness};
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::{AffineRepr, CurveConfig};
use ark_ff::PrimeField;
use crate::gadgets::absorb::relaxed_r1cs_instance_to_sponge_vector;
use crate::gadgets::non_native::util::convert_field_one_to_field_two;
use crate::gadgets::r1cs::r1cs::commit_T;
use crate::nova::cycle_fold::coprocessor::{synthesize, SecondaryCircuit};

pub struct NovaProver<F, G1, G2, C1, C2>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    pub initial_transcript: Transcript<F>,
    pub final_transcript: Transcript<F>,

    /// shape of the main curve
    pub shape: R1CSShape<G1>,

    /// srs for the accumulation
    pub commitment_pp: <C1 as CommitmentScheme<Projective<G1>>>::PP,

    /// the instance to be folded
    pub current_accumulator: (R1CSInstance<G1, C1>, R1CSWitness<G1>),

    /// the running accumulator
    pub running_accumulator: (RelaxedR1CSInstance<G1, C1>, RelaxedR1CSWitness<G1>),

    /// running cycle fold instance
    pub ova_shape: R1CSShape<G2>,
    pub ova_commitment_pp: <C2 as CommitmentScheme<Projective<G2>>>::PP,
    pub cycle_fold_running_instance: RelaxedOvaInstance<G2, C2>,
    pub cycle_fold_running_witness: RelaxedOvaWitness<G2>,
}

impl<F, G1, G2, C1, C2> NovaProver<F, G1, G2, C1, C2>
where
    F: PrimeField + Absorb,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F>,
    G2::BaseField: PrimeField,
    C1: CommitmentScheme<Projective<G1>, PP=Vec<Affine<G1>>, Commitment = Projective<G1>>,
    C2: CommitmentScheme<Projective<G2>, PP=Vec<Affine<G2>>>,
    <G2 as CurveConfig>::ScalarField: Absorb,
{
    pub fn compute_cross_term_error(&self) -> Projective<G1> {
        let (_, com_t) = commit_T(&self.shape, &self.commitment_pp, &self.running_accumulator.0, &self.running_accumulator.1, &self.current_accumulator.0, &self.current_accumulator.1).unwrap();

        com_t
    }

    pub fn compute_beta(&mut self) -> F {
        self.initial_transcript.append_scalars(b"label", relaxed_r1cs_instance_to_sponge_vector(&self.running_accumulator.0).as_slice());
        self.initial_transcript.append_scalars(b"label", relaxed_r1cs_instance_to_sponge_vector(&self.current_accumulator.0).as_slice());

        let beta = self.initial_transcript.challenge_scalar("challenge");

        beta
    }

    pub fn compute_final_accumulator(&self, beta: F) -> (RelaxedR1CSInstance<G1, C1>, RelaxedR1CSWitness<G1>) {
        let (t, com_t) = commit_T(&self.shape, &self.commitment_pp, &self.running_accumulator.0, &self.running_accumulator.1, &self.current_accumulator.0, &self.current_accumulator.1).unwrap();

        // folding two instances
        let folded_instance = self.running_accumulator.0.fold(&self.current_accumulator.0, &com_t, &beta).unwrap();

        // folding two witnesses
        let folded_witness = self.running_accumulator.1.fold(&self.current_accumulator.1, &t, &beta).unwrap();

        (folded_instance, folded_witness)
    }

    pub fn compute_auxiliary_input_T(&self, beta: F) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let g1 = self.compute_cross_term_error();
        let g2 = affine_to_projective(self.running_accumulator.0.commitment_E.into());

        let g_out = (g1 * beta) + g2 ;

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1,
            g2,
            g_out,
            r: convert_field_one_to_field_two::<G1::ScalarField, G1::BaseField>(beta),
            flag: true,
        }, &self.ova_commitment_pp[0..self.ova_shape.num_vars].to_vec()).unwrap()
    }

    pub fn compute_auxiliary_input_W(&self, beta: F) -> (OvaInstance<G2, C2>, OvaWitness<G2>) {
        let g1 = affine_to_projective(self.current_accumulator.0.commitment_W.into());
        let g2 = affine_to_projective(self.running_accumulator.0.commitment_W.into());
        let g_out = (g1 * beta) + g2;

        synthesize::<G1, G2, C2>(SecondaryCircuit {
            g1,
            g2,
            g_out,
            r: convert_field_one_to_field_two::<G1::ScalarField, G1::BaseField>(beta),
            flag: true,
        }, &self.ova_commitment_pp[0..self.ova_shape.num_vars].to_vec()).unwrap()
    }
}

#[cfg(test)]
mod test {
    use ark_ec::short_weierstrass::{Affine, Projective};
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::{G1Projective, G2Projective, ScalarField, C1, C2, G1, G2};
    use crate::hash::pederson::PedersenCommitment;
    use crate::nova::cycle_fold::coprocessor::setup_shape;
    use crate::nova::nova::prover::NovaProver;
    use crate::transcript::transcript::Transcript;
    use crate::commitment::CommitmentScheme;
    use crate::gadgets::r1cs::{RelaxedOvaInstance, RelaxedOvaWitness};
    use crate::gadgets::r1cs::conversion::{get_random_r1cs_instance_witness, get_random_relaxed_r1cs_instance_witness};

    type F = ScalarField;

    #[test]
    fn test_compute_final_accumulator() {
        // the shape of the secondary curve R1CS instance
        let ova_shape = setup_shape::<G1, G2>().unwrap();
        // the main shape
        let (num_constraints, num_io, num_vars) = (10, 3, 7);

        // get commitment_pp
        let ova_commitment_pp: Vec<Affine<G2>> = C2::setup(ova_shape.num_vars + ova_shape.num_constraints, b"test", &());
        let pp: Vec<Affine<G1>> = C1::setup(num_constraints + num_vars, b"test", &());


        let (shape, instance, witness) = get_random_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &pp);

        // assert it's satisfied
        shape.is_satisfied(&instance, &witness, &pp).expect("unsatisfied r1cs");

        // generate a relaxed instance/witness this time
        let (relaxed_shape, relaxed_instance, relaxed_witness) = get_random_relaxed_r1cs_instance_witness::<F, C1, G1>(num_constraints, num_vars, num_io, &pp);

        // assert the shape is equal to the previous shape
        assert_eq!(shape, relaxed_shape);

        // make sure the instance is satisfied
        shape.is_relaxed_satisfied(&relaxed_instance, &relaxed_witness, &pp).expect("unsatisfied r1cs");

        let p: NovaProver<F, G1, G2, C1, C2> = NovaProver {
            initial_transcript: Transcript::new(b"label"),
            final_transcript: Transcript::new(b"label"),
            shape: shape.clone(),
            commitment_pp: pp.clone(),
            current_accumulator: (instance, witness),
            running_accumulator: (relaxed_instance, relaxed_witness),
            ova_shape: ova_shape.clone(),
            ova_commitment_pp,
            cycle_fold_running_instance: RelaxedOvaInstance::new(&ova_shape),
            cycle_fold_running_witness: RelaxedOvaWitness::zero(&ova_shape),
        };

        let beta = F::rand(&mut thread_rng());
        let folded_accumulator = p.compute_final_accumulator(beta);

        shape.is_relaxed_satisfied(&folded_accumulator.0, &folded_accumulator.1, &pp).unwrap();
    }

    #[test]
    fn test_cycle_fold() {

    }
}

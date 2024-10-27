use crate::accumulation_circuit::instance_circuit::AccumulatorInstanceVar;
use crate::accumulation_circuit::verifier_circuit::{AccumulatorVerifier, AccumulatorVerifierVar};
use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::nexus_spartan::partial_verifier::partial_verifier::PartialVerifier;
use crate::nexus_spartan::partial_verifier::partial_verifier_var::PartialVerifierVar;
use crate::nova::cycle_fold::coprocessor_constraints::RelaxedOvaInstanceVar;
use crate::pcs::multilinear_pcs::split_between_x_and_y;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ec::{AffineRepr, CurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{Namespace, SynthesisError};
use itertools::izip;
use std::borrow::Borrow;

type Output<G2, C2, G1, F> = (RelaxedOvaInstanceVar<G2, C2>, AccumulatorInstanceVar<G1>, Vec<FpVar<F>>, Vec<FpVar<F>>);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AugmentedCircuit<G1, G2, C2, E, F>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    F: PrimeField,
{
    pub spartan_partial_verifier: PartialVerifier<F, E>,
    pub kzh_acc_verifier: AccumulatorVerifier<G1, G2, C2, E>,
}

pub struct AugmentedCircuitVar<G1, G2, C2, F>
where
    F: PrimeField + Absorb,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
{
    pub spartan_partial_verifier: PartialVerifierVar<F, G1>,
    pub kzh_acc_verifier: AccumulatorVerifierVar<G1, G2, C2>,
}

impl<G1, G2, C2, E, F> AllocVar<AugmentedCircuit<G1, G2, C2, E, F>, F> for AugmentedCircuitVar<G1, G2, C2, F>
where
    G1: SWCurveConfig + Clone,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField + Absorb,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField>,
    E: Pairing<G1Affine=Affine<G1>, ScalarField=F>,
    F: PrimeField,
{
    fn new_variable<T: Borrow<AugmentedCircuit<G1, G2, C2, E, F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        // Convert to Namespace<F>
        let ns = cs.into();
        // Get the constraint system reference
        let cs = ns.cs();

        // Fetch the instance of `AugmentedCircuit<F>`
        let binding = f()?;
        let data = binding.borrow();

        // Allocate the Spartan partial verifier
        let spartan_partial_verifier = PartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.spartan_partial_verifier),
            mode,
        )?;

        // Allocate the accumulator verifier
        let kzh_acc_verifier = AccumulatorVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.kzh_acc_verifier),
            mode,
        )?;

        Ok(AugmentedCircuitVar {
            spartan_partial_verifier,
            kzh_acc_verifier,
        })
    }
}

impl<G1, G2, C2, F> AugmentedCircuitVar<G1, G2, C2, F>
where
    F: PrimeField + Absorb,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
{
    fn verify<E: Pairing>(&self, transcript: &mut TranscriptVar<F>) -> Output<G2, C2, G1, F> {
        let (final_cycle_fold_instance, final_accumulator_instance) = self.kzh_acc_verifier.accumulate();
        let (rx, ry) = self.spartan_partial_verifier.verify(transcript);

        // ************* do the consistency checks *************
        let length_x = self.kzh_acc_verifier.current_accumulator_instance_var.x_var.len();
        let length_y = self.kzh_acc_verifier.current_accumulator_instance_var.y_var.len();

        let (expected_x_var, expected_y_var) = split_between_x_and_y(length_x, length_y, &ry[1..], FpVar::zero());
        for (e1, e2) in izip!(&self.kzh_acc_verifier.current_accumulator_instance_var.x_var, expected_x_var) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }

        for (e1, e2) in izip!(&self.kzh_acc_verifier.current_accumulator_instance_var.y_var, expected_y_var) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }


        // enforce equal eval_Z_at_ry and accumulator.z_var
        self.spartan_partial_verifier.eval_vars_at_ry.enforce_equal(
            &self.kzh_acc_verifier
                .current_accumulator_instance_var
                .z_var
        ).expect("error while enforcing equality");

        // enforce the commitment in spartan verifier and the accumulator new instance
        NonNativeAffineVar::enforce_equal(
            &self.spartan_partial_verifier.instance.1,
            &self.kzh_acc_verifier.current_accumulator_instance_var.C_var,
        ).expect("error while enforcing equality");

        (final_cycle_fold_instance, final_accumulator_instance, rx, ry)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::accumulation::accumulator::Accumulator;
    use crate::accumulation_circuit::prover::AccumulatorVerifierCircuitProver;
    use crate::constant_for_curves::{ScalarField, E, G1, G2};
    use crate::hash::pederson::PedersenCommitment;
    use crate::nexus_spartan::crr1cs::is_sat;
    use crate::nexus_spartan::crr1cs::produce_synthetic_crr1cs;
    use crate::nexus_spartan::crr1csproof::CRR1CSProof;
    use crate::nexus_spartan::polycommitments::{PolyCommitmentScheme, ToAffine};
    use crate::nova::cycle_fold::coprocessor::setup_shape;
    use crate::pcs::multilinear_pcs::PolyCommit;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use crate::transcript::transcript::Transcript;
    use ark_ff::AdditiveGroup;
    use ark_relations::r1cs::ConstraintSystem;
    use rand::thread_rng;

    type C2 = PedersenCommitment<Projective<G2>>;
    type F = ScalarField;

    #[test]
    fn test_augmented_circuit() {
        test_augmented_circuit_helper();
    }

    /// Take as input `proof_i` and `running_accumulator_{i}` and produce `proof_{i+1}` and `running_accumulator_{i+1}`.
    fn test_augmented_circuit_helper() {
        // ******************************* generate a satisfying instance for Spartan and get structure *******************************

        let num_vars = 1024;
        let num_cons = num_vars;
        let num_inputs = 10;

        // this generates a new instance/witness for spartan as well as PCS parameters
        let (shape, instance, witness, gens) = produce_synthetic_crr1cs::<E, MultilinearPolynomial<F>>(num_cons, num_vars, num_inputs);
        assert!(is_sat(&shape, &instance, &witness, &gens.gens_r1cs_sat).unwrap());

        let (num_cons, num_vars, _num_inputs) = (
            shape.get_num_cons(),
            shape.get_num_vars(),
            shape.get_num_inputs(),
        );

        let mut prover_transcript = Transcript::new(b"example");

        let (proof, rx, ry) = CRR1CSProof::prove(
            &shape,
            &instance,
            witness,
            &gens.gens_r1cs_sat,
            &mut prover_transcript,
        );

        let inst_evals = shape.inst.inst.evaluate(&rx, &ry);

        // ******************************* construct partial verifier circuit *******************************

        let mut prover_transcript = Transcript::new(b"example");
        let mut verifier_transcript = prover_transcript.clone();
        let verifier_transcript_clone = verifier_transcript.clone();
        let partial_verifier = PartialVerifier::initialise(
            &proof,
            num_vars,
            num_cons,
            (instance.input.assignment, {
                let com_w: <E as Pairing>::G1Affine = instance.comm_W.clone().to_affine();
                com_w
            }),
            &inst_evals,
            &mut prover_transcript,
        );

        partial_verifier.verify(&mut verifier_transcript);

        // ******************************* get corresponding accumulator of the spartan verifier *******************************

        // random running accumulator
        let pcs_srs = gens.gens_r1cs_sat.keys.ck.srs.clone();
        let acc_srs = Accumulator::setup(pcs_srs.clone(), &mut thread_rng());
        let running_acc = Accumulator::random_satisfying_accumulator(&acc_srs, &mut thread_rng());

        // generate an opening proof for the witness
        let opening_proof = proof.proof_eval_vars_at_ry.clone();

        // commit to witness again
        let commitment = instance.comm_W.clone();

        let length_x = pcs_srs.get_x_length();
        let length_y = pcs_srs.get_y_length();
        let (x, y) = split_between_x_and_y::<F>(length_x, length_y, &ry[1..], F::ZERO);

        assert!(
            PolyCommit::verify(
                &PolyCommit { srs: pcs_srs.clone() },
                &commitment,
                &opening_proof,
                x.as_slice(),
                y.as_slice(),
                &proof.eval_vars_at_ry
            )
        );

        let acc_instance = Accumulator::new_accumulator_instance_from_fresh_kzh_instance(
            &acc_srs,
            &commitment.C,
            x.as_slice(),
            y.as_slice(),
            &proof.eval_vars_at_ry,
        );

        let acc_witness = Accumulator::new_accumulator_witness_from_fresh_kzh_witness(
            &acc_srs,
            opening_proof,
            x.as_slice(),
            y.as_slice(),
        );

        let current_acc = Accumulator::new_accumulator(&acc_instance, &acc_witness);

        // assert decide
        assert!(
            Accumulator::decide(
                &acc_srs,
                &current_acc,
            )
        );

        // ******************************* Construct the accumulator verifier circuit *******************************

        // the shape of the R1CS instance
        let shape = setup_shape::<G1, G2>().unwrap();

        // get trivial running instance
        let (cycle_fold_running_instance, cycle_fold_running_witness) = AccumulatorVerifierCircuitProver::<G1, G2, C2, E>::get_trivial_cycle_fold_running_instance_witness(&shape);

        // get commitment_pp
        let commitment_pp = AccumulatorVerifierCircuitProver::<G1, G2, C2, E>::get_commitment_pp(&shape);


        let prover: AccumulatorVerifierCircuitProver<G1, G2, C2, E> = AccumulatorVerifierCircuitProver::new(
            &acc_srs,
            commitment_pp,
            running_acc,
            current_acc,
            cycle_fold_running_instance,
            cycle_fold_running_witness,
        );

        // assert it's formated correctly
        prover.is_satisfied();

        // ******************************* Construct the augmented circuit *******************************

        let cs = ConstraintSystem::<ScalarField>::new_ref();
        let partial_verifier_var = PartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(partial_verifier.clone()),
            AllocationMode::Input,
        ).unwrap();
        let acc_verifier_var = AccumulatorVerifierVar::<G1, G2, C2>::new::<E>(cs.clone(), prover);

        let augmented_circuit = AugmentedCircuitVar {
            spartan_partial_verifier: partial_verifier_var,
            kzh_acc_verifier: acc_verifier_var,
        };
        let mut transcript_var = TranscriptVar::from_transcript(cs.clone(), verifier_transcript_clone);

        augmented_circuit.verify::<E>(&mut transcript_var);

        assert!(cs.is_satisfied().unwrap());
        println!("{}", cs.num_constraints());
    }
}
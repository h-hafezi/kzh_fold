use crate::commitment::CommitmentScheme;
use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::hash::poseidon::PoseidonHashVar;
use crate::kzh::kzh3::{KZH3, KZH3SRS};
use crate::kzh::KZH;
use crate::kzh3_verifier_circuit::verifier_circuit::{KZH3Verifier, KZH3VerifierVar};
use crate::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationAccVerifier, MatrixEvaluationAccVerifierVar};
use crate::nexus_spartan::partial_verifier::partial_verifier::SpartanPartialVerifier;
use crate::nexus_spartan::partial_verifier::partial_verifier_var::SpartanPartialVerifierVar;
use crate::transcript::transcript_var::TranscriptVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use itertools::izip;
use rand::thread_rng;
use std::borrow::Borrow;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KZH3AugmentedCircuit<G1, G2, C2, E, F>
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
    pub spartan_partial_verifier: SpartanPartialVerifier<F, E>,
    pub kzh_acc_verifier: KZH3Verifier<G1, G2, C2, E>,
    pub matrix_evaluation_verifier: MatrixEvaluationAccVerifier<F>,
}

pub struct KZH3AugmentedCircuitVar<G1, G2, C2, F>
where
    F: PrimeField + Absorb,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
{
    pub spartan_partial_verifier: SpartanPartialVerifierVar<F, G1>,
    pub kzh_acc_verifier: KZH3VerifierVar<G1, G2, C2>,
    pub matrix_evaluation_verifier: MatrixEvaluationAccVerifierVar<F>,
}

impl<G1, G2, C2, E, F> AllocVar<KZH3AugmentedCircuit<G1, G2, C2, E, F>, F> for KZH3AugmentedCircuitVar<G1, G2, C2, F>
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
    fn new_variable<T: Borrow<KZH3AugmentedCircuit<G1, G2, C2, E, F>>>(
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
        let spartan_partial_verifier = SpartanPartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.spartan_partial_verifier),
            mode,
        )?;

        // Allocate the accumulator verifier
        let kzh_acc_verifier = KZH3VerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.kzh_acc_verifier),
            mode,
        )?;

        // Allocate the accumulator verifier
        let matrix_evaluation_verifier = MatrixEvaluationAccVerifierVar::new_variable(
            cs.clone(),
            || Ok(&data.matrix_evaluation_verifier),
            mode,
        )?;

        Ok(KZH3AugmentedCircuitVar {
            spartan_partial_verifier,
            kzh_acc_verifier,
            matrix_evaluation_verifier,
        })
    }
}

impl<G1, G2, C2, F> KZH3AugmentedCircuitVar<G1, G2, C2, F>
where
    F: PrimeField + Absorb,
    G1::BaseField: PrimeField,
    G1::ScalarField: PrimeField,
    G2: SWCurveConfig<BaseField=F> + Clone,
    G2::BaseField: PrimeField,
    C2: CommitmentScheme<Projective<G2>>,
    G1: SWCurveConfig<BaseField=G2::ScalarField, ScalarField=G2::BaseField> + Clone,
{
    pub fn verify<E: Pairing>(&self, pcs_srs: &KZH3SRS<E>, cs: ConstraintSystemRef<F>, transcript: &mut TranscriptVar<F>, poseidon_num: usize)
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as ark_ec::AffineRepr>::BaseField: PrimeField,
    {
        let (_rx, ry) = self.spartan_partial_verifier.verify(transcript);
        let _ = self.kzh_acc_verifier.accumulate(transcript);

        // also return these later
        let _ = self.matrix_evaluation_verifier.accumulate(transcript);

        // ************* do the consistency checks *************
        let split_input = KZH3::split_input(&pcs_srs, &ry[1..], FpVar::zero());
        for (e1, e2) in izip!(&self.kzh_acc_verifier.current_accumulator_instance_var.x_var, split_input[0].clone()) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }

        for (e1, e2) in izip!(&self.kzh_acc_verifier.current_accumulator_instance_var.y_var, split_input[1].clone()) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }

        for (e1, e2) in izip!(&self.kzh_acc_verifier.current_accumulator_instance_var.z_var, split_input[2].clone()) {
            e1.enforce_equal(&e2).expect("error while enforcing equality");
        }


        // enforce equal eval_Z_at_ry and accumulator.z_var
        self.spartan_partial_verifier.eval_vars_at_ry.enforce_equal(
            &self.kzh_acc_verifier
                .current_accumulator_instance_var
                .output
        ).expect("error while enforcing equality");

        // enforce the commitment in spartan verifier and the accumulator new instance
        NonNativeAffineVar::enforce_equal(
            &self.spartan_partial_verifier.instance.1,
            &self.kzh_acc_verifier.current_accumulator_instance_var.C_var,
        ).expect("error while enforcing equality");

        // pad it with some random poseidon hash
        let mut hash = PoseidonHashVar::new(cs.clone());
        for _ in 0..poseidon_num {
            // get a random element
            let r = FpVar::new_variable(cs.clone(), || Ok(F::rand(&mut thread_rng())), AllocationMode::Witness).unwrap();
            // update sponge with this random element
            hash.update_sponge(vec![r]);
            // output the hash
            let _ = hash.output();
        }
    }
}

#[cfg(test)]
mod test {
    use ark_std::{end_timer, start_timer};
    use crate::constant_for_curves::{ScalarField as F, C2, E, G1, G2};
    use crate::kzh::kzh3::{KZH3, KZH3SRS};
    use crate::kzh::KZH;
    use crate::kzh_fold::kzh3_fold::Accumulator3 as Accumulator;
    use crate::nexus_spartan::commitment_traits::ToAffine;
    use crate::nexus_spartan::committed_relaxed_snark::CRSNARKKey;
    use crate::nexus_spartan::crr1cs::{is_sat, produce_synthetic_crr1cs, CRR1CSInstance, CRR1CSShape, CRR1CSWitness};
    use crate::nexus_spartan::crr1csproof::CRR1CSProof;
    use crate::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationAccVerifier, MatrixEvaluationAccVerifierVar};
    use crate::nexus_spartan::partial_verifier::partial_verifier::SpartanPartialVerifier;
    use crate::nexus_spartan::partial_verifier::partial_verifier_var::SpartanPartialVerifierVar;
    use crate::nova::cycle_fold::coprocessor::setup_shape;
    use crate::transcript::transcript::Transcript;
    use crate::transcript::transcript_var::TranscriptVar;
    use ark_ec::pairing::Pairing;
    use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
    use ark_relations::r1cs::{ConstraintSystem, SynthesisMode};
    use ark_serialize::CanonicalSerialize;
    use rand::thread_rng;
    use crate::kzh3_augmented_circuit::kzh3_augmented_circuit::KZH3AugmentedCircuitVar;
    use crate::kzh3_verifier_circuit::prover::KZH3VerifierCircuitProver;
    use crate::kzh3_verifier_circuit::verifier_circuit::KZH3VerifierVar;
    use crate::kzh::kzh2::KZH2;
    use crate::nexus_spartan::conversion::convert_crr1cs;

    #[test]
    fn test_kzh3_augmented_circuit() {
        let poseidon_num = 0;

        // Directly map poseidon_num to (num_vars, num_inputs)
        let (num_vars, num_inputs) = match poseidon_num {
            0 => (131072, 10),
            150 => (524288, 1685),
            1000 => (1048576, 1529),
            2000 => (2097152, 4623),
            _ => {
                eprintln!("Error: Invalid poseidon_num value!");
                return;
            }
        };

        println!("[*] KZH3: Poseidon num: {}, Num cons: {}", poseidon_num, num_vars);

        let (pcs_srs, spartan_shape, spartan_instance, spartan_proof, rx, ry) = {
            let num_cons = num_vars;


            // this generates a new instance/witness for spartan as well as PCS parameters
            let (spartan_shape, spartan_instance, spartan_witness, spartan_key) = produce_synthetic_crr1cs::<E, KZH3<E>>(num_cons, num_vars, num_inputs);

            assert!(is_sat(&spartan_shape, &spartan_instance, &spartan_witness, &spartan_key.gens_r1cs_sat).unwrap());

            let pcs_srs = spartan_key.gens_r1cs_sat.clone();

            let mut prover_transcript = Transcript::new(b"example");

            // Get `proof_i` and random evaluation point (r_x, r_y)
            let (spartan_proof, rx, ry) = CRR1CSProof::prove(
                &spartan_shape,
                &spartan_instance,
                spartan_witness,
                &spartan_key.gens_r1cs_sat,
                &mut prover_transcript,
            );

            (pcs_srs, spartan_shape, spartan_instance, spartan_proof, rx, ry)
        };

        // fresh transcripts to be used by the prover and verifier
        let mut prover_transcript = Transcript::new(b"example");
        let verifier_transcript_clone = prover_transcript.clone();
        let cs = ConstraintSystem::<F>::new_ref();

        let partial_verifier_var = {
            let mut verifier_transcript = prover_transcript.clone();
            // Get A(r_x, r_y), B(r_x, r_y), C(r_x, r_y)
            let current_A_B_C_evaluations = spartan_shape.inst.inst.evaluate(&rx, &ry);

            let partial_verifier = SpartanPartialVerifier::initialise(
                &spartan_proof,
                spartan_shape.get_num_vars(),
                spartan_shape.get_num_cons(),
                (spartan_instance.input.assignment, {
                    let com_w: <E as Pairing>::G1Affine = spartan_instance.comm_W.clone().to_affine();
                    com_w
                }),
                &current_A_B_C_evaluations,
                &mut prover_transcript,
            );

            partial_verifier.verify(&mut verifier_transcript);

            let partial_verifier_var = SpartanPartialVerifierVar::new_variable(
                cs.clone(),
                || Ok(partial_verifier.clone()),
                AllocationMode::Input,
            ).unwrap();

            partial_verifier_var
        };

        let acc_verifier_var = {
            let acc_srs = Accumulator::setup(pcs_srs.clone(), &mut thread_rng());

            // Get the KZH opening proof from the Spartan proof
            let opening_proof = spartan_proof.proof_eval_vars_at_ry.clone();

            // Commitment to witness polynomial
            let commitment_w = spartan_instance.comm_W.clone();

            let input = &ry[1..];

            // Sanity check: verify the opening proof
            KZH3::verify(
                &pcs_srs,
                &input,
                &spartan_proof.eval_vars_at_ry,
                &commitment_w,
                &opening_proof,
            );

            // Get accumulator from the opening proof
            let acc_instance = Accumulator::proof_to_accumulator_instance(
                &acc_srs,
                &input,
                &spartan_proof.eval_vars_at_ry,
                &commitment_w,
                &opening_proof,
            );

            let acc_witness = Accumulator::proof_to_accumulator_witness(
                &acc_srs,
                commitment_w,
                opening_proof,
                &input,
            );

            let current_acc = Accumulator::new(&acc_instance, &acc_witness);

            // println!("proof size: {}", proof.compressed_size());
            println!("acc size: {}", current_acc.compressed_size());

            // Check that the accumulator is valid
            Accumulator::decide(
                &acc_srs,
                &current_acc,
            );

            // use a random accumulator as the running one
            let running_acc = Accumulator::rand(&acc_srs);

            // the shape of the R1CS instance
            let ova_shape = setup_shape::<G1, G2>().unwrap();

            // get trivial running instance
            let (ova_running_instance, ova_running_witness) = KZH3VerifierCircuitProver::<G1, G2, C2, E, F>::get_trivial_cycle_fold_running_instance_witness(&ova_shape);

            // get commitment_pp
            let ova_commitment_pp = KZH3VerifierCircuitProver::<G1, G2, C2, E, F>::get_commitment_pp(&ova_shape);

            let kzh_acc_verifier_prover: KZH3VerifierCircuitProver<G1, G2, C2, E, F> = KZH3VerifierCircuitProver::new(
                &acc_srs,
                ova_commitment_pp,
                running_acc,
                current_acc.clone(),
                ova_running_instance,
                ova_running_witness,
                prover_transcript,
            );

            // assert it's formated correctly
            kzh_acc_verifier_prover.is_satisfied();

            let acc_verifier_var = KZH3VerifierVar::<G1, G2, C2>::new::<E>(cs.clone(), kzh_acc_verifier_prover);

            acc_verifier_var
        };

        let matrix_evaluation_verifier_var = {
            let matrix_eval_acc_verifier = MatrixEvaluationAccVerifier::random_from_eval_point(
                &spartan_shape,
                rx.clone(),
                ry.clone(),
                &mut thread_rng(),
            );

            let matrix_evaluation_verifier_var = MatrixEvaluationAccVerifierVar::new_variable(
                cs.clone(),
                || Ok(matrix_eval_acc_verifier.clone()),
                AllocationMode::Input,
            ).unwrap();

            matrix_evaluation_verifier_var
        };

        // construct the augmented circuit
        let augmented_circuit = KZH3AugmentedCircuitVar {
            spartan_partial_verifier: partial_verifier_var,
            kzh_acc_verifier: acc_verifier_var,
            matrix_evaluation_verifier: matrix_evaluation_verifier_var,
        };

        let mut transcript_var = TranscriptVar::from_transcript(cs.clone(), verifier_transcript_clone);

        // run the verification function on augmented circuit
        let _ = augmented_circuit.verify::<E>(&pcs_srs, cs.clone(), &mut transcript_var, poseidon_num);

        assert!(cs.is_satisfied().unwrap());
        println!("augmented circuit constraints: {}", cs.num_constraints());

        // Set the mode to Prove before we convert it for spartan
        cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
        cs.finalize();

        ////////// Now run the spartan prover on the augmented circuit /////////////////

        // convert to the corresponding Spartan types
        let shape = CRR1CSShape::<F>::convert::<G1>(cs.clone());

        // get the number the minimum size we need for committing to the constraint system
        let min_num_vars = CRSNARKKey::<E, KZH3<E>>::get_min_num_vars(shape.get_num_cons(), shape.get_num_vars(), shape.get_num_inputs());
        let SRS: KZH3SRS<E> = KZH3::setup(min_num_vars + 1, &mut thread_rng());

        let (instance, witness): (CRR1CSInstance<E, KZH3<E>>, CRR1CSWitness<E, KZH3<E>>) = convert_crr1cs(cs.clone(), &SRS);

        let mut new_prover_transcript = Transcript::new(b"example");
        let (proof, rx, ry) = CRR1CSProof::prove(
            &shape,
            &instance,
            witness,
            &SRS,
            &mut new_prover_transcript,
        );

        // evaluate matrices A B C
        let A_B_C_eval_timer = start_timer!(|| "ABC evals");
        let inst_evals = shape.inst.inst.evaluate(&rx, &ry);
        end_timer!(A_B_C_eval_timer);

        let mut new_verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(
                shape.get_num_vars(),
                shape.get_num_cons(),
                &instance,
                &inst_evals,
                &mut new_verifier_transcript,
            )
            .is_ok());
    }
}

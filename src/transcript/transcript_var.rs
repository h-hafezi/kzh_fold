use crate::gadgets::non_native::util::non_native_to_fpvar;
use crate::hash::poseidon::PoseidonHashVar;
use crate::transcript::transcript::Transcript;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_relations::r1cs::ConstraintSystemRef;

pub struct TranscriptVar<F: PrimeField + Absorb> {
    // This will hold the current state of the transcript
    pub state: FpVar<F>,
    // the poseidon hash
    pub poseidon_hash: PoseidonHashVar<F>,
}

impl<F: Absorb + PrimeField> TranscriptVar<F> {
    pub fn new(cs: ConstraintSystemRef<F>, label: &'static [u8]) -> Self {
        let trans: Transcript<F> = Transcript::new(label);

        TranscriptVar {
            state: FpVar::new_input(cs.clone(), || Ok(trans.state)).unwrap(),
            poseidon_hash: PoseidonHashVar::new(cs.clone()),
        }
    }

    pub fn from_transcript(cs: ConstraintSystemRef<F>, transcript: Transcript<F>) -> TranscriptVar<F> {
        let state = FpVar::new_witness(
            cs.clone(),
            || Ok(transcript.state.clone()),
        ).unwrap();

        let poseidon_hash = PoseidonHashVar::from_poseidon_hash(
            cs.clone(),
            transcript.poseidon_hash,
        );

        TranscriptVar {
            state,
            poseidon_hash,
        }
    }
}


impl<F: PrimeField + Absorb> TranscriptVar<F> {
    pub fn append_message(&mut self, _label: &'static [u8], _message: &[u8]) {
        // do not do anything
    }

    pub fn append_scalar(&mut self, _label: &'static [u8], scalar: &FpVar<F>) {
        self.poseidon_hash.update_sponge(vec![scalar.clone()]);
    }

    pub fn append_scalar_non_native<Q: PrimeField>(&mut self, _label: &'static [u8], scalar: &NonNativeFieldVar<Q, F>) {
        let converted_scalar = non_native_to_fpvar(&scalar);
        self.poseidon_hash.update_sponge(vec![converted_scalar.clone()]);
    }

    pub fn append_scalars(&mut self, _label: &'static [u8], scalars: &[FpVar<F>]) {
        for f in scalars {
            self.append_scalar(_label, f);
        }
    }

    pub fn append_scalars_non_native<Q: PrimeField>(&mut self, _label: &'static [u8], scalars: &[NonNativeFieldVar<Q, F>]) {
        for q in scalars {
            self.append_scalar_non_native(_label, q);
        }
    }

    pub fn challenge_scalar(&mut self, _label: &'static [u8]) -> FpVar<F> {
        let new_state = self.poseidon_hash.output();
        self.state = new_state.clone();
        self.append_scalar(_label, &new_state);
        new_state
    }

    pub fn challenge_vector(&mut self, _label: &'static [u8], len: usize) -> Vec<FpVar<F>> {
        let mut res = Vec::with_capacity(len);
        for _ in 0..len {
            res.push(self.challenge_scalar(_label));
        }
        res
    }
}

pub trait AppendToTranscriptVar<F: PrimeField + Absorb> {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut TranscriptVar<F>);
}


#[cfg(test)]
mod tests {
    use crate::constant_for_curves::{BaseField, ScalarField};
    use crate::transcript::transcript::Transcript;
    use crate::transcript::transcript_var::TranscriptVar;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
    use ark_std::UniformRand;
    use rand::thread_rng;

    type F = ScalarField;
    type Q = BaseField;

    #[test]
    fn test_transcript_vs_transcript_var() {
        // Initialize constraint system
        let cs: ConstraintSystemRef<F> = ConstraintSystem::new_ref();

        // Create random scalar in F
        let mut rng = thread_rng();

        // Initialize the transcript
        let label = b"test_label";
        let mut transcript = Transcript::new(label);

        // Initialize the transcript_var
        let mut transcript_var = TranscriptVar::new(cs.clone(), label);

        for _ in 0..10 {
            let random_scalar: F = F::rand(&mut rng);
            let non_native: Q = Q::rand(&mut rng);
            // Append the same random scalar to both transcript and transcript_var
            transcript.append_scalar(label, &random_scalar);
            transcript.append_scalar_non_native(label, &non_native);

            // the zk version
            let random_scalar_var = FpVar::new_witness(
                cs.clone(),
                || Ok(random_scalar),
            ).unwrap();
            let random_non_native = NonNativeFieldVar::new_witness(
                cs.clone(),
                || Ok(non_native),
            ).unwrap();

            // append to transcript
            transcript_var.append_scalar(label, &random_scalar_var);
            transcript_var.append_scalar_non_native(label, &random_non_native);
        }

        // Compare the final states
        let transcript_state = transcript.challenge_scalar(label);
        let transcript_var_state = transcript_var.challenge_scalar(label).value().unwrap();

        // Assert that the states are equal
        assert_eq!(transcript_state, transcript_var_state, "The transcript states are not equal");
        println!("number of constraints: {}", cs.num_constraints());

        // test new_variable works correctly
        let mut new_transcript_var = TranscriptVar::from_transcript(
            cs.clone(),
            transcript,
        );

        // assert output of new_transcript_var and transcript_var are identical
        assert_eq!(
            transcript_var.challenge_scalar(label).value().unwrap(),
            new_transcript_var.challenge_scalar(label).value().unwrap()
        );
    }
}
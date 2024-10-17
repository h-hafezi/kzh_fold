use ark_crypto_primitives::crh::poseidon::constraints::CRHParametersVar;
use ark_crypto_primitives::sponge::Absorb;
use ark_crypto_primitives::sponge::poseidon::constraints::PoseidonSpongeVar;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_relations::r1cs::ConstraintSystemRef;
use crate::hash::poseidon::{get_poseidon_config, PoseidonHash, PoseidonHashVar};
use crate::transcript::transcript::Transcript;

pub struct TranscriptVar<F: PrimeField + Absorb> {
    // This will hold the current state of the transcript
    state: FpVar<F>,
    // the poseidon hash
    poseidon_hash: PoseidonHashVar<F>,
}

impl<F: Absorb + PrimeField> TranscriptVar<F> {
    pub fn new(cs: ConstraintSystemRef<F>, label: &[u8]) -> Self {
        let trans = Transcript::new(label);
        
        TranscriptVar {
            state: FpVar::new_input(cs.clone(), || Ok(trans.state)).unwrap(),
            poseidon_hash: PoseidonHashVar::new(cs.clone()),
        }
    }
}

impl<F: PrimeField + Absorb> TranscriptVar<F> {
    pub fn append_message(&mut self, _label: &'static [u8], _msg: &[u8]) {
        // I'm not sure if it's important to implement this
    }

    pub fn append_scalar(&mut self, _label: &'static [u8], scalar: &FpVar<F>) {
        self.poseidon_hash.update_sponge(vec![scalar.clone()]);
    }

    pub fn append_scalars(&mut self, _label: &'static [u8], scalars: &[FpVar<F>]) {
        for f in scalars {
            self.append_scalar(_label, f);
        }
    }

    /*pub fn challenge_scalar(&mut self, _label: &'static [u8]) -> F {
        //let new_state = self.poseidon_hash.output();
        //self.state = new_state;
        new_state
    }

    pub fn challenge_vector(&mut self, _label: &'static [u8], len: usize) -> Vec<F> {
        let mut res = Vec::with_capacity(len);
        for _ in 0..len {
            res.push(self.challenge_scalar(_label));
        }
        res
    }
    
     */
}

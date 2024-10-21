use crate::nexus_spartan::unipoly::unipoly::{CompressedUniPoly, UniPoly};
use crate::transcript::transcript_var::{AppendToTranscriptVar, TranscriptVar};
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;

#[derive(Clone)]
pub struct UniPolyVar<F: PrimeField + Absorb> {
    coeffs: Vec<FpVar<F>>,
}

// ax^2 + bx + c stored as vec![c,a]
// ax^3 + bx^2 + cx + d stored as vec![d,b,a]
#[derive(Clone)]
pub struct CompressedUniPolyVar<F: PrimeField + Absorb> {
    coeffs_except_linear_term: Vec<FpVar<F>>,
}

impl<F: PrimeField + Absorb> AllocVar<UniPoly<F>, F> for UniPolyVar<F> {
    fn new_variable<T: Borrow<UniPoly<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        // Fetch the vector of coefficients to allocate as FpVars
        let binding = f()?;
        let values = binding.borrow();

        // Allocate each coefficient as an FpVar
        let coeffs = values
            .coeffs
            .iter()
            .map(|v| FpVar::new_variable(cs.clone(), || Ok(v), mode))
            .collect::<Result<Vec<_>, SynthesisError>>()?;

        Ok(UniPolyVar { coeffs })
    }
}

impl<F: PrimeField + Absorb> R1CSVar<F> for UniPolyVar<F> {
    type Value = UniPoly<F>;

    // This method returns the constraint system that the variables are attached to
    fn cs(&self) -> ConstraintSystemRef<F> {
        let mut result = ConstraintSystemRef::None;
        for coeff in &self.coeffs {
            result = coeff.cs().or(result);
        }
        result
    }

    // This method returns the underlying values of the variables, if available
    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let mut coeffs = Vec::new();
        for c in &self.coeffs {
            coeffs.push(c.value()?);
        }
        Ok(UniPoly { coeffs })
    }
}

impl<F: PrimeField + Absorb> UniPolyVar<F> {
    pub fn eval_at_zero(&self) -> FpVar<F> {
        self.coeffs[0].clone()
    }

    pub fn eval_at_one(&self) -> FpVar<F> {
        let mut result = FpVar::<F>::zero();
        for i in 0..self.coeffs.len() {
            result += self.coeffs[i].clone();
        }
        result
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }
    pub fn evaluate(&self, r: &FpVar<F>) -> FpVar<F> {
        let mut eval = self.coeffs[0].clone();
        let mut power = r.clone();
        for i in 1..self.coeffs.len() {
            eval += power.clone() * self.coeffs[i].clone();
            power *= r;
        }
        eval
    }
}

impl<F: PrimeField + Absorb> CompressedUniPolyVar<F> {
    pub fn new(cs: ConstraintSystemRef<F>, poly: CompressedUniPoly<F>, allocation_mode: AllocationMode) -> Self {
        CompressedUniPolyVar {
            coeffs_except_linear_term: {
                let mut coeffs = Vec::new();
                for fp in poly.coeffs_except_linear_term {
                    let fpvar = FpVar::new_variable(cs.clone(), || Ok(fp.clone()), allocation_mode).unwrap();
                    coeffs.push(fpvar);
                }
                coeffs
            },
        }
    }
}

impl<F: PrimeField + Absorb> CompressedUniPolyVar<F> {
    // we require eval(0) + eval(1) = hint, so we can solve for the linear term as:
    // linear_term = hint - 2 * constant_term - deg2 term - deg3 term
    pub fn decompress(&self, hint: &FpVar<F>) -> UniPolyVar<F> {
        let mut linear_term = hint - self.coeffs_except_linear_term[0].clone() - self.coeffs_except_linear_term[0].clone();
        for i in 1..self.coeffs_except_linear_term.len() {
            linear_term -= self.coeffs_except_linear_term[i].clone();
        }

        let mut coeffs = vec![self.coeffs_except_linear_term[0].clone(), linear_term];
        coeffs.extend({
            let res = self.coeffs_except_linear_term[1..].to_vec();
            res
        });
        assert_eq!(self.coeffs_except_linear_term.len() + 1, coeffs.len());
        UniPolyVar { coeffs }
    }
}


// This has to be consistent with the append_to_transcript for UniPoly
impl<F: PrimeField + Absorb> AppendToTranscriptVar<F> for UniPolyVar<F> {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut TranscriptVar<F>) {
        for i in 0..self.coeffs.len() {
            TranscriptVar::append_scalar(transcript, label, &self.coeffs[i]);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::constant_for_curves::ScalarField;
    use crate::nexus_spartan::unipoly::unipoly::{CompressedUniPoly, UniPoly};
    use crate::nexus_spartan::unipoly::unipoly_var::{CompressedUniPolyVar, UniPolyVar};
    use crate::transcript::transcript::{AppendToTranscript, Transcript};
    use crate::transcript::transcript_var::{AppendToTranscriptVar, TranscriptVar};
    use ark_ff::UniformRand;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::alloc::AllocationMode;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use rand::thread_rng;

    type F = ScalarField;

    #[test]
    fn test_transcript_unipoly_vs_unipolyvar() {
        // Initialize a random number generator
        let mut rng = thread_rng();

        // Generate a random UniPoly
        let random_unipoly = UniPoly {
            coeffs: vec![
                F::rand(&mut rng),
                F::rand(&mut rng),
                F::rand(&mut rng),
            ],
        };

        let label = b"label";
        // Create a transcript and append the random UniPoly
        let mut transcript = Transcript::new(label);
        random_unipoly.append_to_transcript(label, &mut transcript);

        // Create a constraint system
        let cs = ConstraintSystem::<F>::new_ref();

        // Generate the UniPolyVar version from the UniPoly
        let random_unipoly_var = UniPolyVar::new_variable(
            cs.clone(),
            || Ok(random_unipoly.clone()),
            AllocationMode::Witness,
        ).unwrap();

        // test value function
        assert_eq!(random_unipoly, random_unipoly_var.value().unwrap());

        // Create the transcript var and append the UniPolyVar
        let mut transcript_var = TranscriptVar::new(cs.clone(), label);
        random_unipoly_var.append_to_transcript(label, &mut transcript_var);

        // Compare the state of both transcripts
        assert_eq!(transcript.challenge_scalar(label), transcript_var.challenge_scalar(label).value().unwrap());
    }

    #[test]
    fn test_decompress_unipoly_vs_unipolyvar() {
        // Initialize a random number generator
        let mut rng = thread_rng();

        // Generate a random CompressedUniPoly
        let random_compressed_unipoly = CompressedUniPoly {
            coeffs_except_linear_term: vec![
                F::rand(&mut rng),
                F::rand(&mut rng),
                F::rand(&mut rng),
            ],
        };

        // Generate a hint
        let hint = F::rand(&mut rng);

        // Decompress the CompressedUniPoly
        let decompressed_unipoly = random_compressed_unipoly.decompress(&hint);

        // Create a constraint system
        let cs = ConstraintSystem::<F>::new_ref();

        // Create the CompressedUniPolyVar version from the CompressedUniPoly
        let random_compressed_unipoly_var = CompressedUniPolyVar::new(cs.clone(), random_compressed_unipoly.clone(), AllocationMode::Witness);

        // Generate a variable for the hint
        let hint_var = FpVar::new_input(cs.clone(), || Ok(hint)).unwrap();

        // Decompress the CompressedUniPolyVar
        let decompressed_unipoly_var = random_compressed_unipoly_var.decompress(&hint_var);

        // Compare the coefficients of both decompressed versions
        for (coeff, coeff_var) in decompressed_unipoly.coeffs.iter().zip(decompressed_unipoly_var.coeffs.iter()) {
            assert_eq!(coeff, &coeff_var.value().unwrap());
        }
    }
}
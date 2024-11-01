use ark_crypto_primitives::sponge::Absorb;
use ark_ec::AffineRepr;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_poly_commit::Error;
use rand::RngCore;

use crate::nexus_spartan::polycommitments::error::PCSError;
use crate::nexus_spartan::polycommitments::{PCSKeys, PolyCommitmentScheme, ToAffine};
use crate::pcs::multilinear_pcs::{split_between_x_and_y, PCSCommitment, PCSOpeningProof, PCSEngine, PolynomialCommitmentSRS};
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for PCSCommitment<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>) {
        Transcript::append_point::<E>(transcript, label, &self.C);
    }
}

impl<E: Pairing> ToAffine<E> for PCSCommitment<E> {
    fn to_affine(self) -> E::G1Affine {
        self.C
    }
}

impl<F: PrimeField + Absorb, E: Pairing<ScalarField=F>> PolyCommitmentScheme<E> for MultilinearPolynomial<F>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    type SRS = PolynomialCommitmentSRS<E>;
    type PolyCommitmentKey = PCSEngine<E>;
    type EvalVerifierKey = PCSEngine<E>;
    type Commitment = PCSCommitment<E>;
    type PolyCommitmentProof = PCSOpeningProof<E>;

    fn commit(poly: &MultilinearPolynomial<E::ScalarField>, ck: &Self::PolyCommitmentKey) -> Self::Commitment {
        let len = ck.srs.get_x_length() + ck.srs.get_y_length();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        ck.commit(&poly)
    }

    fn prove(C: Option<&Self::Commitment>, poly: &MultilinearPolynomial<E::ScalarField>, r: &[E::ScalarField], ck: &Self::PolyCommitmentKey) -> Self::PolyCommitmentProof {
        let len = ck.srs.get_x_length() + ck.srs.get_y_length();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        let length_x = ck.srs.get_x_length();
        let length_y = ck.srs.get_y_length();
        let (x, _) = split_between_x_and_y::<F>(length_x, length_y, r, F::ZERO);
        ck.open(&poly, C.unwrap().clone(), x.as_slice())
    }

    fn verify(commitment: &Self::Commitment, proof: &Self::PolyCommitmentProof, ck: &Self::EvalVerifierKey, r: &[F], eval: &F) -> Result<(), PCSError> {
        let length_x = ck.srs.get_x_length();
        let length_y = ck.srs.get_y_length();
        let (x, y) = split_between_x_and_y::<F>(length_x, length_y, r, F::ZERO);

        // verify the proof
        assert!(ck.verify(commitment, proof, x.as_slice(), y.as_slice(), eval));

        Ok(())
    }

    fn setup(max_poly_vars: usize, rng: &mut impl RngCore) -> Result<Self::SRS, Error> {
        let x = max_poly_vars / 2;
        let y = max_poly_vars - x;
        let degree_x = 2usize.pow(x as u32);
        let degree_y = 2usize.pow(y as u32);
        Ok(PCSEngine::<E>::setup(degree_x, degree_y, rng))
    }

    fn trim(srs: &Self::SRS) -> PCSKeys<E, Self> {
        PCSKeys {
            ck: PCSEngine { srs: srs.clone() },
            vk: PCSEngine { srs: srs.clone() },
        }
    }
}

#[cfg(test)]
pub mod test {
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{ScalarField, E};
    use crate::nexus_spartan::polycommitments::{PCSKeys, PolyCommitmentScheme};
    use crate::pcs::multilinear_pcs::PolynomialCommitmentSRS;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

    #[test]
    fn test_end_to_end() {
        let srs: PolynomialCommitmentSRS<E> = MultilinearPolynomial::<ScalarField>::setup(5, &mut thread_rng()).unwrap();
        let PCSKeys { vk: _, ck } = MultilinearPolynomial::trim(&srs);

        // random bivariate polynomial
        let polynomial = MultilinearPolynomial::<ScalarField>::rand(3, &mut thread_rng());

        // random points and evaluation
        let x = vec![
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
        ];

        let z = polynomial.evaluate(&x);

        let com = MultilinearPolynomial::commit(&polynomial, &ck);
        let open = MultilinearPolynomial::prove(Option::from(&com), &polynomial, x.as_slice(), &ck);
        MultilinearPolynomial::verify(&com, &open, &ck, x.as_slice(), &z).expect("verification failed");
    }
}

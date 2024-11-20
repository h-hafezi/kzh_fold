use ark_crypto_primitives::sponge::Absorb;
use ark_ec::AffineRepr;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_poly_commit::Error;
use rand::RngCore;
use crate::math::Math;
use crate::nexus_spartan::polycommitments::error::PCSError;
use crate::nexus_spartan::polycommitments::{PolyCommitmentScheme, ToAffine};
use crate::kzh::kzh2::{split_between_x_and_y, KZH2Commitment, KZH2OpeningProof, KZH2, KZH2SRS};
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for KZH2Commitment<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>) {
        Transcript::append_point::<E>(transcript, label, &self.C);
    }
}

impl<E: Pairing> ToAffine<E> for KZH2Commitment<E> {
    fn to_affine(self) -> E::G1Affine {
        self.C
    }
}

impl<F: PrimeField + Absorb, E: Pairing<ScalarField=F>> PolyCommitmentScheme<E> for MultilinearPolynomial<F>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    type SRS = KZH2SRS<E>;
    type Commitment = KZH2Commitment<E>;
    type PolyCommitmentProof = KZH2OpeningProof<E>;

    fn commit(poly: &MultilinearPolynomial<E::ScalarField>, srs: &Self::SRS) -> Self::Commitment {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        KZH2::commit_1(&srs, &poly)
    }

    fn prove(C: Option<&Self::Commitment>, poly: &MultilinearPolynomial<E::ScalarField>, r: &[E::ScalarField], srs: &Self::SRS) -> Self::PolyCommitmentProof {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        let length_x = srs.degree_x.log_2();
        let length_y = srs.degree_y.log_2();
        let (x, _) = split_between_x_and_y::<F>(length_x, length_y, r, F::ZERO);
        KZH2::open_1(&poly, C.unwrap().clone(), x.as_slice())
    }

    fn verify(commitment: &Self::Commitment, proof: &Self::PolyCommitmentProof, srs: &Self::SRS, r: &[F], eval: &F) -> Result<(), PCSError> {
        let length_x = srs.degree_x.log_2();
        let length_y = srs.degree_y.log_2();
        let (x, y) = split_between_x_and_y::<F>(length_x, length_y, r, F::ZERO);

        // verify the proof
        KZH2::verify_1(&srs, commitment, proof, x.as_slice(), y.as_slice(), eval);

        Ok(())
    }

    fn setup(max_poly_vars: usize, rng: &mut impl RngCore) -> Result<Self::SRS, Error> {
        let x = max_poly_vars / 2;
        let y = max_poly_vars - x;
        let degree_x = 2usize.pow(x as u32);
        let degree_y = 2usize.pow(y as u32);
        Ok(KZH2::setup_1(degree_x, degree_y, rng))
    } }

#[cfg(test)]
pub mod test {
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{ScalarField, E};
    use crate::nexus_spartan::polycommitments::{PolyCommitmentScheme};
    use crate::kzh::kzh2::KZH2SRS;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

    #[test]
    fn test_end_to_end() {
        let srs: KZH2SRS<E> = MultilinearPolynomial::<ScalarField>::setup(5, &mut thread_rng()).unwrap();

        // random bivariate polynomial
        let polynomial = MultilinearPolynomial::<ScalarField>::rand(3, &mut thread_rng());

        // random points and evaluation
        let x = vec![
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
        ];

        let z = polynomial.evaluate(&x);

        let com = MultilinearPolynomial::commit(&polynomial, &srs);
        let open = MultilinearPolynomial::prove(Option::from(&com), &polynomial, x.as_slice(), &srs);
        MultilinearPolynomial::verify(&com, &open, &srs, x.as_slice(), &z).expect("verification failed");
    }
}

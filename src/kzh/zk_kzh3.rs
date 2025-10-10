use std::marker::PhantomData;
use std::ops::Mul;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::AffineRepr;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::Rng;
use crate::kzh::{pad_at_start, SparseMultilinearPolynomial, KZH};
use crate::kzh::kzh3::{KZH3Aux, KZH3Commitment, KZH3Opening, KZH3, KZH3SRS};
use crate::math::Math;
use crate::nexus_spartan::commitment_traits::ToAffine;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};

/// Define the new struct that encapsulates the functionality of polynomial commitment
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct zkKZH3<E: Pairing> {
    phantom: PhantomData<E>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct zkKZH3SRS<E: Pairing> {
    pub srs: KZH3SRS<E>,
    pub h: E::G1Affine,
}

#[derive(
    Default,
    Clone,
    Debug,
    PartialEq,
    Eq,
    CanonicalSerialize,
    CanonicalDeserialize,
    Derivative
)]
pub struct zkKZH3Commitment<E: Pairing> {
    /// the commitment C to the polynomial
    pub C: E::G1Affine,
}

impl<E: Pairing> ToAffine<E> for zkKZH3Commitment<E> {
    fn to_affine(&self) -> E::G1Affine {
        self.C
    }
}

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for zkKZH3Commitment<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>) {
        Transcript::append_point::<E>(transcript, label, &self.C);
    }
}

#[derive(
    Default,
    Clone,
    Debug,
    PartialEq,
    Eq,
    CanonicalSerialize,
    CanonicalDeserialize,
    Derivative
)]
pub struct zkKZH3Aux<E: Pairing> {
    /// the randomness used for hiding
    pub tau: E::ScalarField,
    /// auxiliary data which is in fact Pedersen commitments to rows of the polynomial
    pub D_x: Vec<E::G1>,
}

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for zkKZH3Aux<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>) {

    }
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct zkKZH3Opening<E: Pairing> {
    pub com_r: E::G1Affine,
    pub r_eval: E::ScalarField,
    pub rho_prime: E::ScalarField,
    pub D_x: Vec<E::G1>,
    pub D_y: Vec<E::G1>,
    pub C_y: E::G1Affine,
    pub f_star: MultilinearPolynomial<E::ScalarField>,
}


impl<E: Pairing> KZH<E> for zkKZH3<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    type Degree = <KZH3<E> as KZH<E>>::Degree;
    type SRS = zkKZH3SRS<E>;
    type Commitment = zkKZH3Commitment<E>;
    type Aux = zkKZH3Aux<E>;
    type Opening = zkKZH3Opening<E>;

    fn split_input<T: Clone>(srs: &Self::SRS, input: &[T], default: T) -> Vec<Vec<T>> {
        KZH3::split_input(&srs.srs, input, default)
    }

    fn get_degree_from_maximum_supported_degree(n: usize) -> Self::Degree {
        KZH3::<E>::get_degree_from_maximum_supported_degree(n)
    }

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS {
        let h = E::G1Affine::rand(rng);
        let srs = KZH3::setup(maximum_degree, rng);

        zkKZH3SRS {
            srs,
            h,
        }
    }

    fn commit<R: Rng>(srs: &Self::SRS, poly: &MultilinearPolynomial<E::ScalarField>, rng: &mut R) -> (Self::Commitment, Self::Aux) {
        let (com, aux) = KZH3::commit(&srs.srs, poly, rng);
        let tau = E::ScalarField::rand(rng);

        (
            zkKZH3Commitment {
                C: (srs.h.mul(&tau) + com.C).into(),
            },
            zkKZH3Aux {
                tau,
                D_x: aux.D_x,
            }
        )
    }

    fn open<R: Rng>(srs: &Self::SRS, input: &[E::ScalarField], com: &Self::Commitment, aux: &Self::Aux, poly: &MultilinearPolynomial<E::ScalarField>, rng: &mut R) -> Self::Opening {
        let num_non_zero = 3 * (srs.srs.degree_y * srs.srs.degree_x * srs.srs.degree_z).square_root();

        // pad the polynomial
        let len = srs.srs.degree_x.log_2() + srs.srs.degree_y.log_2() + srs.srs.degree_z.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        // generate random polynomial r
        let r_poly = SparseMultilinearPolynomial::rand(rng, len, num_non_zero).to_dense();

        // commit to r_poly, this step can be improved by using the sparse polynomial directly to compute this
        let (com_r, aux_r) = zkKZH3::commit(srs, &r_poly, rng);

        // evaluate r_poly at point input
        let padded_input = pad_at_start(input, len);
        let r_eval = r_poly.evaluate(padded_input.as_slice());

        // It is better this transcript is passed as an argument
        let mut transcript = Transcript::<E::ScalarField>::new(b"label");
        transcript.append_point::<E>(b"label", &com_r.C);
        transcript.append_scalar(b"label", &r_eval);
        let alpha = transcript.challenge_scalar(b"alpha");

        let rho_prime = alpha * aux.tau + aux_r.tau;

        let C_lin = (com.C.mul(alpha) + com_r.C - srs.h.mul(rho_prime)).into();

        let new_poly = poly.scale_by_r(alpha) + r_poly;

        let mut new_aux = KZH3Aux { D_x: aux.D_x.clone() };
        new_aux.scale_by_r(&alpha);
        let new_aux = new_aux + KZH3Aux { D_x: aux_r.D_x };

        let open = KZH3::open(&srs.srs, input, &KZH3Commitment { C: C_lin }, &new_aux, &new_poly, rng);

        zkKZH3Opening {
            com_r: com_r.C,
            r_eval,
            rho_prime,
            D_x: open.D_x,
            D_y: open.D_y,
            C_y: open.C_y,
            f_star: open.f_star,
        }
    }

    fn verify(srs: &Self::SRS, input: &[E::ScalarField], output: &E::ScalarField, com: &Self::Commitment, open: &Self::Opening) {
        // It is better this transcript is passed as an argument
        let mut transcript = Transcript::<E::ScalarField>::new(b"label");
        transcript.append_point::<E>(b"label", &open.com_r);
        transcript.append_scalar(b"label", &open.r_eval);
        let alpha = transcript.challenge_scalar(b"alpha");

        let C_lin = (com.C.mul(alpha) + open.com_r - srs.h.mul(open.rho_prime)).into();

        let y_lin = alpha * output + open.r_eval;

        KZH3::verify(&srs.srs, input, &y_lin, &KZH3Commitment { C: C_lin }, &KZH3Opening { D_x: open.D_x.clone(), D_y: open.D_y.clone(), C_y: open.C_y.clone(), f_star: open.f_star.clone() })
    }
}

#[cfg(test)]
pub mod test {
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::{E, ScalarField as F};
    use crate::kzh::KZH;
    use crate::kzh::zk_kzh3::{zkKZH3, zkKZH3SRS};
    use crate::math::Math;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

    #[test]
    fn test_end_to_end() {
        let (degree_x, degree_y, degree_z) = (4usize, 4usize, 4usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        let input: Vec<F> = (0..num_vars)
            .map(|_| F::rand(&mut thread_rng()))
            .collect();

        // build the srs
        let srs: zkKZH3SRS<E> = zkKZH3::setup((degree_x * degree_y * degree_z).log_2(), &mut thread_rng());

        // build a random polynomials
        let polynomial: MultilinearPolynomial<F> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        // evaluate polynomial
        let eval = polynomial.evaluate(input.as_slice());

        // commit to the polynomial
        let (com, aux) = zkKZH3::commit(&srs, &polynomial, &mut thread_rng());

        // open it
        let open = zkKZH3::open(&srs, input.as_slice(), &com, &aux, &polynomial, &mut thread_rng());

        // verify the commit
        zkKZH3::verify(&srs, input.as_slice(), &eval, &com, &open);
    }
}


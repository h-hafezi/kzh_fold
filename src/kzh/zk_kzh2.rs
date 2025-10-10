use crate::kzh::pad_at_start;
use std::marker::PhantomData;
use std::ops::Mul;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::AffineRepr;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::Rng;
use crate::kzh::{SparseMultilinearPolynomial, KZH};
use crate::kzh::kzh2::{KZH2Aux, KZH2Commitment, KZH2Opening, KZH2, KZH2SRS};
use crate::math::Math;
use crate::nexus_spartan::commitment_traits::ToAffine;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};

/// Define the new struct that encapsulates the functionality of polynomial commitment
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct zkKZH2<E: Pairing> {
    phantom: PhantomData<E>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct zkKZH2SRS<E: Pairing> {
    pub srs: KZH2SRS<E>,
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
pub struct zkKZH2Commitment<E: Pairing> {
    /// the commitment C to the polynomial
    pub C: E::G1Affine,
}

impl<E: Pairing> ToAffine<E> for zkKZH2Commitment<E> {
    fn to_affine(&self) -> E::G1Affine {
        self.C
    }
}

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for zkKZH2Commitment<E>
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
pub struct zkKZH2Aux<E: Pairing> {
    /// the randomness used for hiding
    pub tau: E::ScalarField,
    /// auxiliary data which is in fact Pedersen commitments to rows of the polynomial
    pub aux: Vec<E::G1>,
}

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for zkKZH2Aux<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>) {

    }
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct zkKZH2Opening<E: Pairing> {
    pub com_r: E::G1Affine,
    pub r_eval: E::ScalarField,
    pub rho_prime: E::ScalarField,
    pub D_x: Vec<E::G1Affine>,
    pub f_star: MultilinearPolynomial<E::ScalarField>,
}


impl<E: Pairing> KZH<E> for zkKZH2<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    type Degree = <KZH2<E> as KZH<E>>::Degree;
    type SRS = zkKZH2SRS<E>;
    type Commitment = zkKZH2Commitment<E>;
    type Aux = zkKZH2Aux<E>;
    type Opening = zkKZH2Opening<E>;

    fn split_input<T: Clone>(srs: &Self::SRS, input: &[T], default: T) -> Vec<Vec<T>> {
        KZH2::split_input(&srs.srs, input, default)
    }

    fn get_degree_from_maximum_supported_degree(n: usize) -> Self::Degree {
        KZH2::<E>::get_degree_from_maximum_supported_degree(n)
    }

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS {
        let h = E::G1Affine::rand(rng);
        let srs = KZH2::setup(maximum_degree, rng);

        zkKZH2SRS {
            srs,
            h,
        }
    }

    fn commit<R: Rng>(srs: &Self::SRS, poly: &MultilinearPolynomial<E::ScalarField>, rng: &mut R) -> (Self::Commitment, Self::Aux) {
        let (com, aux) = KZH2::commit(&srs.srs, poly, rng);
        let tau = E::ScalarField::rand(rng);

        (
            zkKZH2Commitment {
                C: (srs.h.mul(&tau) + com.C).into(),
            },
            zkKZH2Aux {
                tau,
                aux: aux.aux,
            }
        )
    }

    fn open<R: Rng>(srs: &Self::SRS, input: &[E::ScalarField], com: &Self::Commitment, aux: &Self::Aux, poly: &MultilinearPolynomial<E::ScalarField>, rng: &mut R) -> Self::Opening {
        let num_non_zero = 2 * (srs.srs.degree_y * srs.srs.degree_x).square_root();

        // pad the polynomial
        let len = srs.srs.degree_x.log_2() + srs.srs.degree_y.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        // generate random polynomial r
        let r_poly = SparseMultilinearPolynomial::rand(rng, len, num_non_zero).to_dense();

        // commit to r_poly, this step can be improved by using the sparse polynomial directly to compute this
        let (com_r, aux_r) = zkKZH2::commit(srs, &r_poly, rng);

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

        let mut new_aux = KZH2Aux { aux: aux.aux.clone() };
        new_aux.scale_by_r(&alpha);
        let new_aux = new_aux + KZH2Aux { aux: aux_r.aux };

        let open = KZH2::open(&srs.srs, input, &KZH2Commitment { C: C_lin }, &new_aux, &new_poly, rng);

        zkKZH2Opening {
            com_r: com_r.C,
            r_eval,
            rho_prime,
            D_x: open.D_x,
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
        
        KZH2::verify(&srs.srs, input, &y_lin, &KZH2Commitment { C: C_lin }, &KZH2Opening { D_x: open.D_x.clone(), f_star: open.f_star.clone() })
    }
}



#[cfg(test)]
pub mod test {
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::{E, ScalarField as F};
    use crate::kzh::KZH;
    use crate::kzh::zk_kzh2::{zkKZH2, zkKZH2SRS};
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

    #[test]
    fn test_end_to_end() {
        let srs: zkKZH2SRS<E> = zkKZH2::setup(10, &mut thread_rng());

        // random bivariate polynomial
        let polynomial = MultilinearPolynomial::rand(3 + 5, &mut thread_rng());

        let x: Vec<_> = std::iter::repeat_with(|| F::rand(&mut thread_rng()))
            .take(3)
            .collect();
        let y: Vec<_> = std::iter::repeat_with(|| F::rand(&mut thread_rng()))
            .take(5)
            .collect();

        // Concatenate x and y into input
        let input: Vec<_> = x.into_iter().chain(y.into_iter()).collect();

        let z = polynomial.evaluate(&input);

        // commit to the polynomial
        let (com, aux) = zkKZH2::commit(&srs, &polynomial, &mut thread_rng());

        // open the commitment
        let open = zkKZH2::open(&srs, input.as_slice(), &com, &aux, &polynomial, &mut thread_rng());

        // re compute x and y verify the proof
        zkKZH2::verify(&srs, input.as_slice(), &z, &com, &open);
    }
}


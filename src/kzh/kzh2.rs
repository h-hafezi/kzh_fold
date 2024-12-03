use std::marker::PhantomData;
use crate::kzh::KZH;
use crate::math::Math;
use crate::nexus_spartan::commitment_traits::ToAffine;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ff::{AdditiveGroup, PrimeField, Zero};
use ark_serialize::Valid;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::{Rng};
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use std::ops::{Add, Mul};

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH2SRS<E: Pairing> {
    /// degree_x = 2 ^ length of x variable
    pub degree_x: usize,
    /// degree_y = 2 ^ length of y variable
    pub degree_y: usize,

    pub H_xy: Vec<Vec<E::G1Affine>>,
    pub H_y: Vec<E::G1Affine>,

    pub V_x: Vec<E::G2>,

    pub V_prime: E::G2,
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
pub struct KZH2Commitment<E: Pairing> {
    /// the commitment C to the polynomial
    pub C: E::G1Affine,
    /// auxiliary data which is in fact Pedersen commitments to rows of the polynomial
    pub aux: Vec<E::G1>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH2Opening<E: Pairing> {
    pub D_x: Vec<E::G1Affine>,
    pub f_star_poly: MultilinearPolynomial<E::ScalarField>,
}

/// Define the new struct that encapsulates the functionality of polynomial commitment
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH2<E: Pairing> {
    phantom: PhantomData<E>,
}

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


impl<E: Pairing> KZH<E> for KZH2<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    type Degree = (usize, usize);
    type SRS = KZH2SRS<E>;
    type Commitment = KZH2Commitment<E>;
    type Opening = KZH2Opening<E>;

    /// the function receives an input r and splits into two sub-vectors x and y to be used for PCS
    /// It's used later when we have a constant SRS, and we pad the polynomial so we can commit to it via SRS
    /// This function in fact pads to polynomial inputs by appends necessary zeros and split the input into x and y input
    fn split_input(srs: &Self::SRS, input: &[E::ScalarField]) -> Vec<Vec<E::ScalarField>> {
        let total_length = srs.degree_x.log_2() + srs.degree_y.log_2();

        // If r is smaller than the required length, extend it with zeros at the beginning
        let mut extended_r = input.to_vec();
        if input.len() < total_length {
            let mut zeros = vec![E::ScalarField::ZERO; total_length - input.len()];
            zeros.extend(extended_r);  // Prepend zeros to the beginning
            extended_r = zeros;
        }

        // Split the vector into two parts
        let r_x = extended_r[..srs.degree_x.log_2()].to_vec();
        let r_y = extended_r[srs.degree_x.log_2()..total_length].to_vec();

        vec![r_x, r_y]
    }

    fn get_degree_from_maximum_supported_degree(n: usize) -> (usize, usize) {
        match n % 2 {
            0 => (n / 2, n / 2),
            1 => (n / 2, n / 2 + 1),
            _ => unreachable!(),
        }
    }

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS {
        let (degree_x, degree_y): (usize, usize) = Self::get_degree_from_maximum_supported_degree(maximum_degree);

        let degree_x = 1 << degree_x;
        let degree_y = 1 << degree_y;

        // sample G_0, G_1, ..., G_m generators from group one
        let G1_generator_vec: Vec<_> = (0..degree_y).map(|_| E::G1Affine::rand(rng)).collect();

        // sample V, generator for group two
        let G2_generator = E::G2Affine::rand(rng);

        // sample trapdoors tau_0, tau_1, ..., tau_n, alpha
        let tau: Vec<E::ScalarField> = (0..degree_x).map(|_| E::ScalarField::rand(rng)).collect();

        let alpha = E::ScalarField::rand(rng);

        // generate matrix_H
        let matrix_H: Vec<Vec<_>> = (0..degree_x).into_par_iter()
            .map(|i| {
                let mut row = Vec::new();
                for j in 0..degree_y {
                    let g = G1_generator_vec[j].mul(tau[i]);
                    row.push(g.into());
                }
                row
            }).collect();

        let vec_H: Vec<_> = (0..degree_y).map(|j| G1_generator_vec[j].mul(alpha).into()).collect();
        let vec_V: Vec<_> = (0..degree_x).map(|j| G2_generator.mul(tau[j])).collect();

        // generate V_prime
        let V_prime = G2_generator.mul(alpha);

        // return the output
        KZH2SRS {
            degree_x,
            degree_y,
            H_xy: matrix_H,
            H_y: vec_H,
            V_x: vec_V,
            V_prime,
        }
    }

    fn commit(srs: &Self::SRS, poly: &MultilinearPolynomial<E::ScalarField>) -> Self::Commitment {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        KZH2Commitment {
            C: {
                // Collect all points and scalars into single vectors
                let mut base = Vec::new();
                let mut scalar = Vec::new();

                for i in 0..srs.degree_x {
                    // Collect points from matrix_H
                    base.extend_from_slice(srs.H_xy[i].as_slice());
                    // Collect corresponding scalars from partial evaluations
                    scalar.extend_from_slice(poly.get_partial_evaluation_for_boolean_input(i, srs.degree_y).as_slice());
                }

                E::G1::msm_unchecked(&base, &scalar).into_affine()
            },
            aux: (0..srs.degree_x)
                .into_par_iter() // Parallelize the D^{(x)} computation
                .map(|i| {
                    E::G1::msm_unchecked(
                        srs.H_y.as_slice(),
                        poly.get_partial_evaluation_for_boolean_input(i, srs.degree_y).as_slice(),
                    )
                })
                .collect::<Vec<_>>(),
        }
    }

    fn open(
        srs: &Self::SRS,
        input: &[E::ScalarField],
        com: &Self::Commitment,
        poly: &MultilinearPolynomial<E::ScalarField>,
    ) -> Self::Opening {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        let split_input = Self::split_input(&srs, input);

        KZH2Opening {
            D_x: com.aux.clone().into_iter().map(|g| g.into()).collect(),
            f_star_poly: poly.partial_evaluation(split_input[0].as_slice()),
        }
    }

    fn verify(srs: &Self::SRS, input: &[E::ScalarField], output: &E::ScalarField, com: &Self::Commitment, open: &Self::Opening) {
        let split_input = Self::split_input(&srs, input);

        // Step 1: pairing check
        // Combine the pairings into a single multi-pairing
        let g1_elems: Vec<_> = std::iter::once(com.C.clone())
            .chain(open.D_x.iter().map(|g1| (E::G1Affine::zero() - g1).into()))
            .collect();

        let g2_elems: Vec<_> = std::iter::once(srs.V_prime.clone())
            .chain(srs.V_x.iter().cloned())
            .collect();

        // Perform the combined pairing check
        E::multi_pairing(&g1_elems, &g2_elems).check().unwrap();

        // Step 2: MSM check
        let negated_eq_evals: Vec<_> = EqPolynomial::new(split_input[0].clone())
            .evals()
            .into_iter()
            .map(|scalar| -scalar)
            .collect();

        let scalars: Vec<_> = open.f_star_poly.evaluation_over_boolean_hypercube
            .iter()
            .chain(negated_eq_evals.iter())
            .cloned()
            .collect();

        let bases: Vec<_> = srs.H_y.iter()
            .chain(open.D_x.iter())
            .cloned()
            .collect();

        assert!(E::G1::msm_unchecked(&bases, &scalars).is_zero());

        // Step 3: complete poly eval
        let y_expected = open.f_star_poly.evaluate(split_input[1].as_slice());
        assert_eq!(y_expected, *output);
    }
}

/// the function receives an input r and splits into two sub-vectors x and y to be used for PCS
/// It's used later when we have a constant SRS, and we pad the polynomial so we can commit to it via SRS
/// This function in fact pads to polynomial inputs by appends necessary zeros and split the input into x and y input
pub fn split_between_x_and_y<T: Clone>(x_length: usize, y_length: usize, r: &[T], zero: T) -> (Vec<T>, Vec<T>) {
    let total_length = x_length + y_length;

    // If r is smaller than the required length, extend it with zeros at the beginning
    let mut extended_r = r.to_vec();
    if r.len() < total_length {
        let mut zeros = vec![zero; total_length - r.len()];
        zeros.extend(extended_r);  // Prepend zeros to the beginning
        extended_r = zeros;
    }

    // Split the vector into two parts
    let r_x = extended_r[..x_length].to_vec();
    let r_y = extended_r[x_length..total_length].to_vec();

    (r_x, r_y)
}


impl<E: Pairing> KZH2Commitment<E> {
    /// Scales the commitment and its auxiliary elements by a scalar `r`
    pub fn scale_by_r(&mut self, r: &E::ScalarField) {
        // Scale the main commitment C by r
        let scaled_C = self.C.mul(r); // G1Affine -> G1Projective when multiplied by scalar

        // Scale each element in the aux vector by r
        let scaled_aux: Vec<E::G1> = self.aux.iter()
            .map(|element| element.mul(r))  // Multiply each element in aux by r
            .collect();

        // Update the commitment with the scaled values
        self.C = scaled_C.into_affine();  // Convert back to G1Affine after multiplication
        self.aux = scaled_aux;
    }
}


impl<E: Pairing> Add for KZH2Commitment<E> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        // Ensure both commitments have the same size in aux vectors
        assert_eq!(self.aux.len(), other.aux.len(), "Aux vectors must have the same length");

        // Add the main commitment points C
        let new_C = (self.C + other.C).into_affine();

        // Add the corresponding elements in the aux vector
        let new_aux: Vec<E::G1> = self.aux.iter()
            .zip(other.aux.iter())
            .map(|(a, b)| *a + *b)
            .collect();

        // Return a new Commitment with the resulting sums
        KZH2Commitment {
            C: new_C,
            aux: new_aux,
        }
    }
}

#[cfg(test)]
pub mod test {
    use ark_serialize::CanonicalSerialize;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{ScalarField as F, E};
    use crate::kzh::kzh2::{KZH2, KZH2SRS};
    use crate::kzh::KZH;
    use crate::math::Math;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

    #[test]
    fn test_end_to_end() {
        let srs: KZH2SRS<E> = KZH2::setup(10, &mut thread_rng());

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
        let com = KZH2::commit(&srs, &polynomial);

        // open the commitment
        let open = KZH2::open(&srs, input.as_slice(), &com, &polynomial);

        // re compute x and y verify the proof
        KZH2::verify(&srs, input.as_slice(), &z, &com, &open);
    }


    /// Given f(x) and g(x) and their KZH commitments F and G.
    /// This test computes p(x) = f(x) + r * g(x),
    /// and checks that its commitment is P = F + r*G
    ///
    /// Prover sends F,G
    /// Verifier responds with r, rho
    /// Prover sends p(rho), f(rho), g(rho), proof_P_at_rho
    /// Verifier checks that p(rho) = f(rho) + r * g(rho)
    /// and the proof verifies using P = F + r * G
    #[test]
    fn test_homomorphism() {
        let num_vars = 8; // degree_x.log_2() + degree_y.log_2()

        let srs: KZH2SRS<E> = KZH2::setup(8, &mut thread_rng());

        let f_x: MultilinearPolynomial<F> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());
        let g_x: MultilinearPolynomial<F> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        let F = KZH2::commit(&srs, &f_x);
        let G = KZH2::commit(&srs, &g_x);

        // Verifier's challenge: for poly batching
        let r = F::rand(&mut thread_rng());
        // Verifier's challenge: evaluation point
        let rho = vec![F::rand(&mut thread_rng()); num_vars];

        // Compute p(x) = f(x) + r * g(x)
        let mut r_times_g_x = g_x.clone();
        r_times_g_x.scalar_mul(&r);
        let p_x = f_x.clone() + r_times_g_x;
        let P = KZH2::commit(&srs, &p_x);

        // Open p_x at rho
        let proof_P_at_rho = KZH2::open(&srs, rho.as_slice(), &P, &p_x);
        let p_at_rho = p_x.evaluate(&rho);

        // Verifier:
        assert_eq!(p_at_rho, f_x.evaluate(&rho) + r * g_x.evaluate(&rho));

        // Verifier: compute P = F + r*G
        let mut r_times_G = G.clone();
        r_times_G.scale_by_r(&r);
        let P_verifier = F + r_times_G;

        KZH2::verify(&srs, rho.as_slice(), &p_at_rho, &P_verifier, &proof_P_at_rho);
    }
}


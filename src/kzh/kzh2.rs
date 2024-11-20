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
use rand::{Rng, RngCore};
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
pub struct KZH2OpeningProof<E: Pairing> {
    pub vec_D: Vec<E::G1Affine>,
    pub f_star_poly: MultilinearPolynomial<E::ScalarField>,
}

/// Define the new struct that encapsulates the functionality of polynomial commitment
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH2<E: Pairing> {
    phatom: PhantomData<E>,
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
    type Opening = KZH2OpeningProof<E>;

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

    fn get_degree_from_maximum_supported_degree(n: usize) -> Self::Degree {
        (n / 2 + 1, n / 2 + 1)
    }

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS {
        let (degree_x, degree_y): (usize, usize) = Self::get_degree_from_maximum_supported_degree(maximum_degree);

        let degree_x = 1 << degree_x;
        let degree_y = 1 << degree_y;

        // sample G_0, G_1, ..., G_m generators from group one
        let G1_generator_vec = {
            let mut elements = Vec::new();
            for _ in 0..degree_y {
                elements.push(E::G1Affine::rand(rng));
            }
            elements
        };

        // sample V, generator for group two
        let G2_generator = E::G2Affine::rand(rng);

        // sample trapdoors tau_0, tau_1, ..., tau_n, alpha
        let tau = {
            let mut elements = Vec::new();
            for _ in 0..degree_x {
                elements.push(E::ScalarField::rand(rng));
            }
            elements
        };

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
            })
            .collect();

        // generate vec_H
        let vec_H = {
            let mut vec_h = Vec::new();
            for j in 0..degree_y {
                vec_h.push(G1_generator_vec[j].mul(alpha).into());
            }
            vec_h
        };

        // generate vec_V
        let vec_V = {
            let mut vec_h = Vec::new();
            for j in 0..degree_x {
                vec_h.push(G2_generator.mul(tau[j]));
            }
            vec_h
        };

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

        KZH2OpeningProof {
            vec_D: {
                let mut vec = Vec::new();
                for g in com.aux.clone() {
                    vec.push(g.into());
                }
                vec
            },
            f_star_poly: poly.partial_evaluation(split_input[0].as_slice()),
        }
    }

    fn verify(srs: &Self::SRS, input: &[E::ScalarField], output: &E::ScalarField, com: &Self::Commitment, open: &Self::Opening) {
        let split_input = Self::split_input(&srs, input);

        // Step 1: pairing check
        // Combine the pairings into a single multi-pairing
        let mut g1_elems: Vec<E::G1Affine> = Vec::with_capacity(1 + open.vec_D.len());
        g1_elems.push(com.C.clone());
        for g1 in &open.vec_D {
            let g1_neg: E::G1Affine = (E::G1Affine::zero() - g1).into();
            g1_elems.push(g1_neg);
        }

        let mut g2_elems = Vec::with_capacity(1 + srs.V_x.len());
        g2_elems.push(srs.V_prime.clone());
        g2_elems.extend_from_slice(&srs.V_x);

        // Perform the combined pairing check
        let pairing_product = E::multi_pairing(&g1_elems, &g2_elems);
        pairing_product.check().unwrap();

        // Step 2: MSM check
        // Combine the two MSMs into one
        let mut negated_eq_evals = EqPolynomial::new(split_input[0].clone()).evals();
        for scalar in &mut negated_eq_evals {
            *scalar = -*scalar;
        }

        let mut scalars = Vec::with_capacity(
            open.f_star_poly.evaluation_over_boolean_hypercube.len() + negated_eq_evals.len(),
        );
        scalars.extend_from_slice(&open.f_star_poly.evaluation_over_boolean_hypercube);
        scalars.extend_from_slice(&negated_eq_evals);

        let mut bases = Vec::with_capacity(srs.H_y.len() + open.vec_D.len());
        bases.extend_from_slice(&srs.H_y);
        bases.extend_from_slice(&open.vec_D);

        let msm_result = E::G1::msm_unchecked(&bases, &scalars);
        assert!(msm_result.is_zero());


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


/// all functions assume that poly size is already compatible with srs size, if not there's an interface that does padding
/// in the beginning in kzh.rs
impl<E: Pairing> KZH2<E> {
    /// set up the PCS srs
    pub fn setup_1<T: RngCore>(degree_x: usize, degree_y: usize, rng: &mut T) -> KZH2SRS<E> {
        // sample G_0, G_1, ..., G_m generators from group one
        let G1_generator_vec = {
            let mut elements = Vec::new();
            for _ in 0..degree_y {
                elements.push(E::G1Affine::rand(rng));
            }
            elements
        };

        // sample V, generator for group two
        let G2_generator = E::G2Affine::rand(rng);

        // sample trapdoors tau_0, tau_1, ..., tau_n, alpha
        let tau = {
            let mut elements = Vec::new();
            for _ in 0..degree_x {
                elements.push(E::ScalarField::rand(rng));
            }
            elements
        };

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
            })
            .collect();

        // generate vec_H
        let vec_H = {
            let mut vec_h = Vec::new();
            for j in 0..degree_y {
                vec_h.push(G1_generator_vec[j].mul(alpha).into());
            }
            vec_h
        };

        // generate vec_V
        let vec_V = {
            let mut vec_h = Vec::new();
            for j in 0..degree_x {
                vec_h.push(G2_generator.mul(tau[j]));
            }
            vec_h
        };

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

    pub fn commit_1(srs: &KZH2SRS<E>, poly: &MultilinearPolynomial<E::ScalarField>) -> KZH2Commitment<E> {
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

    /// Creates a KZH proof for p(x,y) = z.
    /// This function does not actually need y, so we only get the left half of the eval point.
    pub fn open_1(poly: &MultilinearPolynomial<E::ScalarField>, com: KZH2Commitment<E>, x: &[E::ScalarField]) -> KZH2OpeningProof<E> {
        KZH2OpeningProof {
            vec_D: {
                let mut vec = Vec::new();
                for g in com.aux {
                    vec.push(g.into());
                }
                vec
            },
            f_star_poly: poly.partial_evaluation(x),
        }
    }

    pub fn verify_1(srs: &KZH2SRS<E>,
                                C: &KZH2Commitment<E>,
                                proof: &KZH2OpeningProof<E>,
                                x: &[E::ScalarField],
                                y: &[E::ScalarField],
                                z: &E::ScalarField,
    ) {
        // Step 1: pairing check
        // Combine the pairings into a single multi-pairing
        let mut g1_elems: Vec<E::G1Affine> = Vec::with_capacity(1 + proof.vec_D.len());
        g1_elems.push(C.C.clone());
        for g1 in &proof.vec_D {
            let g1_neg: E::G1Affine = (E::G1Affine::zero() - g1).into();
            g1_elems.push(g1_neg);
        }

        let mut g2_elems = Vec::with_capacity(1 + srs.V_x.len());
        g2_elems.push(srs.V_prime.clone());
        g2_elems.extend_from_slice(&srs.V_x);

        // Perform the combined pairing check
        let pairing_product = E::multi_pairing(&g1_elems, &g2_elems);
        pairing_product.check().unwrap();

        // Step 2: MSM check
        // Combine the two MSMs into one
        let mut negated_eq_evals = EqPolynomial::new(x.to_vec()).evals();
        for scalar in &mut negated_eq_evals {
            *scalar = -*scalar;
        }

        let mut scalars = Vec::with_capacity(
            proof.f_star_poly.evaluation_over_boolean_hypercube.len() + negated_eq_evals.len(),
        );
        scalars.extend_from_slice(&proof.f_star_poly.evaluation_over_boolean_hypercube);
        scalars.extend_from_slice(&negated_eq_evals);

        let mut bases = Vec::with_capacity(srs.H_y.len() + proof.vec_D.len());
        bases.extend_from_slice(&srs.H_y);
        bases.extend_from_slice(&proof.vec_D);

        let msm_result = E::G1::msm_unchecked(&bases, &scalars);
        assert!(msm_result.is_zero());


        // Step 3: complete poly eval
        let y_expected = proof.f_star_poly.evaluate(y);
        assert_eq!(y_expected, *z);
    }
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
    use std::cmp::min;

    use ark_ec::pairing::Pairing;
    use ark_ff::AdditiveGroup;
    use ark_serialize::CanonicalSerialize;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{ScalarField, E};
    use crate::kzh::kzh2::{split_between_x_and_y, KZH2, KZH2SRS};
    use crate::kzh::KZH;
    use crate::math::Math;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;

    #[test]
    fn test_setup() {
        let degree_y = 4usize;
        let degree_x = 4usize;
        let srs: KZH2SRS<E> = KZH2::setup_1(degree_x, degree_y, &mut thread_rng());

        // asserting the sizes
        assert_eq!(srs.degree_y, degree_y);
        assert_eq!(srs.degree_x, degree_x);
        assert_eq!(srs.H_y.len(), degree_y);
        assert_eq!(srs.V_x.len(), degree_x);
        assert_eq!(srs.H_xy.len(), degree_x);
        assert_eq!(srs.H_xy[0].len(), degree_y);

        // checking pairing equalities
        // e(H[j, i], V[i]) = e(G_i^{tau_j}, V^{tau_i}) = e(H[i, i], V[j])
        for i in 0..min(degree_y, degree_x) {
            for j in 0..min(degree_y, degree_x) {
                let p1 = E::pairing(srs.H_xy[j][i], srs.V_x[i]);
                let p2 = E::pairing(srs.H_xy[i][i], srs.V_x[j]);
                assert_eq!(p1, p2);
            }
        }
    }

    #[test]
    fn test_end_to_end() {
        let (degree_x, degree_y) = (8usize, 32usize);
        let srs: KZH2SRS<E> = KZH2::setup_1(degree_x, degree_y, &mut thread_rng());

        // testing srs functions
        assert_eq!(3, srs.degree_x);
        assert_eq!(5, srs.degree_y);

        // ********************** test the input padding function split_between_x_and_y **********************
        let mut r = vec![
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
        ];
        r.extend(vec![ScalarField::ZERO; 3]);
        let x = r[0..3].to_vec();
        let y = r[3..].to_vec();

        // do the split and assert equality
        let length_x = srs.degree_x;
        let length_y = srs.degree_y;
        let (x_new, y_new) = split_between_x_and_y::<ScalarField>(
            length_x,
            length_y,
            r.as_slice(),
            ScalarField::ZERO,
        );
        assert_eq!(x_new, x);
        assert_eq!(y_new, y);

        // ********************** test the input padding function split_between_x_and_y **********************

        // random bivariate polynomial
        let polynomial = MultilinearPolynomial::rand(3 + 5, &mut thread_rng());

        // random points and evaluation
        let x = vec![
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
        ];
        let y = vec![
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
        ];
        // concat inputs x and y, to evaluate the function
        let input = {
            let mut res = vec![];
            res.extend(x.clone());
            res.extend(y.clone());
            res
        };

        let z = polynomial.evaluate(&input);

        // commit to the polynomial
        let com = KZH2::commit_1(&srs, &polynomial);

        // open the commitment
        let open = KZH2::open_1(&polynomial, com.clone(), &x);

        // re compute x and y verify the proof
        KZH2::verify_1(&srs, &com, &open, &x, &y, &z);
    }

    #[test]
    fn test_end_to_end_2() {
        let (degree_x, degree_y) = (8usize, 32usize);
        let srs: KZH2SRS<E> = KZH2::setup_1(degree_x, degree_y, &mut thread_rng());

        // testing srs functions
        assert_eq!(3, srs.degree_x.log_2());
        assert_eq!(5, srs.degree_y.log_2());

        // random bivariate polynomial
        let polynomial = MultilinearPolynomial::rand(3 + 5, &mut thread_rng());

        // concat inputs x and y, to evaluate the function
        // random points and evaluation
        let x = vec![
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
        ];
        let y = vec![
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
            ScalarField::rand(&mut thread_rng()),
        ];
        // concat inputs x and y, to evaluate the function
        let input = {
            let mut res = vec![];
            res.extend(x.clone());
            res.extend(y.clone());
            res
        };

        let z = polynomial.evaluate(&input);

        // commit to the polynomial
        let com = KZH2::commit(&srs, &polynomial);
        // commit to the polynomial
        let com1 = KZH2::commit_1(&srs, &polynomial);
        assert_eq!(com, com1);
        // open the commitment
        let open = KZH2::open(&srs, input.as_slice(), &com, &polynomial);
        // open the commitment
        let open1 = KZH2::open_1(&polynomial, com.clone(), &x);
        assert_eq!(open, open1);

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
        let degree_x = 16usize;
        let degree_y = 16usize;
        let num_vars = 8; // degree_x.log_2() + degree_y.log_2()

        let srs: KZH2SRS<E> = KZH2::setup_1(degree_x, degree_y, &mut thread_rng());

        let f_x: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());
        let g_x: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        let F = KZH2::commit_1(&srs, &f_x);
        let G = KZH2::commit_1(&srs, &g_x);

        // Verifier's challenge: for poly batching
        let r = ScalarField::rand(&mut thread_rng());
        // Verifier's challenge: evaluation point
        let rho = vec![ScalarField::rand(&mut thread_rng()); num_vars];
        // Split rho in half
        assert_eq!(rho.len() % 2, 0);
        let mid = rho.len() / 2;
        let (rho_first_half, rho_second_half) = rho.split_at(mid);

        // Compute p(x) = f(x) + r * g(x)
        let mut r_times_g_x = g_x.clone();
        r_times_g_x.scalar_mul(&r);
        let p_x = f_x.clone() + r_times_g_x;
        let P = KZH2::commit_1(&srs, &p_x);

        // Open p_x at rho
        let proof_P_at_rho = KZH2::open_1(&p_x, P.clone(), &rho_first_half);
        let p_at_rho = p_x.evaluate(&rho);

        // Verifier:
        assert_eq!(p_at_rho, f_x.evaluate(&rho) + r * g_x.evaluate(&rho));

        // Verifier: compute P = F + r*G
        let mut r_times_G = G.clone();
        r_times_G.scale_by_r(&r);
        let P_verifier = F + r_times_G;

        KZH2::verify_1(&srs, &P_verifier, &proof_P_at_rho, rho_first_half, rho_second_half, &p_at_rho);
    }

    #[test]
    fn count_witness() {
        let degrees = vec![(4, 4), (8, 8), (16, 16), (32, 32), (64, 64)];
        for (degree_x, degree_y) in degrees {
            let srs: KZH2SRS<E> = KZH2::setup_1(degree_x, degree_y, &mut thread_rng());
            // random bivariate polynomial
            let polynomial = MultilinearPolynomial::rand(
                srs.degree_x + srs.degree_y,
                &mut thread_rng(),
            );
            let com = KZH2::commit_1(&srs, &polynomial);

            // random points and evaluation
            let x = {
                let mut res = Vec::new();
                for _ in 0..srs.degree_x {
                    res.push(ScalarField::rand(&mut thread_rng()));
                }
                res
            };

            let open = KZH2::open_1(&polynomial, com.clone(), &x);
            let degree = degree_x * degree_y;
            println!("witness length in bytes: {} for degree {degree}",
                     open.compressed_size()
            );
        }
    }
}


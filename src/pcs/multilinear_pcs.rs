use std::iter::Sum;
use std::ops::{Add, Mul};

use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ff::AdditiveGroup;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::RngCore;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use crate::nexus_spartan::math::Math;
use crate::polynomial::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::MultilinearPolynomial;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct SRS<E: Pairing> {
    pub degree_x: usize,
    pub degree_y: usize,
    pub matrix_H: Vec<Vec<E::G1Affine>>,
    pub vec_H: Vec<E::G1Affine>,
    pub vec_V: Vec<E::G2>,
    pub V_prime: E::G2,
}

#[derive(Default, Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct Commitment<E: Pairing> {
    pub C: E::G1Affine,
    pub aux: Vec<E::G1>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct OpeningProof<E: Pairing> {
    pub vec_D: Vec<E::G1Affine>,
    pub f_star_poly: MultilinearPolynomial<E::ScalarField>,
}

// Define the new struct that encapsulates the functionality of polynomial commitment
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct PolyCommit<E: Pairing> {
    pub srs: SRS<E>,
}

impl<E: Pairing> SRS<E> {
    pub fn get_x_length(&self) -> usize {
        self.degree_x.log_2()
    }

    pub fn get_y_length(&self) -> usize {
        self.degree_y.log_2()
    }

    pub fn split_between_x_and_y(&self, r: &[E::ScalarField]) -> (Vec<E::ScalarField>, Vec<E::ScalarField>) {
        let x_length = self.get_x_length();
        let y_length = self.get_y_length();
        let total_length = x_length + y_length;

        // If r is smaller than the required length, extend it with zeros
        let mut extended_r = r.to_vec();
        if r.len() < total_length {
            extended_r.extend(vec![E::ScalarField::ZERO; total_length - r.len()]);
        }

        // Split the vector into two parts
        let r_x = extended_r[..x_length].to_vec();
        let r_y = extended_r[x_length..total_length].to_vec();

        (r_x, r_y)
    }
}

impl<E: Pairing> PolyCommit<E> {
    pub fn setup<T: RngCore>(degree_x: usize, degree_y: usize, rng: &mut T) -> SRS<E> {
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
        return SRS {
            degree_x,
            degree_y,
            matrix_H,
            vec_H,
            vec_V,
            V_prime,
        };
    }

    pub fn commit(&self, poly: &MultilinearPolynomial<E::ScalarField>) -> Commitment<E> {
        Commitment {
            C: E::G1::sum((0..self.srs.degree_x)
                .map(|i| {
                    E::G1::msm_unchecked(
                        self.srs.matrix_H[i].as_slice(),
                        poly.get_partial_evaluation_for_boolean_input(i, self.srs.degree_y).as_slice(),
                    )
                })
                .collect::<Vec<_>>()
                .iter()
            ).into_affine(),
            aux: (0..self.srs.degree_x)
                .map(|i| {
                    E::G1::msm_unchecked(
                        self.srs.vec_H.as_slice(),
                        poly.get_partial_evaluation_for_boolean_input(i, self.srs.degree_y).as_slice(),
                    )
                })
                .collect::<Vec<_>>(),
        }
    }

    /// Creates a KZH proof for p(x,y) = z.
    /// This function does not actually need y, so we only get the left half of the eval point.
    pub fn open(&self, poly: &MultilinearPolynomial<E::ScalarField>, com: Commitment<E>, x: &[E::ScalarField]) -> OpeningProof<E> {
        OpeningProof {
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

    pub fn verify(&self,
                  C: &Commitment<E>,
                  proof: &OpeningProof<E>,
                  x: &[E::ScalarField],
                  y: &[E::ScalarField],
                  z: &E::ScalarField,
    ) -> bool {
        // first condition
        let pairing_rhs = E::multi_pairing(proof.vec_D.clone(), &self.srs.vec_V);
        let pairing_lhs = E::pairing(&C.C, &self.srs.V_prime);

        // second condition
        let msm_lhs = E::G1::msm_unchecked(&self.srs.vec_H, proof
            .f_star_poly
            .evaluation_over_boolean_hypercube.as_slice(),
        );
        let msm_rhs = E::G1::msm_unchecked(&proof.vec_D, &EqPolynomial::new(x.to_vec()).evals());

        // third condition
        let y_expected = proof.f_star_poly.evaluate(y);

        // checking all three conditions
        return (pairing_lhs == pairing_rhs) && (msm_lhs == msm_rhs) && (y_expected == *z);
    }
}


impl<E: Pairing> Commitment<E> {
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


impl<E: Pairing> Add for Commitment<E> {
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
        Commitment {
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
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{E, ScalarField};
    use crate::pcs::multilinear_pcs::{PolyCommit, SRS};
    use crate::polynomial::multilinear_poly::MultilinearPolynomial;

    #[test]
    fn test_setup() {
        let degree_y = 4usize;
        let degree_x = 4usize;
        let srs: SRS<E> = PolyCommit::<E>::setup(degree_x, degree_y, &mut thread_rng());

        // asserting the sizes
        assert_eq!(srs.degree_y, degree_y);
        assert_eq!(srs.degree_x, degree_x);
        assert_eq!(srs.vec_H.len(), degree_y);
        assert_eq!(srs.vec_V.len(), degree_x);
        assert_eq!(srs.matrix_H.len(), degree_x);
        assert_eq!(srs.matrix_H[0].len(), degree_y);

        // checking pairing equalities
        // e(H[j, i], V[i]) = e(G_i^{tau_j}, V^{tau_i}) = e(H[i, i], V[j])
        for i in 0..min(degree_y, degree_x) {
            for j in 0..min(degree_y, degree_x) {
                let p1 = E::pairing(srs.matrix_H[j][i], srs.vec_V[i]);
                let p2 = E::pairing(srs.matrix_H[i][i], srs.vec_V[j]);
                assert_eq!(p1, p2);
            }
        }
    }

    #[test]
    fn test_end_to_end() {
        let degree_x = 8usize;
        let degree_y = 32usize;
        let srs: SRS<E> = PolyCommit::<E>::setup(degree_x, degree_y, &mut thread_rng());

        // testing srs functions
        assert_eq!(3, srs.get_x_length());
        assert_eq!(5, srs.get_y_length());

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
        let (x_new, y_new) = srs.split_between_x_and_y(r.as_slice());
        assert_eq!(x_new, x);
        assert_eq!(y_new, y);


        // define the polynomial commitment
        let poly_commit: PolyCommit<E> = PolyCommit { srs };

        // random bivariate polynomial
        let mut polynomial = MultilinearPolynomial::rand(3 + 5, &mut thread_rng());

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
        let concat = {
            let mut res = vec![];
            res.extend(x.clone());
            res.extend(y.clone());
            res
        };

        let z = polynomial.evaluate(&concat);

        // commit to the polynomial
        let com = poly_commit.commit(&polynomial);

        // open the commitment
        let open = poly_commit.open(&polynomial, com.clone(), &x);

        // verify the proof
        assert!(poly_commit.verify(&com, &open, &x, &y, &z));
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

        let srs: SRS<E> = PolyCommit::<E>::setup(degree_x, degree_y, &mut thread_rng());

        // define the polynomial commitment
        let poly_commit: PolyCommit<E> = PolyCommit { srs };

        let f_x: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());
        let g_x: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        let F = poly_commit.commit(&f_x);
        let G = poly_commit.commit(&g_x);

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
        let P = poly_commit.commit(&p_x);

        // Open p_x at rho
        let proof_P_at_rho = poly_commit.open(&p_x, P.clone(), &rho_first_half);
        let p_at_rho = p_x.evaluate(&rho);

        // Verifier:
        assert_eq!(p_at_rho, f_x.evaluate(&rho) + r * g_x.evaluate(&rho));

        // Verifier: compute P = F + r*G
        let mut r_times_G = G.clone();
        r_times_G.scale_by_r(&r);
        let P_verifier = F + r_times_G;

        assert!(poly_commit.verify(&P_verifier, &proof_P_at_rho, rho_first_half, rho_second_half, &p_at_rho));
    }
}


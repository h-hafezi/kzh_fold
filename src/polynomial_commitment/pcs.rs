/// A sqrt(n) polynomial commitment scheme
use std::{mem, ops::Mul};
use std::iter::Sum;
use std::ops::AddAssign;
use std::time::Instant;

use ark_ec::{AffineRepr, CurveGroup};
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ff::{AdditiveGroup, PrimeField, UniformRand};
use rand::RngCore;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use crate::polynomial::bivariate_polynomial::bivariate_poly::BivariatePolynomial;
use crate::polynomial::bivariate_polynomial::lagrange_basis::LagrangeBasis;
use crate::polynomial::bivariate_polynomial::univariate_poly::UnivariatePolynomial;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SRS<E: Pairing> {
    pub degree_x: usize,
    pub degree_y: usize,
    pub matrix_H: Vec<Vec<E::G1Affine>>,
    pub vec_H: Vec<E::G1Affine>,
    pub vec_V: Vec<E::G2>,
    pub V_prime: E::G2,
}

impl<E: Pairing> SRS<E> {
    pub fn size_of(&self) -> usize {
        let mut total_size = 0;

        // Size of the degree_x and degree_y fields
        total_size += mem::size_of::<usize>() * 2;

        // Size of the matrix_H Vec<Vec<E::G1Affine>>
        total_size += mem::size_of::<Vec<Vec<E::G1Affine>>>();
        for inner_vec in &self.matrix_H {
            total_size += mem::size_of::<Vec<E::G1Affine>>();
            total_size += inner_vec.capacity() * mem::size_of::<E::G1Affine>();
        }

        // Size of the vec_H Vec<E::G1Affine>
        total_size += mem::size_of::<Vec<E::G1Affine>>();
        total_size += self.vec_H.capacity() * mem::size_of::<E::G1Affine>();

        // Size of the vec_V Vec<E::G2>
        total_size += mem::size_of::<Vec<E::G2>>();
        total_size += self.vec_V.capacity() * mem::size_of::<E::G2>();

        // Size of the V_prime field
        total_size += mem::size_of::<E::G2>();

        total_size
    }
}


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<E: Pairing> {
    pub C: E::G1Affine,
    pub aux: Vec<E::G1>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpeningProof<E: Pairing> {
    pub vec_D: Vec<E::G1Affine>,
    pub f_star_poly: UnivariatePolynomial<E::ScalarField>,
}

// Define the new struct that encapsulates the functionality of polynomial commitment
pub struct PolyCommit<E: Pairing> {
    pub srs: SRS<E>,
}

pub trait PolyCommitTrait<E: Pairing> {
    fn setup<T: RngCore>(n: usize, m: usize, rng: &mut T) -> SRS<E>;

    fn commit(&self, poly: &BivariatePolynomial<E::ScalarField>) -> Commitment<E>;

    fn open(&self,
            poly: &BivariatePolynomial<E::ScalarField>,
            com: Commitment<E>,
            b: &E::ScalarField,
    ) -> OpeningProof<E>;

    fn verify(&self,
              lagrange_x: LagrangeBasis<E::ScalarField>,
              C: &Commitment<E>,
              proof: &OpeningProof<E>,
              b: &E::ScalarField,
              c: &E::ScalarField,
              y: &E::ScalarField,
    ) -> bool;
}

impl<E: Pairing> PolyCommitTrait<E> for PolyCommit<E> {
    fn setup<T: RngCore>(degree_x: usize, degree_y: usize, rng: &mut T) -> SRS<E> {
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

    fn commit(&self, poly: &BivariatePolynomial<E::ScalarField>) -> Commitment<E> {
        Commitment {
            C: E::G1::sum((0..self.srs.degree_x).into_par_iter()
                          .map(|i| {
                              let start = i * poly.degree_y;
                              let end = start + poly.degree_y;
                              E::G1::msm_unchecked(
                                  self.srs.matrix_H[i].as_slice(),
                                  &poly.evaluations[start..end], // Slice the flattened vector to get the row
                              )
                          })
                          .collect::<Vec<_>>()
                          .iter()
            ).into_affine(),

            aux: (0..self.srs.degree_x).into_par_iter()
                .map(|i| {
                    let start = i * poly.degree_y;
                    let end = start + poly.degree_y;
                    E::G1::msm_unchecked(
                        self.srs.vec_H.as_slice(),
                        &poly.evaluations[start..end], // Slice the flattened vector to get the row
                    )
                })
                .collect::<Vec<_>>(),
        }
    }

    /// Create opening proof for f(b,c)
    fn open(&self, poly: &BivariatePolynomial<E::ScalarField>, com: Commitment<E>, b: &E::ScalarField) -> OpeningProof<E> {
        OpeningProof {
            vec_D: {
                let mut vec = Vec::new();
                for g in com.aux {
                    vec.push(g.into());
                }
                vec
            },
            f_star_poly: poly.partially_evaluate_at_x(b),
        }
    }

    fn verify(&self,
              lagrange_x: LagrangeBasis<E::ScalarField>,
              C: &Commitment<E>,
              proof: &OpeningProof<E>,
              b: &E::ScalarField,
              c: &E::ScalarField,
              y: &E::ScalarField,
    ) -> bool {
        // first condition
        let pairing_rhs = E::multi_pairing(proof.vec_D.clone(), &self.srs.vec_V);
        // TODO: optimize further (do one pairing less) by taking its inverse and moving it to the multi_pairing
        let pairing_lhs = E::pairing(&C.C, &self.srs.V_prime);
        // second condition
        let msm_lhs = E::G1::msm_unchecked(&self.srs.vec_H, &proof.f_star_poly.evaluations);
        let l_b = lagrange_x.evaluate(b);
        let msm_rhs = E::G1::msm_unchecked(proof.vec_D.as_slice(), &l_b);
        // third condition
        let y_expected = proof.f_star_poly.evaluate(c);
        // checking all three conditions
        return (pairing_lhs == pairing_rhs) && (msm_lhs == msm_rhs) && (y_expected == *y);
    }
}


#[cfg(test)]
pub mod test {
    use std::cmp::min;

    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::UniformRand;
    use rand::thread_rng;
    use crate::constant_for_curves::{E, ScalarField};
    use super::*;

    #[test]
    fn test_setup() {
        let degree_y = 4usize;
        let degree_x = 4usize;
        let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut thread_rng());
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
        let degree_x = 4usize;
        let degree_y = 16usize;
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        let srs: SRS<E> = PolyCommit::setup(degree_x, degree_y, &mut thread_rng());
        // define the polynomial commitment
        let poly_commit = PolyCommit { srs };
        // random bivariate polynomial
        let polynomial = BivariatePolynomial::random(&mut thread_rng(), domain_x, domain_y, degree_x, degree_y);
        // random points and evaluation
        let b = ScalarField::rand(&mut thread_rng());
        let c = ScalarField::rand(&mut thread_rng());
        let y = polynomial.evaluate(&b, &c);
        // commit to the polynomial
        let com = poly_commit.commit(&polynomial);
        // open the commitment
        let open = poly_commit.open(&polynomial, com.clone(), &b);
        // verify the proof
        let verify = poly_commit.verify(LagrangeBasis { domain: domain_x }, &com, &open, &b, &c, &y);
        assert!(verify);
    }
}

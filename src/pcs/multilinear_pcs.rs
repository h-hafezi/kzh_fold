use std::iter::Sum;
use std::ops::Mul;

use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_std::UniformRand;
use rand::RngCore;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use crate::polynomial::multilinear_polynomial::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_polynomial::math::Math;
use crate::polynomial::multilinear_polynomial::multilinear_poly::MultilinearPolynomial;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SRS<E: Pairing> {
    pub degree_x: usize,
    pub degree_y: usize,
    pub matrix_H: Vec<Vec<E::G1Affine>>,
    pub vec_H: Vec<E::G1Affine>,
    pub vec_V: Vec<E::G2>,
    pub V_prime: E::G2,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<E: Pairing> {
    pub C: E::G1Affine,
    pub aux: Vec<E::G1>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpeningProof<E: Pairing> {
    pub vec_D: Vec<E::G1Affine>,
    pub f_star_poly: MultilinearPolynomial<E::ScalarField, E>,
}

// Define the new struct that encapsulates the functionality of polynomial commitment
pub struct PolyCommit<E: Pairing> {
    pub srs: SRS<E>,
}

pub trait PolyCommitTrait<E: Pairing> {
    fn setup<T: RngCore>(n: usize, m: usize, rng: &mut T) -> SRS<E>;

    fn commit(&self, poly: &MultilinearPolynomial<E::ScalarField, E>) -> Commitment<E>;

    fn open(&self,
            poly: &MultilinearPolynomial<E::ScalarField, E>,
            com: Commitment<E>,
            x: &Vec<E::ScalarField>,
    ) -> OpeningProof<E>;

    fn verify(&self,
              C: &Commitment<E>,
              proof: &OpeningProof<E>,
              x: &Vec<E::ScalarField>,
              y: &Vec<E::ScalarField>,
              z: &E::ScalarField,
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

    fn commit(&self, poly: &MultilinearPolynomial<E::ScalarField, E>) -> Commitment<E> {
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

    fn open(&self, poly: &MultilinearPolynomial<E::ScalarField, E>, com: Commitment<E>, x: &Vec<E::ScalarField>) -> OpeningProof<E> {
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

    fn verify(&self,
              C: &Commitment<E>,
              proof: &OpeningProof<E>,
              x: &Vec<E::ScalarField>,
              y: &Vec<E::ScalarField>,
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
        let msm_rhs = E::G1::msm_unchecked(proof.vec_D.as_slice(), &EqPolynomial::new(x.clone()).evals());

        // third condition
        let y_expected = proof.f_star_poly.evaluate(y);

        // checking all three conditions
        return (pairing_lhs == pairing_rhs) && (msm_lhs == msm_rhs) && (y_expected == *z);
    }
}

#[cfg(test)]
pub mod test {
    use std::cmp::min;

    use ark_ec::pairing::Pairing;
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::{E, ScalarField};
    use crate::pcs::multilinear_pcs::{PolyCommit, PolyCommitTrait, SRS};
    use crate::polynomial::multilinear_polynomial::multilinear_poly::MultilinearPolynomial;

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

        // define the polynomial commitment
        let poly_commit: PolyCommit<E> = PolyCommit { srs };

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
}


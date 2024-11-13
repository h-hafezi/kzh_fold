use crate::math::Math;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use ark_ec::{VariableBaseMSM};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::Rng;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use std::ops::Mul;
use ark_ff::Field;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH3SRS<E: Pairing> {
    /// degree_x = 2 ^ length of x variable
    pub degree_x: usize,
    /// degree_y = 2 ^ length of y variable
    pub degree_y: usize,
    /// degree_z = 2 ^ length of z variable
    pub degree_z: usize,

    pub H_xyz: Vec<E::G1Affine>,
    pub H_yz: Vec<E::G1Affine>,
    pub H_z: Vec<E::G1Affine>,
    pub V_x: Vec<E::G2Affine>,
    pub V_y: Vec<E::G2Affine>,
    pub V_z: Vec<E::G2Affine>,
    pub v: E::G2Affine,
}

pub struct KZH3Opening<E: Pairing> {
    D_y: Vec<E::G1>,
    f_star: MultilinearPolynomial<E::ScalarField>,
}

pub struct KZH3Commitment<E: Pairing> {
    D_x: Vec<E::G1>,
    c: E::G1Affine,
}


impl<E: Pairing> KZH3SRS<E> {
    pub fn setup<R: Rng>(degree_x: usize,
                         degree_y: usize,
                         degree_z: usize,
                         rng: &mut R,
    ) -> KZH3SRS<E> {
        let g = E::G1Affine::rand(rng);
        let v = E::G2Affine::rand(rng);

        let tau_x = {
            let mut elements = Vec::new();
            for _ in 0..degree_x {
                elements.push(E::ScalarField::ONE);
            }
            elements
        };
        let tau_y = {
            let mut elements = Vec::new();
            for _ in 0..degree_y {
                elements.push(E::ScalarField::ONE);
            }
            elements
        };
        let tau_z = {
            let mut elements = Vec::new();
            for _ in 0..degree_z {
                elements.push(E::ScalarField::ONE);
            }
            elements
        };

        let H_xyz = {
            let mut res = Vec::new();
            // i = i_x + degree_x * i_y + (degree_x + degree_y) * i_z
            for i in 0..degree_x * degree_y * degree_z {
                let (i_x, i_y, i_z) = decompose_index(i, degree_y, degree_z);
                let h_xyz = g.mul(tau_x[i_x] * tau_y[i_y] * tau_z[i_z]).into();
                res.push(h_xyz);
            }
            res
        };

        let H_yz = {
            let mut res = Vec::new();
            // i = i_y + degree_y * i_z
            for i in 0..degree_y * degree_z {
                let (i_y, i_z) = (i / degree_z, i % degree_z);
                let h_yz = g.mul(tau_y[i_y] * tau_z[i_z]).into();
                res.push(h_yz);
            }
            res
        };

        let H_z = {
            let mut res = Vec::new();
            for i in 0..degree_z {
                let h_z = g.mul(tau_z[i]).into();
                res.push(h_z);
            }
            res
        };

        for i_y in 0..degree_y {
            for i_z in 0..degree_z {
                let i = i_z + i_y * (degree_z);
                assert_eq!(H_yz[i], H_z[i_z].mul(tau_y[i_y]).into());
            }
        }

        let V_x = {
            let mut res = Vec::new();
            for i in 0..degree_x {
                let v_x = v.mul(tau_x[i]).into();
                res.push(v_x);
            }
            res
        };

        let V_y = {
            let mut res = Vec::new();
            for i in 0..degree_y {
                let v_y = v.mul(tau_y[i]).into();
                res.push(v_y);
            }
            res
        };

        let V_z = {
            let mut res = Vec::new();
            for i in 0..degree_z {
                let v_z = v.mul(tau_z[i]).into();
                res.push(v_z);
            }
            res
        };

        KZH3SRS {
            degree_x,
            degree_y,
            degree_z,
            H_xyz,
            H_yz,
            H_z,
            V_x,
            V_y,
            V_z,
            v,
        }
    }
}


pub fn commit<E: Pairing>(srs: &KZH3SRS<E>, polynomial: &MultilinearPolynomial<E::ScalarField>) -> KZH3Commitment<E> {
    assert_eq!(srs.degree_x * srs.degree_y * srs.degree_z, polynomial.len);

    KZH3Commitment {
        D_x: (0..srs.degree_x)
            .into_par_iter()
            .map(|i| {
                E::G1::msm_unchecked(
                    srs.H_yz.as_slice(),
                    polynomial.get_partial_evaluation_for_boolean_input(i, srs.degree_y * srs.degree_z).as_slice(),
                )
            })
            .collect::<Vec<_>>(),
        c: E::G1::msm(&srs.H_xyz, &polynomial.evaluation_over_boolean_hypercube).unwrap().into(),
    }
}

pub fn open<E: Pairing>(srs: &KZH3SRS<E>,
                        polynomial: &MultilinearPolynomial<E::ScalarField>,
                        x: &[E::ScalarField],
                        y: &[E::ScalarField],
) -> KZH3Opening<E> {
    // Compute f'=f(x,Y,Z); partial evaluation of f using x as the evaluation point. And then compute
    // D_y terms  by committing to each eval slice of f'(Y,Z)
    let f_prime= polynomial.partial_evaluation(x);

    // this takes O(2/3n) but can be precomputed during commitment too
    let D_y = (0..srs.degree_y)
        .into_par_iter()
        .map(|i| {
            E::G1::msm_unchecked(
                srs.H_z.as_slice(),
                f_prime.get_partial_evaluation_for_boolean_input(i, srs.degree_z).as_slice(),
            )
        })
        .collect::<Vec<_>>();

    assert_eq!(x.len(), srs.degree_x.log_2(), "wrong length");
    assert_eq!(y.len(), srs.degree_y.log_2(), "wrong length");

    // compute the partial evaluation of the polynomial
    let f_star = polynomial.partial_evaluation({
        let mut res = Vec::new();
        res.extend_from_slice(&x);
        res.extend_from_slice(&y);
        res
    }.as_slice());

    KZH3Opening {
        D_y,
        f_star,
    }
}

pub fn verify<E: Pairing>(srs: &KZH3SRS<E>,
                          c: &KZH3Commitment<E>,
                          x: &[E::ScalarField],
                          y: &[E::ScalarField],
                          z: &[E::ScalarField],
                          open: &KZH3Opening<E>,
                          output: E::ScalarField,
) {
    // making sure D_x is well-formatted
    let lhs = E::multi_pairing(&c.D_x, &srs.V_x).0;
    let rhs = E::pairing(c.c, &srs.v).0;

    assert_eq!(lhs, rhs);

    // making sure D_y is well formatted
    let new_c = E::G1::msm(
        {
            let mut res = Vec::new();
            for e in &c.D_x {
                res.push(e.clone().into());
            }
            res
        }.as_slice(),
        EqPolynomial::new(x.to_vec()).evals().as_slice(),
    ).unwrap();

    let lhs = E::multi_pairing(&open.D_y, &srs.V_y).0;
    let rhs = E::pairing(new_c, &srs.v).0;

    assert_eq!(lhs, rhs);

    // making sure f^star is well formatter
    let lhs = E::G1::msm(
        srs.H_z.as_slice(),
        open.f_star.evaluation_over_boolean_hypercube.as_slice(),
    ).unwrap();
    let rhs = E::G1::msm(
        {
            let mut res = Vec::new();
            for e in &open.D_y {
                res.push(e.clone().into());
            }
            res
        }.as_slice(),
        EqPolynomial::new(y.to_vec()).evals().as_slice(),
    ).unwrap();

    assert_eq!(lhs, rhs);

    // making sure the output of f_star and the given output are consistent
    assert_eq!(open.f_star.evaluate(z), output);
}

fn decompose_index(i: usize, degree_y: usize, degree_z: usize) -> (usize, usize, usize) {
    // Compute i_z first, as it is the highest order term
    let i_x = i / (degree_y * degree_z);

    // Compute the remainder after removing the contribution of i_z
    let remainder = i % (degree_y * degree_z);

    // Compute i_y next, as it is the middle order term
    let i_y = remainder / degree_z;

    // Finally, compute i_x as the lowest order term
    let i_z = remainder % degree_z;

    (i_x, i_y, i_z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use rand::thread_rng;

    #[test]
    fn test_decompose_index() {
        let i = 1234; // Example index
        let degree_z = 8;  // Must be a power of 2
        let degree_y = 16; // Must be a power of 2

        let (i_x, i_y, i_z) = decompose_index(i, degree_y, degree_z);

        // Check that recomposing gives the original index
        let recomposed_i = i_z + degree_z * i_y + degree_z * degree_y * i_x;
        assert_eq!(i, recomposed_i, "Decomposition and recomposition did not match");
    }

    #[test]
    fn pcs_test() {
        let (degree_x, degree_y, degree_z) = (4usize, 8usize, 16usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        let x :Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..degree_x.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        let y :Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..degree_y.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        let z :Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..degree_z.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        // build the srs
        let srs: KZH3SRS<E> = KZH3SRS::setup(degree_x, degree_y, degree_z, &mut thread_rng());

        // build a random polynomials
        let polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        // evaluate polynomial
        let eval = polynomial.evaluate({
            let mut res = Vec::new();
            res.extend_from_slice(&x);
            res.extend_from_slice(&y);
            res.extend_from_slice(&z);
            res
        }.as_slice());

        // commit to the polynomial
        let c = commit(&srs, &polynomial);

        // open it
        let open = open(&srs, &polynomial, x.as_slice(), y.as_slice());

        // verify the commit
        verify(&srs, &c, x.as_slice(), y.as_slice(), z.as_slice(), &open, eval);
    }
}

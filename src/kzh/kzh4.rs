use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::Rng;
use std::ops::Mul;
use crate::math::Math;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH4SRS<E: Pairing> {
    /// degree_x = 2 ^ length of x variable
    pub degree_x: usize,
    /// degree_y = 2 ^ length of y variable
    pub degree_y: usize,
    /// degree_z = 2 ^ length of z variable
    pub degree_z: usize,
    /// degree_t = 2 ^ length of t variable
    pub degree_t: usize,

    pub H_xyzt: Vec<E::G1Affine>,
    pub H_yzt: Vec<E::G1Affine>,
    pub H_zt: Vec<E::G1Affine>,
    pub H_t: Vec<E::G1Affine>,

    pub V_x: Vec<E::G2Affine>,
    pub V_y: Vec<E::G2Affine>,
    pub V_z: Vec<E::G2Affine>,
    pub V_t: Vec<E::G2Affine>,

    pub v: E::G2Affine,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct KZH4Opening<E: Pairing> {
    D_x: Vec<E::G1>,
    D_y: Vec<E::G1>,
    D_z: Vec<E::G1>,
    f_star: MultilinearPolynomial<E::ScalarField>,
}

pub struct KZH4Commitment<E: Pairing> {
    c: E::G1Affine,
}

impl<E: Pairing> KZH4SRS<E> {
    pub fn setup<R: Rng>(degree_x: usize,
                         degree_y: usize,
                         degree_z: usize,
                         degree_t: usize,
                         rng: &mut R,
    ) -> KZH4SRS<E> {
        let g = E::G1Affine::rand(rng);
        let v = E::G2Affine::rand(rng);

        let tau_x = {
            let mut elements = Vec::new();
            for _ in 0..degree_x {
                elements.push(E::ScalarField::rand(rng));
            }
            elements
        };
        let tau_y = {
            let mut elements = Vec::new();
            for _ in 0..degree_y {
                elements.push(E::ScalarField::rand(rng));
            }
            elements
        };
        let tau_z = {
            let mut elements = Vec::new();
            for _ in 0..degree_z {
                elements.push(E::ScalarField::rand(rng));
            }
            elements
        };
        let tau_t = {
            let mut elements = Vec::new();
            for _ in 0..degree_t {
                elements.push(E::ScalarField::rand(rng));
            }
            elements
        };

        let H_xyzt = {
            let mut res = Vec::new();
            for i in 0..degree_x * degree_y * degree_z * degree_t {
                let (i_x, i_y, i_z, i_t) = decompose_index(i, degree_y, degree_z, degree_t);
                let h_xyzt = g.mul(tau_x[i_x] * tau_y[i_y] * tau_z[i_z] * tau_t[i_t]).into();
                res.push(h_xyzt);
            }
            res
        };

        let H_yzt = {
            let mut res = Vec::new();
            for i in 0..degree_y * degree_z * degree_t {
                let (i_y, i_z, i_t) = {
                    let i_y = i / (degree_z * degree_t);
                    let remainder = i % (degree_z * degree_t);
                    let i_z = remainder / degree_t;
                    let i_t = remainder % degree_t;

                    (i_y, i_z, i_t)
                };
                let h_yzt = g.mul(tau_y[i_y] * tau_z[i_z] * tau_t[i_t]).into();
                res.push(h_yzt);
            }
            res
        };

        let H_zt = {
            let mut res = Vec::new();
            for i in 0..degree_z * degree_t {
                let (i_z, i_t) = {
                    let i_z = i / degree_t;
                    let i_t = i % degree_t;

                    (i_z, i_t)
                };
                let h_zt = g.mul(tau_z[i_z] * tau_t[i_t]).into();
                res.push(h_zt);
            }
            res
        };

        let H_t = {
            let mut res = Vec::new();
            for i in 0..degree_t {
                let h_t = g.mul(tau_t[i]).into();
                res.push(h_t);
            }
            res
        };


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

        let V_t = {
            let mut res = Vec::new();
            for i in 0..degree_t {
                let v_t = v.mul(tau_t[i]).into();
                res.push(v_t);
            }
            res
        };


        KZH4SRS {
            degree_x,
            degree_y,
            degree_z,
            degree_t,
            H_xyzt,
            H_yzt,
            H_zt,
            H_t,
            V_x,
            V_y,
            V_z,
            V_t,
            v,
        }
    }
}

pub fn commit<E: Pairing>(srs: &KZH4SRS<E>, polynomial: &MultilinearPolynomial<E::ScalarField>) -> KZH4Commitment<E> {
    assert_eq!(srs.degree_x * srs.degree_y * srs.degree_z * srs.degree_t, polynomial.len);

    KZH4Commitment {
        c: E::G1::msm(&srs.H_xyzt, &polynomial.evaluation_over_boolean_hypercube).unwrap().into(),
    }
}

pub fn open<E: Pairing>(srs: &KZH4SRS<E>,
                        polynomial: &MultilinearPolynomial<E::ScalarField>,
                        x: &[E::ScalarField],
                        y: &[E::ScalarField],
                        z: &[E::ScalarField],
) -> KZH4Opening<E> {
    let D_x = (0..srs.degree_x)
        .into_iter()
        .map(|i| {
            E::G1::msm_unchecked(
                srs.H_yzt.as_slice(),
                polynomial.get_partial_evaluation_for_boolean_input(i, srs.degree_y * srs.degree_z * srs.degree_t).as_slice(),
            )
        })
        .collect::<Vec<_>>();

    let D_y = (0..srs.degree_y)
        .into_iter()
        .map(|i| {
            E::G1::msm_unchecked(
                srs.H_zt.as_slice(),
                polynomial.partial_evaluation(x).get_partial_evaluation_for_boolean_input(i, srs.degree_z * srs.degree_t).as_slice(),
            )
        })
        .collect::<Vec<_>>();

    let D_z = (0..srs.degree_z)
        .into_iter()
        .map(|i| {
            E::G1::msm_unchecked(
                srs.H_t.as_slice(),
                polynomial.partial_evaluation(
                    {
                        let mut res = Vec::new();
                        res.extend_from_slice(x);
                        res.extend_from_slice(y);
                        res
                    }.as_slice()
                ).get_partial_evaluation_for_boolean_input(i, srs.degree_t).as_slice(),
            )
        })
        .collect::<Vec<_>>();


    assert_eq!(x.len(), srs.degree_x.log_2(), "wrong length");
    assert_eq!(y.len(), srs.degree_y.log_2(), "wrong length");
    assert_eq!(z.len(), srs.degree_z.log_2(), "wrong length");

    // compute the partial evaluation of the polynomial
    let f_star = polynomial.partial_evaluation({
        let mut res = Vec::new();
        res.extend_from_slice(&x);
        res.extend_from_slice(&y);
        res.extend_from_slice(&z);
        res
    }.as_slice());

    KZH4Opening {
        D_x,
        D_y,
        D_z,
        f_star,
    }
}


pub fn verify<E: Pairing>(srs: &KZH4SRS<E>,
                          c: &KZH4Commitment<E>,
                          x: &[E::ScalarField],
                          y: &[E::ScalarField],
                          z: &[E::ScalarField],
                          t: &[E::ScalarField],
                          open: &KZH4Opening<E>,
                          output: E::ScalarField,
) {
    // making sure D_x is well-formatted
    let lhs = E::multi_pairing(&open.D_x, &srs.V_x).0;
    let rhs = E::pairing(c.c, &srs.v).0;

    assert_eq!(lhs, rhs);

    // making sure D_y is well formatted
    let new_c = E::G1::msm(
        {
            let mut res = Vec::new();
            for e in &open.D_x {
                res.push(e.clone().into());
            }
            res
        }.as_slice(),
        EqPolynomial::new(x.to_vec()).evals().as_slice(),
    ).unwrap();

    let lhs = E::multi_pairing(&open.D_y, &srs.V_y).0;
    let rhs = E::pairing(new_c, &srs.v).0;

    assert_eq!(lhs, rhs);

    // making sure D_z is well formatted
    let new_c = E::G1::msm(
        {
            let mut res = Vec::new();
            for e in &open.D_y {
                res.push(e.clone().into());
            }
            res
        }.as_slice(),
        EqPolynomial::new(y.to_vec()).evals().as_slice(),
    ).unwrap();

    let lhs = E::multi_pairing(&open.D_z, &srs.V_z).0;
    let rhs = E::pairing(new_c, &srs.v).0;

    assert_eq!(lhs, rhs);

    // making sure f^star is well formatter
    let lhs = E::G1::msm(
        srs.H_t.as_slice(),
        open.f_star.evaluation_over_boolean_hypercube.as_slice(),
    ).unwrap();
    let rhs = E::G1::msm(
        {
            let mut res = Vec::new();
            for e in &open.D_z {
                res.push(e.clone().into());
            }
            res
        }.as_slice(),
        EqPolynomial::new(z.to_vec()).evals().as_slice(),
    ).unwrap();

    assert_eq!(lhs, rhs);

    // making sure the output of f_star and the given output are consistent
    assert_eq!(open.f_star.evaluate(t), output);
}


fn decompose_index(i: usize, degree_y: usize, degree_z: usize, degree_t: usize) -> (usize, usize, usize, usize) {
    // Compute i_z first, as it is the highest order term
    let i_x = i / (degree_y * degree_z * degree_t);

    // Compute the remainder after removing the contribution of i_z
    let remainder = i % (degree_y * degree_z * degree_t);

    // Compute i_y next, as it is the middle order term
    let i_y = remainder / (degree_z * degree_t);

    // Finally, compute i_x as the lowest order term
    let remainder = remainder % (degree_z * degree_t);

    let i_z = remainder / degree_t;

    let i_t = remainder % degree_t;

    (i_x, i_y, i_z, i_t)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use rand::thread_rng;

    #[test]
    fn pcs_test() {
        let (degree_x, degree_y, degree_z, degree_t) = (4usize, 8usize, 16usize, 8usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2() +degree_t.log_2();

        let x: Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..degree_x.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        let y: Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..degree_y.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        let z: Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..degree_z.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        let t: Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..degree_t.log_2() {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };

        // build the srs
        let srs: KZH4SRS<E> = KZH4SRS::setup(degree_x, degree_y, degree_z, degree_t, &mut thread_rng());

        // build a random polynomials
        let polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        // evaluate polynomial
        let eval = polynomial.evaluate({
            let mut res = Vec::new();
            res.extend_from_slice(&x);
            res.extend_from_slice(&y);
            res.extend_from_slice(&z);
            res.extend_from_slice(&t);
            res
        }.as_slice());

        // commit to the polynomial
        let c = commit(&srs, &polynomial);

        // open it
        let open = open(&srs, &polynomial, x.as_slice(), y.as_slice(), z.as_slice());

        // verify the commit
        verify(&srs, &c, x.as_slice(), y.as_slice(), z.as_slice(), t.as_slice(), &open, eval);
    }
}

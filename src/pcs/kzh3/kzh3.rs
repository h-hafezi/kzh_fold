use crate::math::Math;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
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
    D_x: Vec<E::G1>,
    D_xy: Vec<E::G1>,
    f_star: MultilinearPolynomial<E::ScalarField>,
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
                let (i_x, i_y, i_z) = decompose_index(i, degree_x, degree_y);
                let h_xyz = g.mul(tau_x[i_x] * tau_y[i_y] * tau_z[i_z]).into();
                res.push(h_xyz);
            }
            res
        };

        let H_yz = {
            let mut res = Vec::new();
            // i = i_y + degree_y * i_z
            for i in 0..degree_y * degree_z {
                let (i_y, i_z) = (i % degree_y, i / degree_y);
                let h_yz = g.mul(tau_y[i_y] * tau_z[i_z]).into();
                res.push(h_yz);
            }
            res
        };

        for i in 0..degree_x * degree_y * degree_z {
            let (i_x, i_y, i_z) = decompose_index(i, degree_x, degree_y);
            let kir = i_x + i_y * degree_x + i_z * (degree_x * degree_y);
            assert_eq!(kir, i);
            let i_yz = i_y + degree_y * i_z;
            assert_eq!(H_xyz[i], H_yz[i_yz].mul(tau_x[i_x]).into());
        }


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
                let i = i_y + i_z * (degree_y);
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

        for i in 0..degree_x * degree_y * degree_z {
            let (i_x, i_y, i_z) = decompose_index(i, degree_x, degree_y);
            let i_yz = i_y + degree_y * i_z;
            assert_eq!(E::pairing(H_xyz[i], v).0, E::pairing(H_yz[i_yz], V_x[i_x]).0);
        }

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


pub fn commit<E: Pairing>(srs: &KZH3SRS<E>, polynomial: &MultilinearPolynomial<E::ScalarField>) -> E::G1Affine {
    assert_eq!(srs.degree_x * srs.degree_y * srs.degree_z, polynomial.len);

    E::G1::msm(&srs.H_xyz, &polynomial.evaluation_over_boolean_hypercube).unwrap().into()
}

pub fn open<E: Pairing>(srs: &KZH3SRS<E>,
                        polynomial: &MultilinearPolynomial<E::ScalarField>,
                        x: &[E::ScalarField],
                        y: &[E::ScalarField],
) -> KZH3Opening<E> {
    let D_x = (0..srs.degree_x)
        .into_par_iter()
        .map(|i| {
            E::G1::msm_unchecked(
                srs.H_yz.as_slice(),
                polynomial.get_partial_evaluation_for_boolean_input(i, srs.degree_y * srs.degree_z).as_slice(),
            )
        })
        .collect::<Vec<_>>();

    let D_xy = (0..srs.degree_x * srs.degree_y)
        .into_par_iter()
        .map(|i| {
            E::G1::msm_unchecked(
                srs.H_z.as_slice(),
                polynomial.get_partial_evaluation_for_boolean_input(i, srs.degree_z).as_slice(),
            )
        })
        .collect::<Vec<_>>();

    assert_eq!(x.len(), srs.degree_x.log_2(), "wrong length");
    assert_eq!(y.len(), srs.degree_y.log_2(), "wrong length");

    let concat = {
        let mut res = Vec::new();
        res.extend_from_slice(x);
        res.extend_from_slice(y);

        res
    };

    let f_star = polynomial.partial_evaluation(concat.as_slice());

    KZH3Opening {
        D_x,
        D_xy,
        f_star,
    }
}

pub fn verify<E: Pairing>(srs: &KZH3SRS<E>,
                          c: E::G1Affine,
                          x: &[E::ScalarField],
                          y: &[E::ScalarField],
                          z: &[E::ScalarField],
                          open: &KZH3Opening<E>,
                          output: E::ScalarField,
) {
    let lhs = E::multi_pairing(&open.D_x, &srs.V_x).0;
    let rhs = E::pairing(c, &srs.v).0;

    assert_eq!(lhs, rhs);

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
    //let lhs = E::multi_pairing(&open.D_xy, &srs.V_y).0;
    //let rhs = E::pairing(new_c, &srs.v).0;
    //assert_eq!(lhs, rhs);

    //let lhs = E::multi_pairing(&open.D_xy, &srs.)
}

fn decompose_index(i: usize, degree_x: usize, degree_y: usize) -> (usize, usize, usize) {
    // Compute i_z first, as it is the highest order term
    let i_z = i / (degree_x * degree_y);

    // Compute the remainder after removing the contribution of i_z
    let remainder = i % (degree_x * degree_y);

    // Compute i_y next, as it is the middle order term
    let i_y = remainder / degree_x;

    // Finally, compute i_x as the lowest order term
    let i_x = remainder % degree_x;

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
        let degree_x = 8;  // Must be a power of 2
        let degree_y = 16; // Must be a power of 2

        let (i_x, i_y, i_z) = decompose_index(i, degree_x, degree_y);

        // Output the decomposed values for verification
        println!("i_x: {}, i_y: {}, i_z: {}", i_x, i_y, i_z);

        // Check that recomposing gives the original index
        let recomposed_i = i_x + degree_x * i_y + degree_x * degree_y * i_z;
        assert_eq!(i, recomposed_i, "Decomposition and recomposition did not match");
    }

    #[test]
    fn pcs_test() {
        let (degree_x, degree_y, degree_z) = (16usize, 8usize, 2usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();
        let x = vec![ScalarField::rand(&mut thread_rng()); degree_x.log_2()];
        let y = vec![ScalarField::rand(&mut thread_rng()); degree_y.log_2()];
        let z = vec![ScalarField::rand(&mut thread_rng()); degree_z.log_2()];
        let srs: KZH3SRS<E> = KZH3SRS::setup(degree_x, degree_y, degree_z, &mut thread_rng());
        let polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        let c = commit(&srs, &polynomial);

        let open = open(&srs, &polynomial, x.as_slice(), y.as_slice());

        verify(&srs, c, x.as_slice(), y.as_slice(), z.as_slice(), &open, ScalarField::rand(&mut thread_rng()));
    }
}

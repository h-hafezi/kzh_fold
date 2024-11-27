use crate::kzh::KZH;
use crate::math::Math;
use crate::nexus_spartan::commitment_traits::ToAffine;
use crate::polynomial::eq_poly::eq_poly::EqPolynomial;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, VariableBaseMSM};
use ark_ff::{AdditiveGroup, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::Rng;
use std::marker::PhantomData;
use std::ops::Mul;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH4<E: Pairing> {
    phantom: PhantomData<E>,
}

impl<E: Pairing> KZH<E> for KZH4<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    type Degree = (usize, usize, usize, usize);
    type SRS = KZH4SRS<E>;
    type Commitment = KZH4Commitment<E>;
    type Opening = KZH4Opening<E>;

    fn split_input(srs: &Self::SRS, input: &[E::ScalarField]) -> Vec<Vec<E::ScalarField>> {
        let total_length = srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2() + srs.degree_t.log_2();

        // If r is smaller than the required length, extend it with zeros at the beginning
        let mut extended_r = input.to_vec();
        if input.len() < total_length {
            let mut zeros = vec![E::ScalarField::ZERO; total_length - input.len()];
            zeros.extend(extended_r);  // Prepend zeros to the beginning
            extended_r = zeros;
        }

        // Split the vector into two parts
        let r_x = extended_r[..srs.degree_x.log_2()].to_vec();
        let r_y = extended_r[srs.degree_x.log_2()..srs.degree_x.log_2() + srs.degree_y.log_2()].to_vec();
        let r_z = extended_r[srs.degree_x.log_2() + srs.degree_y.log_2()..srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2()].to_vec();
        let r_t = extended_r[srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2()..].to_vec();

        vec![r_x, r_y, r_z, r_t]
    }

    fn get_degree_from_maximum_supported_degree(n: usize) -> Self::Degree {
        (n / 4 + 1, n / 4 + 1, n / 4 + 1, n / 4 + 1)
    }

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS {
        let (degree_x, degree_y, degree_z, degree_t) = Self::get_degree_from_maximum_supported_degree(maximum_degree);

        let degree_x = 1 << degree_x;
        let degree_y = 1 << degree_y;
        let degree_z = 1 << degree_z;
        let degree_t = 1 << degree_t;

        let g = E::G1Affine::rand(rng);
        let v = E::G2Affine::rand(rng);

        let tau_x: Vec<E::ScalarField> = (0..degree_x).map(|_| E::ScalarField::rand(rng)).collect();
        let tau_y: Vec<E::ScalarField> = (0..degree_y).map(|_| E::ScalarField::rand(rng)).collect();
        let tau_z: Vec<E::ScalarField> = (0..degree_z).map(|_| E::ScalarField::rand(rng)).collect();
        let tau_t: Vec<E::ScalarField> = (0..degree_t).map(|_| E::ScalarField::rand(rng)).collect();

        let H_xyzt: Vec<_> = (0..degree_x * degree_y * degree_z * degree_t)
            .map(|i| {
                let (i_x, i_y, i_z, i_t) = decompose_index(i, degree_y, degree_z, degree_t);
                g.mul(tau_x[i_x] * tau_y[i_y] * tau_z[i_z] * tau_t[i_t]).into()
            }).collect();

        let H_yzt: Vec<_> = (0..degree_y * degree_z * degree_t)
            .map(|i| {
                let i_y = i / (degree_z * degree_t);
                let remainder = i % (degree_z * degree_t);
                let i_z = remainder / degree_t;
                let i_t = remainder % degree_t;

                g.mul(tau_y[i_y] * tau_z[i_z] * tau_t[i_t]).into()
            }).collect();

        let H_zt: Vec<_> = (0..degree_z * degree_t)
            .map(|i| {
                let i_z = i / degree_t;
                let i_t = i % degree_t;
                g.mul(tau_z[i_z] * tau_t[i_t]).into()
            })
            .collect();

        let H_t: Vec<_> = (0..degree_t)
            .map(|i| g.mul(tau_t[i]).into())
            .collect();

        let V_x: Vec<_> = (0..degree_x).map(|i| v.mul(tau_x[i]).into()).collect();
        let V_y: Vec<_> = (0..degree_y).map(|i| v.mul(tau_y[i]).into()).collect();
        let V_z: Vec<_> = (0..degree_z).map(|i| v.mul(tau_z[i]).into()).collect();
        let V_t: Vec<_> = (0..degree_t).map(|i| v.mul(tau_t[i]).into()).collect();

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

    fn commit(srs: &Self::SRS, poly: &MultilinearPolynomial<E::ScalarField>) -> Self::Commitment {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2()+srs.degree_t.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        assert_eq!(srs.degree_x * srs.degree_y * srs.degree_z * srs.degree_t, poly.len);

        KZH4Commitment {
            C: E::G1::msm(&srs.H_xyzt, &poly.evaluation_over_boolean_hypercube).unwrap().into(),
        }
    }

    fn open(srs: &Self::SRS, input: &[E::ScalarField], _com: &Self::Commitment, poly: &MultilinearPolynomial<E::ScalarField>) -> Self::Opening {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2()+srs.degree_t.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        let split_input = Self::split_input(&srs, input);

        let D_x = (0..srs.degree_x)
            .into_iter()
            .map(|i| {
                E::G1::msm_unchecked(
                    srs.H_yzt.as_slice(),
                    poly.get_partial_evaluation_for_boolean_input(i, srs.degree_y * srs.degree_z * srs.degree_t).as_slice(),
                )
            })
            .collect::<Vec<_>>();

        let D_y = (0..srs.degree_y)
            .into_iter()
            .map(|i| {
                E::G1::msm_unchecked(
                    srs.H_zt.as_slice(),
                    poly.partial_evaluation(split_input[0].as_slice()).get_partial_evaluation_for_boolean_input(i, srs.degree_z * srs.degree_t).as_slice(),
                )
            })
            .collect::<Vec<_>>();

        let D_z = (0..srs.degree_z)
            .into_iter()
            .map(|i| {
                E::G1::msm_unchecked(
                    srs.H_t.as_slice(),
                    poly.partial_evaluation(
                        {
                            let mut res = Vec::new();
                            res.extend_from_slice(split_input[0].as_slice());
                            res.extend_from_slice(split_input[1].as_slice());
                            res
                        }.as_slice()
                    ).get_partial_evaluation_for_boolean_input(i, srs.degree_t).as_slice(),
                )
            })
            .collect::<Vec<_>>();


        assert_eq!(split_input[0].len(), srs.degree_x.log_2(), "wrong length");
        assert_eq!(split_input[1].len(), srs.degree_y.log_2(), "wrong length");
        assert_eq!(split_input[2].len(), srs.degree_z.log_2(), "wrong length");

        // compute the partial evaluation of the polynomial
        let f_star = poly.partial_evaluation({
            let mut res = Vec::new();
            res.extend_from_slice(split_input[0].as_slice());
            res.extend_from_slice(split_input[1].as_slice());
            res.extend_from_slice(split_input[2].as_slice());
            res
        }.as_slice());

        KZH4Opening {
            D_x,
            D_y,
            D_z,
            f_star,
        }
    }

    fn verify(srs: &Self::SRS, input: &[E::ScalarField], output: &E::ScalarField, com: &Self::Commitment, open: &Self::Opening) {
        let split_input = Self::split_input(&srs, input);

        // making sure D_x is well-formatted
        let lhs = E::multi_pairing(&open.D_x, &srs.V_x).0;
        let rhs = E::pairing(com.C, &srs.v).0;

        assert_eq!(lhs, rhs);

        // making sure D_y is well formatted
        let new_c = E::G1::msm(
            &open.D_x.iter().map(|e| e.clone().into()).collect::<Vec<_>>().as_slice(),
            EqPolynomial::new(split_input[0].clone()).evals().as_slice(),
        ).unwrap();

        let lhs = E::multi_pairing(&open.D_y, &srs.V_y).0;
        let rhs = E::pairing(new_c, &srs.v).0;

        assert_eq!(lhs, rhs);

        // making sure D_z is well formatted
        let new_c = E::G1::msm(
            &open.D_y.iter().map(|e| e.clone().into()).collect::<Vec<_>>().as_slice(),
            EqPolynomial::new(split_input[1].clone()).evals().as_slice(),
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
            &open.D_z.iter().map(|e| e.clone().into()).collect::<Vec<_>>().as_slice(),
            EqPolynomial::new(split_input[2].clone()).evals().as_slice(),
        ).unwrap();

        assert_eq!(lhs, rhs);

        // making sure the output of f_star and the given output are consistent
        assert_eq!(open.f_star.evaluate(split_input[3].as_slice()), *output);
    }
}

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

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH4Opening<E: Pairing> {
    pub D_x: Vec<E::G1>,
    pub D_y: Vec<E::G1>,
    pub D_z: Vec<E::G1>,
    pub f_star: MultilinearPolynomial<E::ScalarField>,
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
pub struct KZH4Commitment<E: Pairing> {
    pub C: E::G1Affine,
}

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for KZH4Commitment<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>) {
        Transcript::append_point::<E>(transcript, label, &self.C);
    }
}

impl<E: Pairing> ToAffine<E> for KZH4Commitment<E> {
    fn to_affine(self) -> E::G1Affine {
        self.C
    }
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
    use crate::constant_for_curves::{ScalarField as F, E};
    use rand::thread_rng;

    #[test]
    fn pcs_test() {
        let (degree_x, degree_y, degree_z, degree_t) = (4usize, 8usize, 16usize, 8usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2() + degree_t.log_2();

        let input: Vec<F> = {
            let mut res = Vec::new();
            for _ in 0..num_vars {
                res.push(F::rand(&mut thread_rng()));
            }
            res
        };

        // build the srs
        let srs: KZH4SRS<E> = KZH4::setup(num_vars, &mut thread_rng());

        // build a random polynomials
        let polynomial: MultilinearPolynomial<F> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        // evaluate polynomial
        let eval = polynomial.evaluate(input.as_slice());

        // commit to the polynomial
        let com = KZH4::commit(&srs, &polynomial);

        // open it
        let open = KZH4::open(&srs, input.as_slice(), &com, &polynomial);

        // verify the commit
        KZH4::verify(&srs, input.as_slice(), &eval, &com, &open);
    }
}

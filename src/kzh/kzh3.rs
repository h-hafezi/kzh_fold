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
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use std::marker::PhantomData;
use std::ops::Mul;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH3<E: Pairing> {
    phantom: PhantomData<E>,
}

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

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct KZH3Opening<E: Pairing> {
    D_y: Vec<E::G1>,
    f_star: MultilinearPolynomial<E::ScalarField>,
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
pub struct KZH3Commitment<E: Pairing> {
    D_x: Vec<E::G1>,
    C: E::G1Affine,
}

impl<E: Pairing> KZH<E> for KZH3<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    type Degree = (usize, usize, usize);
    type SRS = KZH3SRS<E>;
    type Commitment = KZH3Commitment<E>;
    type Opening = KZH3Opening<E>;

    fn split_input(srs: &Self::SRS, input: &[E::ScalarField]) -> Vec<Vec<E::ScalarField>> {
        let total_length = srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2();

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
        let r_z = extended_r[srs.degree_x.log_2() + srs.degree_y.log_2()..].to_vec();

        vec![r_x, r_y, r_z]
    }

    fn get_degree_from_maximum_supported_degree(n: usize) -> (usize, usize, usize) {
        match n % 3 {
            0 => (n / 3, n / 3, n / 3),
            1 => (n / 3 + 1, n / 3, n / 3),
            2 => (n / 3 + 1, n / 3 + 1, n / 3),
            _ => unreachable!(),
        }
    }

    fn setup<R: Rng>(maximum_degree: usize, rng: &mut R) -> Self::SRS {
        let (degree_x, degree_y, degree_z) = Self::get_degree_from_maximum_supported_degree(maximum_degree);

        let (degree_x, degree_y, degree_z) = (1 << degree_x, 1 << degree_y, 1 << degree_z);

        let (g, v) = (E::G1Affine::rand(rng), E::G2Affine::rand(rng));

        let tau_x: Vec<E::ScalarField> = (0..degree_x).map(|_| E::ScalarField::rand(rng)).collect();
        let tau_y: Vec<E::ScalarField> = (0..degree_y).map(|_| E::ScalarField::rand(rng)).collect();
        let tau_z: Vec<E::ScalarField> = (0..degree_z).map(|_| E::ScalarField::rand(rng)).collect();


        let H_xyz: Vec<_> = (0..degree_x * degree_y * degree_z)
            .map(|i| {
                let (i_x, i_y, i_z) = decompose_index(i, degree_y, degree_z);
                g.mul(tau_x[i_x] * tau_y[i_y] * tau_z[i_z]).into()
            }).collect();

        let H_yz: Vec<_> = (0..degree_y * degree_z)
            .map(|i| {
                let (i_y, i_z) = (i / degree_z, i % degree_z);
                g.mul(tau_y[i_y] * tau_z[i_z]).into()
            }).collect();

        let H_z: Vec<_> = (0..degree_z)
            .map(|i| g.mul(tau_z[i]).into())
            .collect();

        for i_y in 0..degree_y {
            for i_z in 0..degree_z {
                let i = i_z + i_y * (degree_z);
                assert_eq!(H_yz[i], H_z[i_z].mul(tau_y[i_y]).into());
            }
        }

        let V_x: Vec<_> = (0..degree_x)
            .map(|i| v.mul(tau_x[i]).into())
            .collect();

        let V_y: Vec<_> = (0..degree_y)
            .map(|i| v.mul(tau_y[i]).into())
            .collect();

        let V_z: Vec<_> = (0..degree_z)
            .map(|i| v.mul(tau_z[i]).into())
            .collect();

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

    fn commit(srs: &Self::SRS, poly: &MultilinearPolynomial<E::ScalarField>) -> Self::Commitment {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        KZH3Commitment {
            D_x: (0..srs.degree_x)
                .into_par_iter()
                .map(|i| {
                    E::G1::msm_unchecked(
                        srs.H_yz.as_slice(),
                        poly.get_partial_evaluation_for_boolean_input(i, srs.degree_y * srs.degree_z).as_slice(),
                    )
                })
                .collect::<Vec<_>>(),
            C: E::G1::msm(&srs.H_xyz, &poly.evaluation_over_boolean_hypercube).unwrap().into(),
        }
    }

    fn open(srs: &Self::SRS, input: &[E::ScalarField], _: &Self::Commitment, poly: &MultilinearPolynomial<E::ScalarField>) -> Self::Opening {
        let len = srs.degree_x.log_2() + srs.degree_y.log_2() + srs.degree_z.log_2();
        let poly = poly.extend_number_of_variables(len);
        assert_eq!(poly.num_variables, len);
        assert_eq!(poly.len, 1 << poly.num_variables);
        assert_eq!(poly.evaluation_over_boolean_hypercube.len(), poly.len);

        let split_input = Self::split_input(&srs, input);

        // Compute f'=f(x,Y,Z); partial evaluation of f using x as the evaluation point. And then compute
        // D_y terms  by committing to each eval slice of f'(Y,Z)
        let f_prime = poly.partial_evaluation(split_input[0].as_slice());

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

        assert_eq!(split_input[0].len(), srs.degree_x.log_2(), "wrong length");
        assert_eq!(split_input[1].len(), srs.degree_y.log_2(), "wrong length");

        // compute the partial evaluation of the polynomial
        let f_star = poly.partial_evaluation({
            let mut res = Vec::new();
            res.extend_from_slice(split_input[0].as_slice());
            res.extend_from_slice(split_input[1].as_slice());
            res
        }.as_slice());

        KZH3Opening {
            D_y,
            f_star,
        }
    }

    fn verify(srs: &Self::SRS, input: &[E::ScalarField], output: &E::ScalarField, com: &Self::Commitment, open: &Self::Opening) {
        let split_input = Self::split_input(&srs, input);

        // making sure D_x is well-formatted
        let lhs = E::multi_pairing(&com.D_x, &srs.V_x).0;
        let rhs = E::pairing(com.C, &srs.v).0;

        assert_eq!(lhs, rhs);

        // making sure D_y is well formatted
        let new_c = E::G1::msm(
            &com.D_x.iter().map(|e| e.clone().into()).collect::<Vec<_>>().as_slice(),
            EqPolynomial::new(split_input[0].clone()).evals().as_slice(),
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
            &open.D_y.iter().map(|e| e.clone().into()).collect::<Vec<_>>().as_slice(),
            EqPolynomial::new(split_input[1].clone()).evals().as_slice(),
        ).unwrap();

        assert_eq!(lhs, rhs);

        // making sure the output of f_star and the given output are consistent
        assert_eq!(open.f_star.evaluate(split_input[2].as_slice()), *output);
    }
}

impl<E: Pairing, F: PrimeField + Absorb> AppendToTranscript<F> for KZH3Commitment<E>
where
    E: Pairing<ScalarField=F>,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript<F>) {
        Transcript::append_point::<E>(transcript, label, &self.C);
    }
}

impl<E: Pairing> ToAffine<E> for KZH3Commitment<E> {
    fn to_affine(self) -> E::G1Affine {
        self.C
    }
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
    use crate::constant_for_curves::{ScalarField as F, E};
    use crate::kzh::kzh3::{decompose_index, KZH3, KZH3SRS};
    use crate::kzh::KZH;
    use crate::math::Math;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use ark_std::UniformRand;
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
        let (degree_x, degree_y, degree_z) = (4usize, 4usize, 4usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        let input: Vec<F> = (0..num_vars)
            .map(|_| F::rand(&mut thread_rng()))
            .collect();

        // build the srs
        let srs: KZH3SRS<E> = KZH3::setup((degree_x * degree_y * degree_z).log_2(), &mut thread_rng());

        // build a random polynomials
        let polynomial: MultilinearPolynomial<F> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        // evaluate polynomial
        let eval = polynomial.evaluate(input.as_slice());

        // commit to the polynomial
        let c = KZH3::commit(&srs, &polynomial);

        // open it
        let open = KZH3::open(&srs, input.as_slice(), &c, &polynomial);

        // verify the commit
        KZH3::verify(&srs, input.as_slice(), &eval, &c, &open);
    }
}

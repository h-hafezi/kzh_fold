use crate::gadgets::non_native::util::convert_affine_to_scalars;
use crate::kzh::kzh3::{KZH3Commitment, KZH3Opening, KZH3, KZH3SRS};
use crate::kzh::KZH;
use crate::kzh_fold::eq_tree::EqTree;
use crate::kzh_fold::generate_random_elements;
use crate::kzh_fold::kzh2_fold::Acc2Instance;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::Transcript;
use crate::utils::inner_product;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, VariableBaseMSM};
use ark_ff::{AdditiveGroup, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::RngCore;
use std::ops::{Add, Mul, Neg};

type Acc3Proof<E: Pairing> = (E::G1Affine, E::G1Affine, E::G1Affine, E::G1Affine);


#[derive(Clone, Debug)]
pub struct Acc3SRS<E: Pairing> {
    // vector of size 2 * degree_x - 1
    pub k_x: Vec<E::G1Affine>,
    // vector of size 2 * degree_y - 1
    pub k_y: Vec<E::G1Affine>,
    // vector of size 2 * degree_z - 1
    pub k_z: Vec<E::G1Affine>,

    pub k_prime: E::G1Affine,

    pub pc_srs: KZH3SRS<E>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Acc3Instance<E: Pairing> {
    pub C: E::G1Affine,
    pub C_y: E::G1Affine,
    pub T: E::G1Affine,
    pub E: Acc3Proof<E>,

    // vector of length log2(degree_x)
    pub x: Vec<E::ScalarField>,
    // vector of length log2(degree_y)
    pub y: Vec<E::ScalarField>,
    // vector of length log2(degree_z)
    pub z: Vec<E::ScalarField>,
    // result of poylnomial evaluation
    pub output: E::ScalarField,
}

impl<E: Pairing> Acc3Instance<E>
where
    E::ScalarField: PrimeField,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
{
    /// returns a vector of E::Scalar that can be used to Poseidon hash
    /// the tricky part is the affine points which are transformed via convert_affine_to_scalars
    /// which basically converts each coordinate from base field into scalar field
    pub fn to_sponge_field_elements(&self) -> Vec<E::ScalarField> {
        let mut dest = Vec::new();

        // Use the closure for C, T, and E
        let (c_x, c_y) = convert_affine_to_scalars::<E>(self.C);
        let (t_x, t_y) = convert_affine_to_scalars::<E>(self.T);

        // Extend the destination vector with the computed values
        dest.extend(vec![c_x, c_y, t_x, t_y]);

        for E in vec![&self.E.0, &self.E.1, &self.E.2, &self.E.3] {
            let (e_x, e_y) = convert_affine_to_scalars::<E>(*E);
            // Extend the destination vector with the computed values
            dest.extend(vec![e_x, e_y]);
        }

        // Extend with other scalar fields
        dest.extend(self.x.clone());
        dest.extend(self.y.clone());
        dest.extend(self.z.clone());
        dest.push(self.output.clone());

        dest
    }
}


#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize, Derivative)]
pub struct Acc3Witness<E: Pairing> {
    pub D_x: Vec<E::G1>,
    pub D_y: Vec<E::G1>,
    pub tree_x: EqTree<E::ScalarField>,
    pub tree_y: EqTree<E::ScalarField>,
    pub tree_z: EqTree<E::ScalarField>,
    pub f_star: MultilinearPolynomial<E::ScalarField>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Accumulator3<E: Pairing> {
    pub witness: Acc3Witness<E>,
    pub instance: Acc3Instance<E>,
}

impl<E: Pairing> Accumulator3<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
{
    pub fn setup<R: RngCore>(pc_srs: KZH3SRS<E>, rng: &mut R) -> Acc3SRS<E> {
        Acc3SRS {
            pc_srs: pc_srs.clone(),
            k_x: generate_random_elements::<E, R>(2 * pc_srs.degree_x - 1, rng),
            k_y: generate_random_elements::<E, R>(2 * pc_srs.degree_y - 1, rng),
            k_z: generate_random_elements::<E, R>(2 * pc_srs.degree_z - 1, rng),
            k_prime: E::G1Affine::rand(rng),
        }
    }

    pub fn new(instance: &Acc3Instance<E>, witness: &Acc3Witness<E>) -> Accumulator3<E> {
        Accumulator3 {
            witness: witness.clone(),
            instance: instance.clone(),
        }
    }

    /// the fiat-shamir challenge is computed as part the transcript operations via hashing two accumulator instances and proof Q
    pub fn compute_fiat_shamir_challenge(
        transcript: &mut Transcript<E::ScalarField>,
        instance_1: &Acc3Instance<E>,
        instance_2: &Acc2Instance<E>,
        Q: Acc3Proof<E>,
    ) -> E::ScalarField {
        // add the instances to the transcript
        transcript.append_scalars(b"instance 1", instance_1.to_sponge_field_elements().as_slice());
        transcript.append_scalars(b"instance 2", instance_2.to_sponge_field_elements().as_slice());

        // convert the proof Q into scalar field elements and add to the transcript
        let (p1, p2) = convert_affine_to_scalars::<E>(Q.0);
        transcript.append_scalars(b"Q", &[p1, p2]);
        let (p1, p2) = convert_affine_to_scalars::<E>(Q.1);
        transcript.append_scalars(b"Q", &[p1, p2]);
        let (p1, p2) = convert_affine_to_scalars::<E>(Q.2);
        transcript.append_scalars(b"Q", &[p1, p2]);
        let (p1, p2) = convert_affine_to_scalars::<E>(Q.3);
        transcript.append_scalars(b"Q", &[p1, p2]);

        // return the challenge
        transcript.challenge_scalar(b"challenge scalar")
    }
}

// impl function to convert proof into accumulator
impl<E: Pairing> Accumulator3<E> {
    pub fn proof_to_accumulator_instance(
        srs: &Acc3SRS<E>,
        input: &[E::ScalarField],
        output: &E::ScalarField,
        com: &KZH3Commitment<E>,
        open: &KZH3Opening<E>,
    ) -> Acc3Instance<E>
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        let split_input = KZH3::split_input(&srs.pc_srs, input);

        let T = {
            let tree_x = EqTree::new(split_input[0].as_slice());
            let tree_y = EqTree::new(split_input[1].as_slice());
            let tree_z = EqTree::new(split_input[2].as_slice());

            let mut T: E::G1 = E::G1::ZERO;
            T = T.add(E::G1::msm_unchecked(srs.k_x.as_slice(), tree_x.nodes.as_slice()));
            T = T.add(E::G1::msm_unchecked(srs.k_y.as_slice(), tree_y.nodes.as_slice()));
            T = T.add(E::G1::msm_unchecked(srs.k_z.as_slice(), tree_z.nodes.as_slice()));

            assert_eq!(srs.k_x.len(), tree_x.nodes.len(), "invalid size of vector x");
            assert_eq!(srs.k_y.len(), tree_y.nodes.len(), "invalid size of vector y");
            assert_eq!(srs.k_z.len(), tree_y.nodes.len(), "invalid size of vector z");

            T.into()
        };

        Acc3Instance {
            C: com.C,
            C_y: open.C_y,
            T,
            E: (
                E::G1Affine::zero(),
                E::G1Affine::zero(),
                E::G1Affine::zero(),
                E::G1Affine::zero(),
            ),
            x: split_input[0].clone(),
            y: split_input[1].clone(),
            z: split_input[2].clone(),
            output: *output,
        }
    }

    pub fn proof_to_accumulator_witness(
        srs: &Acc3SRS<E>,
        com: KZH3Commitment<E>,
        proof: KZH3Opening<E>,
        input: &[E::ScalarField],
    ) -> Acc3Witness<E>
    where <E as ark_ec::pairing::Pairing>::ScalarField: ark_crypto_primitives::sponge::Absorb,
          <<E as ark_ec::pairing::Pairing>::G1Affine as ark_ec::AffineRepr>::BaseField: ark_ff::PrimeField
    {
        // asserting the sizes are correct
        let split_input = KZH3::split_input(&srs.pc_srs, input);

        Acc3Witness {
            D_x: com.D_x,
            D_y: proof.D_y,
            tree_x: EqTree::new(split_input[0].as_slice()),
            tree_y: EqTree::new(split_input[1].as_slice()),
            tree_z: EqTree::new(split_input[2].as_slice()),
            f_star: proof.f_star,
        }
    }
}

// deciding functions
impl<E: Pairing> Accumulator3<E> {
    pub fn dec_1(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
        let error_tree_x = acc.witness.tree_x.difference(acc.instance.x.as_slice());
        let error_tree_y = acc.witness.tree_y.difference(acc.instance.y.as_slice());
        let error_tree_z = acc.witness.tree_z.difference(acc.instance.z.as_slice());

        let mut res: E::G1 = E::G1::ZERO;
        res = res.add(E::G1::msm_unchecked(srs.k_x.as_slice(), error_tree_x.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_y.as_slice(), error_tree_y.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_z.as_slice(), error_tree_z.nodes.as_slice()));

        res.into()
    }

    pub fn dec_2(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
        let instance = &acc.instance;
        let witness = &acc.witness;

        let e_prime: E::ScalarField = inner_product(
            &witness.f_star.evaluation_over_boolean_hypercube,
            &acc.witness.tree_z.get_leaves(),
        ) - instance.output;

        srs.k_prime.mul(e_prime).into()
    }

    pub fn dec_3(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
        let rhs = E::G1::msm_unchecked(
            srs.pc_srs.H_z.as_slice(),
            acc.witness
                .f_star
                .evaluation_over_boolean_hypercube
                .as_slice(),
        );

        let lhs = E::G1::msm_unchecked(
            acc.witness.D_y.iter().map(|g| g.clone().into()).collect::<Vec<_>>().as_slice(),
            acc.witness.tree_y.get_leaves(),
        );

        rhs.add(lhs.neg()).into()
    }

    pub fn dec_4(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
        let lhs = E::G1::msm_unchecked(
            acc.witness.D_x.iter().map(|g| g.clone().into()).collect::<Vec<_>>().as_slice(),
            acc.witness.tree_x.get_leaves(),
        );

        acc.instance.C_y.add(lhs.neg()).into()
    }

    pub fn decide(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) {
        let instance = &acc.instance;
        let witness = &acc.witness;

        // first condition
        let pairing_lhs = E::multi_pairing(&witness.D_x, &srs.pc_srs.V_x);
        let pairing_rhs = E::pairing(instance.C, srs.pc_srs.v);

        assert_eq!(pairing_lhs, pairing_rhs, "first condition fails");

        // second condition
        let ip_rhs = instance.T;
        let ip_lhs = {
            // Concatenate bases and scalars
            let mut combined_bases = Vec::with_capacity(
                srs.k_x.len() + srs.k_y.len() + srs.k_z.len()
            );
            let mut combined_scalars = Vec::with_capacity(
                witness.tree_x.nodes.len() + witness.tree_y.nodes.len() + witness.tree_z.nodes.len(),
            );

            combined_bases.extend_from_slice(srs.k_x.as_slice());
            combined_bases.extend_from_slice(srs.k_y.as_slice());
            combined_bases.extend_from_slice(srs.k_z.as_slice());

            combined_scalars.extend_from_slice(witness.tree_x.nodes.as_slice());
            combined_scalars.extend_from_slice(witness.tree_y.nodes.as_slice());
            combined_scalars.extend_from_slice(witness.tree_z.nodes.as_slice());

            // Perform a single MSM
            E::G1::msm_unchecked(combined_bases.as_slice(), combined_scalars.as_slice())
        };

        assert_eq!(ip_rhs, ip_lhs.into(), "second condition fails");

        // third condition
        assert_eq!(
            (Self::dec_1(srs, acc), Self::dec_2(srs, acc), Self::dec_3(srs, acc), Self::dec_4(srs, acc)),
            acc.instance.E,
            "third condition fails"
        );

        // forth condition
        let pairing_lhs = E::multi_pairing(&witness.D_y, &srs.pc_srs.V_y);
        let pairing_rhs = E::pairing(instance.C_y, srs.pc_srs.v);
        assert_eq!(pairing_lhs, pairing_rhs, "forth condition fails");
    }
}

#[cfg(test)]
mod test {
    use crate::constant_for_curves::{ScalarField as F, E};
    use crate::kzh::kzh3::{KZH3, KZH3SRS};
    use crate::kzh::KZH;
    use crate::kzh_fold::kzh_3_fold::Accumulator3;
    use crate::math::Math;
    use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
    use ark_std::UniformRand;
    use rand::thread_rng;

    #[test]
    fn test() {
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
        let output = polynomial.evaluate(input.as_slice());

        // commit to the polynomial
        let c = KZH3::commit(&srs, &polynomial);

        // open it
        let open = KZH3::open(&srs, input.as_slice(), &c, &polynomial);

        // verify the commit
        KZH3::verify(&srs, input.as_slice(), &output, &c, &open);

        let acc_srs = Accumulator3::setup(srs, &mut thread_rng());

        let acc_instance = Accumulator3::proof_to_accumulator_instance(
            &acc_srs,
            input.as_slice(),
            &output,
            &c,
            &open
        );

        let acc_witness = Accumulator3::proof_to_accumulator_witness(
            &acc_srs,
            c,
            open,
            input.as_slice(),
        );

        let acc = Accumulator3::new(&acc_instance, &acc_witness);

        Accumulator3::decide(&acc_srs, &acc);
    }
}

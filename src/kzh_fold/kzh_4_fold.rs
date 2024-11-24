/*use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{AdditiveGroup, Field, PrimeField, Zero};
use ark_poly::EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use ark_std::UniformRand;
use rand::{Rng, RngCore};
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use std::ops::{Add, Mul, Neg, Sub};

use crate::gadgets::non_native::util::convert_affine_to_scalars;
use crate::kzh::kzh4::{KZH4Opening, KZH4SRS};
use crate::kzh::KZH;
use crate::kzh_fold::eq_tree::EqTree;
use crate::kzh_fold::{generate_random_elements, generic_linear_combination};
use crate::math::Math;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use crate::utils::inner_product;

#[derive(Clone, Debug)]
pub struct Acc4SRS<E: Pairing> {
    // vector of size 2 * degree_x - 1
    pub k_x: Vec<E::G1Affine>,

    // vector of size 2 * degree_y - 1
    pub k_y: Vec<E::G1Affine>,

    // vector of size 2 * degree_z - 1
    pub k_z: Vec<E::G1Affine>,

    // vector of size 2 * degree_t - 1
    pub k_t: Vec<E::G1Affine>,

    pub k_prime: E::G1Affine,

    pub pc_srs: KZH4SRS<E>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Acc4Instance<E: Pairing> {
    pub C: E::G1Affine,
    pub T: E::G1Affine,
    pub E: E::G1Affine,

    // vector of length log2(degree_x)
    pub x: Vec<E::ScalarField>,

    // vector of length log2(degree_y)
    pub y: Vec<E::ScalarField>,

    // vector of length log2(degree_z)
    pub z: Vec<E::ScalarField>,

    // vector of length log2(degree_t)
    pub t: Vec<E::ScalarField>,

    pub output: E::ScalarField,
}

impl<E: Pairing> Acc4Instance<E>
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
        let (e_x, e_y) = convert_affine_to_scalars::<E>(self.E);

        // Extend the destination vector with the computed values
        dest.extend(vec![c_x, c_y, t_x, t_y, e_x, e_y]);

        // Extend with other scalar fields
        dest.extend(self.x.clone());
        dest.extend(self.y.clone());
        dest.extend(self.z.clone());
        dest.extend(self.t.clone());
        dest.push(self.output);

        dest
    }
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Acc4Witness<E: Pairing> {
    /// size of degree_x
    pub D_x: Vec<E::G1Affine>,

    /// size of degree_y
    pub D_y: Vec<E::G1Affine>,

    /// size of degree_z
    pub D_z: Vec<E::G1Affine>,

    /// MLP of degree_t
    pub f_star_poly: MultilinearPolynomial<E::ScalarField>,

    pub tree_x: EqTree<E::ScalarField>,
    pub tree_y: EqTree<E::ScalarField>,
    pub tree_z: EqTree<E::ScalarField>,
    pub tree_t: EqTree<E::ScalarField>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Accumulator4<E: Pairing> {
    pub witness: Acc4Witness<E>,
    pub instance: Acc4Instance<E>,
}

impl<E: Pairing> Accumulator4<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
{
    pub fn setup<R: RngCore>(pc_srs: KZH4SRS<E>, rng: &mut R) -> Acc4SRS<E> {
        Acc4SRS {
            pc_srs: pc_srs.clone(),
            k_x: generate_random_elements::<E, R>(2 * pc_srs.degree_x - 1, rng),
            k_y: generate_random_elements::<E, R>(2 * pc_srs.degree_y - 1, rng),
            k_z: generate_random_elements::<E, R>(2 * pc_srs.degree_z - 1, rng),
            k_t: generate_random_elements::<E, R>(2 * pc_srs.degree_t - 1, rng),
            k_prime: E::G1Affine::rand(rng),
        }
    }

    pub fn new(instance: &Acc4Instance<E>, witness: &Acc4Witness<E>) -> Accumulator4<E> {
        Accumulator4 {
            witness: witness.clone(),
            instance: instance.clone(),
        }
    }

    /// the fiat-shamir challenge is computed as part the transcript operations via hashing two accumulator instances and proof Q
    pub fn compute_fiat_shamir_challenge(
        transcript: &mut Transcript<E::ScalarField>,
        instance_1: &Acc4Instance<E>,
        instance_2: &Acc4Instance<E>,
        Q: E::G1Affine,
    ) -> E::ScalarField {
        // add the instances to the transcript
        transcript.append_scalars(b"instance 1", instance_1.to_sponge_field_elements().as_slice());
        transcript.append_scalars(b"instance 2", instance_2.to_sponge_field_elements().as_slice());

        // convert the proof Q into scalar field elements and add to the transcript
        let (p1, p2) = convert_affine_to_scalars::<E>(Q);
        transcript.append_scalars(b"Q", &[p1, p2]);

        // return the challenge
        transcript.challenge_scalar(b"challenge scalar")
    }

    /// Given public data for the opening p(x, y) = z, return an accumulator instance
    pub fn new_accumulator_instance_from_fresh_kzh_instance(
        srs: &Acc4SRS<E>,
        C: &E::G1Affine,
        x: &[E::ScalarField],
        y: &[E::ScalarField],
        z: &[E::ScalarField],
        t: &[E::ScalarField],
        output: &E::ScalarField,
    ) -> Acc4Instance<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.pc_srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.pc_srs.degree_y, "invalid size of vector y");
        assert_eq!(1 << z.len(), srs.pc_srs.degree_z, "invalid size of vector z");
        assert_eq!(1 << t.len(), srs.pc_srs.degree_t, "invalid size of vector t");

        let tree_x = EqTree::new(x);
        let tree_y = EqTree::new(y);
        let tree_z = EqTree::new(z);
        let tree_t = EqTree::new(t);


        let mut T: E::G1 = E::G1::ZERO;
        T = T.add(E::G1::msm_unchecked(srs.k_x.as_slice(), tree_x.nodes.as_slice()));
        T = T.add(E::G1::msm_unchecked(srs.k_y.as_slice(), tree_y.nodes.as_slice()));
        T = T.add(E::G1::msm_unchecked(srs.k_z.as_slice(), tree_z.nodes.as_slice()));
        T = T.add(E::G1::msm_unchecked(srs.k_t.as_slice(), tree_t.nodes.as_slice()));

        assert_eq!(srs.k_x.len(), tree_x.nodes.len(), "invalid size of vector x");
        assert_eq!(srs.k_y.len(), tree_y.nodes.len(), "invalid size of vector y");
        assert_eq!(srs.k_z.len(), tree_z.nodes.len(), "invalid size of vector z");
        assert_eq!(srs.k_t.len(), tree_t.nodes.len(), "invalid size of vector t");


        Acc4Instance {
            C: *C,
            T: T.into(),
            E: E::G1Affine::zero(),
            x: x.to_vec(),
            y: y.to_vec(),
            z: z.to_vec(),
            t: t.to_vec(),
            output: output.clone(),
        }
    }

    pub fn new_accumulator_witness_from_fresh_kzh_witness(
        srs: &Acc4SRS<E>,
        proof: KZH4Opening<E>,
        x: &[E::ScalarField],
        y: &[E::ScalarField],
        z: &[E::ScalarField],
        t: &[E::ScalarField],
    ) -> Acc4Witness<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.pc_srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.pc_srs.degree_y, "invalid size of vector y");
        assert_eq!(1 << z.len(), srs.pc_srs.degree_z, "invalid size of vector z");
        assert_eq!(1 << t.len(), srs.pc_srs.degree_t, "invalid size of vector t");
        assert_eq!(proof.D_x.len(), srs.pc_srs.degree_x, "invalid proof size");
        assert_eq!(proof.D_y.len(), srs.pc_srs.degree_y, "invalid proof size");
        assert_eq!(proof.D_z.len(), srs.pc_srs.degree_z, "invalid proof size");

        // return a fresh instance by simply computing the two EqTrees
        Acc4Witness {
            D_x: proof.D_x.iter().map(|el| el.into_affine()).collect(),
            D_y: proof.D_y.iter().map(|el| el.into_affine()).collect(),
            D_z: proof.D_z.iter().map(|el| el.into_affine()).collect(),
            f_star_poly: proof.f_star,
            tree_x: EqTree::new(x),
            tree_y: EqTree::new(y),
            tree_z: EqTree::new(z),
            tree_t: EqTree::new(t),
        }
    }

    pub fn prove(
        srs: &Acc4SRS<E>,
        acc_1: &Accumulator4<E>,
        acc_2: &Accumulator4<E>,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> (Acc4Instance<E>, Acc4Witness<E>, E::G1Affine)
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
    {
        // unwrap the instances and witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        // compute the quotient variable Q
        let Q: E::G1Affine = Self::helper_function_Q(srs, acc_1, acc_2);

        // the transcript is cloned because is used twice, once to compute the fiat_shamir challenge beta directly
        // once to call Self::verify() to get the new accumulated instance, despite that they generate the same randomness
        // and one can do fewer hashes by computing accumulated instance directly, it's for the sake of cleaner code.
        let mut transcript_clone = transcript.clone();

        // get challenge beta
        let beta = Accumulator4::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, Q);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        // get the accumulated new_instance
        let new_instance = Self::verify(srs, instance_1, instance_2, Q, &mut transcript_clone);

        // get the accumulated witness
        let new_witness = Acc4Witness {
            D_x: generic_linear_combination(
                &witness_1.D_x,
                &witness_2.D_x,
                |d_1, d_2| d_1.mul(one_minus_beta).add(d_2.mul(beta)).into_affine(),
            ),
            D_y: generic_linear_combination(
                &witness_1.D_y,
                &witness_2.D_y,
                |d_1, d_2| d_1.mul(one_minus_beta).add(d_2.mul(beta)).into_affine(),
            ),
            D_z: generic_linear_combination(
                &witness_1.D_z,
                &witness_2.D_z,
                |d_1, d_2| d_1.mul(one_minus_beta).add(d_2.mul(beta)).into_affine(),
            ),
            f_star_poly: MultilinearPolynomial::linear_combination(
                &witness_1.f_star_poly,
                &witness_2.f_star_poly,
                |a, b| a * one_minus_beta + b * beta,
            ),
            tree_x: EqTree::linear_combination(
                &witness_1.tree_x,
                &witness_2.tree_x,
                |a, b| a * one_minus_beta + b * beta,
            ),
            tree_y: EqTree::linear_combination(
                &witness_1.tree_y,
                &witness_2.tree_y,
                |a, b| a * one_minus_beta + b * beta,
            ),
            tree_z: EqTree::linear_combination(
                &witness_1.tree_z,
                &witness_2.tree_z,
                |a, b| a * one_minus_beta + b * beta,
            ),
            tree_t: EqTree::linear_combination(
                &witness_1.tree_t,
                &witness_2.tree_t,
                |a, b| a * one_minus_beta + b * beta,
            ),
        };

        (new_instance, new_witness, Q)
    }

    pub fn verify(
        srs: &Acc4SRS<E>,
        instance_1: &Acc4Instance<E>,
        instance_2: &Acc4Instance<E>,
        Q: E::G1Affine,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> Acc4Instance<E> {
        // compute the fiat-shamir challenge
        let beta = Accumulator4::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, Q);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        let new_error_term: E::G1Affine = {
            let mut res = instance_1.E.mul(one_minus_beta);
            res = res.add(instance_2.E.mul(beta));
            res.add(Q.mul(one_minus_beta * beta)).into()
        };

        Acc4Instance {
            C: {
                let res = instance_1.C.mul(one_minus_beta);
                res.add(instance_2.C.mul(beta)).into()
            },
            T: {
                let res = instance_1.T.mul(one_minus_beta);
                res.add(instance_2.T.mul(beta)).into()
            },
            x: generic_linear_combination(
                &instance_1.x,
                &instance_2.x,
                |e1, e2| e1 * one_minus_beta + e2 * beta,
            ),
            y: generic_linear_combination(
                &instance_1.y,
                &instance_2.y,
                |e1, e2| e1 * one_minus_beta + e2 * beta,
            ),
            z: generic_linear_combination(
                &instance_1.z,
                &instance_2.z,
                |e1, e2| e1 * one_minus_beta + e2 * beta,
            ),
            t: generic_linear_combination(
                &instance_1.t,
                &instance_2.t,
                |e1, e2| e1 * one_minus_beta + e2 * beta,
            ),
            output: instance_1.output * one_minus_beta + instance_2.output * beta,
            E: new_error_term,
        }
    }

    pub fn helper_function_Q(srs: &Acc4SRS<E>, acc_1: &Accumulator4<E>, acc_2: &Accumulator4<E>) -> E::G1Affine {
        // unwrap the instances/witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        assert_eq!(witness_1.f_star_poly.num_variables, witness_2.f_star_poly.num_variables);
        assert_eq!(witness_1.tree_x.depth, witness_2.tree_x.depth);
        assert_eq!(witness_1.tree_y.depth, witness_2.tree_y.depth);
        assert_eq!(witness_1.tree_z.depth, witness_2.tree_z.depth);
        assert_eq!(witness_1.tree_t.depth, witness_2.tree_t.depth);
        assert_eq!(witness_1.D_x.len(), witness_2.D_x.len());
        assert_eq!(witness_1.D_y.len(), witness_2.D_y.len());
        assert_eq!(witness_1.D_z.len(), witness_2.D_z.len());


        // value of 2 in the field
        let two = E::ScalarField::from(2u128);

        // build the accumulator from linear combination to run helper_function_V on it
        let temp_acc = Accumulator4 {
            witness: Acc4Witness {
                D_x: generic_linear_combination(
                    &witness_2.D_x,
                    &witness_1.D_x,
                    |d2, d1| d2.mul(two).sub(d1).into_affine(),
                ),
                D_y: generic_linear_combination(
                    &witness_2.D_y,
                    &witness_1.D_y,
                    |d2, d1| d2.mul(two).sub(d1).into_affine(),
                ),
                D_z: generic_linear_combination(
                    &witness_2.D_z,
                    &witness_1.D_z,
                    |d2, d1| d2.mul(two).sub(d1).into_affine(),
                ),
                f_star_poly: MultilinearPolynomial::linear_combination(
                    &witness_2.f_star_poly,
                    &witness_1.f_star_poly,
                    |e2, e1| e2 * two - e1,
                ),
                tree_x: EqTree::linear_combination(
                    &witness_2.tree_x,
                    &witness_1.tree_x,
                    |w2, w1| w2 * two - w1,
                ),
                tree_y: EqTree::linear_combination(
                    &witness_2.tree_y,
                    &witness_1.tree_y,
                    |w2, w1| w2 * two - w1,
                ),
                tree_z: EqTree::linear_combination(
                    &witness_2.tree_z,
                    &witness_1.tree_z,
                    |w2, w1| w2 * two - w1,
                ),
                tree_t: EqTree::linear_combination(
                    &witness_2.tree_t,
                    &witness_1.tree_t,
                    |w2, w1| w2 * two - w1,
                ),
            },
            instance: Acc4Instance {
                // not used by helper function, so we simply pass them as zero or any other random element
                C: E::G1Affine::zero(),
                T: E::G1Affine::zero(),
                E: E::G1Affine::zero(),

                // used by the helper function
                x: generic_linear_combination(
                    &instance_2.x,
                    &instance_1.x,
                    |e2, e1| e2 * two - e1,
                ),
                y: generic_linear_combination(
                    &instance_2.y,
                    &instance_1.y,
                    |e2, e1| e2 * two - e1,
                ),
                z: generic_linear_combination(
                    &instance_2.z,
                    &instance_1.z,
                    |e2, e1| e2 * two - e1,
                ),
                t: generic_linear_combination(
                    &instance_2.t,
                    &instance_1.t,
                    |e2, e1| e2 * two - e1,
                ),
                output: two * instance_2.output - instance_1.output,
            },
        };

        let mut res = Self::helper_function_decide(&srs, &temp_acc);
        res = res.add(instance_1.E).into();
        res = res.sub(instance_2.E).into();
        res = res.sub(instance_2.E).into();

        // -1/2 in the scalar field
        let minus_one_over_two: E::ScalarField = two.neg().inverse().unwrap();
        res.mul(minus_one_over_two).into()
    }

    pub fn decide(srs: &Acc4SRS<E>, acc: &Accumulator4<E>) -> bool {
        let instance = &acc.instance;
        let witness = &acc.witness;

        // first condition
        let pairing_lhs = E::multi_pairing(&witness.D_x, &srs.pc_srs.V_x);
        let pairing_rhs = E::pairing(instance.C, srs.pc_srs.v);

        // third condition
        let verify_lhs = Self::helper_function_decide(srs, acc);
        let verify_rhs = instance.E;


        // return and of all three conditions
        (pairing_lhs == pairing_rhs) && (verify_lhs == verify_rhs)
    }

    pub fn helper_function_decide(srs: &Acc4SRS<E>, acc: &Accumulator4<E>) -> E::G1Affine {
        let instance = &acc.instance;
        let witness = &acc.witness;

        let e_prime: E::ScalarField = inner_product(
            &witness.f_star_poly.evaluation_over_boolean_hypercube,
            &acc.witness.tree_t.get_leaves(),
        ) - &instance.output;

        // Optimize E_G computation by doing one big MSM
        let E_G = {
            // Concatenate the scalar vectors
            let scalars: Vec<E::ScalarField> = {
                // Prepare the right-hand-side scalars to be negated
                let tree_x_leaves = acc.witness.tree_x.get_leaves();
                let negated_tree_x_leaves: Vec<E::ScalarField> = tree_x_leaves.iter().map(|&x| -x).collect();

                let tree_y_leaves = acc.witness.tree_y.get_leaves();
                let negated_tree_y_leaves: Vec<E::ScalarField> = tree_y_leaves.iter().map(|&y| -y).collect();

                let tree_z_leaves = acc.witness.tree_z.get_leaves();
                let negated_tree_z_leaves: Vec<E::ScalarField> = tree_z_leaves.iter().map(|&z| -z).collect();

                let mut res = Vec::new();
                res.extend_from_slice(
                    witness
                        .f_star_poly
                        .evaluation_over_boolean_hypercube
                        .as_slice(),
                );
                res.extend_from_slice(&negated_tree_x_leaves);
                res.extend_from_slice(&negated_tree_y_leaves);
                res.extend_from_slice(&negated_tree_z_leaves);

                res
            };


            // Concatenate the base point vectors
            let bases: Vec<E::G1Affine> = {
                let mut res = Vec::new();
                res.extend_from_slice(srs.pc_srs.H_t.as_slice());
                res.extend_from_slice(witness.D_x.as_slice());
                res.extend_from_slice(witness.D_y.as_slice());
                res.extend_from_slice(witness.D_z.as_slice());

                res
            };

            E::G1::msm_unchecked(&bases, &scalars)
        };

        let error_tree_x = acc.witness.tree_x.difference(acc.instance.x.as_slice());
        let error_tree_y = acc.witness.tree_y.difference(acc.instance.y.as_slice());
        let error_tree_z = acc.witness.tree_z.difference(acc.instance.z.as_slice());
        let error_tree_t = acc.witness.tree_t.difference(acc.instance.t.as_slice());

        let mut res: E::G1 = E_G.clone();
        res = res.add(E::G1::msm_unchecked(srs.k_x.as_slice(), error_tree_x.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_y.as_slice(), error_tree_y.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_z.as_slice(), error_tree_z.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_t.as_slice(), error_tree_t.nodes.as_slice()));
        res.add(srs.k_prime.mul(e_prime)).into()
    }
}


pub mod test {
    use ark_ec::pairing::Pairing;
    use rand::thread_rng;

    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use crate::kzh::kzh4::KZH4;

    #[test]
    fn test_accumulator_end_to_end() {
        let (degree_x, degree_y, degree_z, degree_t) = (4usize, 2usize, 16usize, 8usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2() + degree_t.log_2();

        let input: Vec<ScalarField> = {
            let mut res = Vec::new();
            for _ in 0..num_vars {
                res.push(ScalarField::rand(&mut thread_rng()));
            }
            res
        };
        // build the srs
        let srs: KZH4SRS<E> = KZH4::setup(num_vars, &mut thread_rng());
        let acc_srs: Acc4SRS<E> = Accumulator4::setup(srs.clone(), &mut thread_rng());

        let split_input = KZH4::split_input(&srs, input.as_slice());

        // build a random polynomials
        let polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        // evaluate polynomial
        let output = polynomial.evaluate(input.as_slice());

        // commit to the polynomial
        let com = KZH4::commit(&srs, &polynomial);

        // open it
        let open = KZH4::open(&srs, input.as_slice(), &com, &polynomial);

        let acc_instance = Accumulator4::new_accumulator_instance_from_fresh_kzh_instance(
            &acc_srs,
            &com.C,
            split_input[0].as_slice(),
            split_input[1].as_slice(),
            split_input[2].as_slice(),
            split_input[3].as_slice(),
            &output,
        );

        let acc_witness = Accumulator4::new_accumulator_witness_from_fresh_kzh_witness(
            &acc_srs,
            open,
            split_input[0].as_slice(),
            split_input[1].as_slice(),
            split_input[2].as_slice(),
            split_input[3].as_slice(),
        );

        let acc = Accumulator4::new(&acc_instance, &acc_witness);

        assert!(Accumulator4::decide(&acc_srs, &acc));
    }
}

 */


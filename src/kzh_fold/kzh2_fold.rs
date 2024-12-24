use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{AdditiveGroup, Field, PrimeField, Zero};
use ark_poly::EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use ark_std::UniformRand;
use ark_std::{end_timer, start_timer};
use rand::{Rng, RngCore};
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use std::ops::{Add, Mul, Neg, Sub};

use crate::kzh_fold::eq_tree::EqTree;
use crate::kzh_fold::{generate_random_elements, generic_linear_combination};
use crate::gadgets::non_native::util::convert_affine_to_scalars;
use crate::kzh::KZH;
use crate::math::Math;
use crate::kzh::kzh2::{KZH2, KZH2Opening, KZH2SRS};
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use crate::utils::inner_product;

#[derive(Clone, Debug)]
pub struct Acc2SRS<E: Pairing> {
    // vector of size 2 * degree_x - 1
    pub k_x: Vec<E::G1Affine>,

    // vector of size 2 * degree_y - 1
    pub k_y: Vec<E::G1Affine>,

    pub k_prime: E::G1Affine,
    pub pc_srs: KZH2SRS<E>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Acc2Instance<E: Pairing> {
    pub C: E::G1Affine,
    pub T: E::G1Affine,
    pub E: E::G1Affine,

    // vector of length log2(degree_x)
    pub x: Vec<E::ScalarField>,

    // vector of length log2(degree_y)
    pub y: Vec<E::ScalarField>,

    pub z: E::ScalarField,
}

impl<E: Pairing> Acc2Instance<E>
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
        dest.push(self.z);

        dest
    }
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Acc2Witness<E: Pairing> {
    /// size of degree_x
    pub D_x: Vec<E::G1Affine>,

    /// MLP of degree_y
    pub f_star_poly: MultilinearPolynomial<E::ScalarField>,

    pub tree_x: EqTree<E::ScalarField>,

    pub tree_y: EqTree<E::ScalarField>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Accumulator2<E: Pairing> {
    pub witness: Acc2Witness<E>,
    pub instance: Acc2Instance<E>,
}

impl<E: Pairing> Accumulator2<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
{
    pub fn setup<R: RngCore>(pc_srs: KZH2SRS<E>, rng: &mut R) -> Acc2SRS<E> {
        Acc2SRS {
            pc_srs: pc_srs.clone(),
            k_x: generate_random_elements::<E, R>(2 * pc_srs.degree_x - 1, rng),
            k_y: generate_random_elements::<E, R>(2 * pc_srs.degree_y - 1, rng),
            k_prime: E::G1Affine::rand(rng),
        }
    }

    pub fn new(instance: &Acc2Instance<E>, witness: &Acc2Witness<E>) -> Accumulator2<E> {
        Accumulator2 {
            witness: witness.clone(),
            instance: instance.clone(),
        }
    }

    /// the fiat-shamir challenge is computed as part the transcript operations via hashing two accumulator instances and proof Q
    pub fn compute_fiat_shamir_challenge(transcript: &mut Transcript<E::ScalarField>, instance_1: &Acc2Instance<E>, instance_2: &Acc2Instance<E>, Q: E::G1Affine) -> E::ScalarField {
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
    pub fn proof_to_accumulator_instance(
        srs: &Acc2SRS<E>,
        C: &E::G1Affine,
        x: &[E::ScalarField],
        y: &[E::ScalarField],
        z: &E::ScalarField,
    ) -> Acc2Instance<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.pc_srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.pc_srs.degree_y, "invalid size of vector y");

        let tree_x = EqTree::new(x);
        let tree_y = EqTree::new(y);

        let mut T: E::G1 = E::G1::ZERO;
        T = T.add(E::G1::msm_unchecked(srs.k_x.as_slice(), tree_x.nodes.as_slice()));
        T = T.add(E::G1::msm_unchecked(srs.k_y.as_slice(), tree_y.nodes.as_slice()));

        assert_eq!(srs.k_x.len(), tree_x.nodes.len(), "invalid size of vector x");
        assert_eq!(srs.k_y.len(), tree_y.nodes.len(), "invalid size of vector y");


        Acc2Instance {
            C: *C,
            T: T.into(),
            E: E::G1Affine::zero(),
            x: x.to_vec(),
            y: y.to_vec(),
            z: z.clone(),
        }
    }

    pub fn proof_to_accumulator_witness(
        srs: &Acc2SRS<E>,
        proof: KZH2Opening<E>,
        x: &[E::ScalarField],
        y: &[E::ScalarField],
    ) -> Acc2Witness<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.pc_srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.pc_srs.degree_y, "invalid size of vector y");
        assert_eq!(proof.D_x.len(), srs.pc_srs.degree_x, "invalid proof size");

        // return a fresh instance by simply computing the two EqTrees
        Acc2Witness {
            D_x: proof.D_x,
            f_star_poly: proof.f_star_poly,
            tree_x: EqTree::new(x),
            tree_y: EqTree::new(y),
        }
    }

    pub fn prove(
        srs: &Acc2SRS<E>,
        acc_1: &Accumulator2<E>,
        acc_2: &Accumulator2<E>,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> (Acc2Instance<E>, Acc2Witness<E>, E::G1Affine)
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
        let beta = Accumulator2::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, Q);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        // get the accumulated new_instance
        let new_instance = Self::verify(srs, instance_1, instance_2, Q, &mut transcript_clone);

        // get the accumulated witness
        let new_witness = Acc2Witness {
            D_x: generic_linear_combination(
                &witness_1.D_x,
                &witness_2.D_x,
                |d_1, d_2| d_1.mul(one_minus_beta).add(d_2.mul(beta)).into_affine()
            ),
            f_star_poly: MultilinearPolynomial::linear_combination(
                &witness_1.f_star_poly,
                &witness_2.f_star_poly,
                |a, b| a * one_minus_beta + b * beta
            ),
            tree_x: EqTree::linear_combination(
                &witness_1.tree_x,
                &witness_2.tree_x,
                |a, b| a * one_minus_beta + b * beta
            ),
            tree_y: EqTree::linear_combination(
                &witness_1.tree_y,
                &witness_2.tree_y,
                |a, b| a * one_minus_beta + b * beta
            ),
        };

        (new_instance, new_witness, Q)
    }

    pub fn verify(
        srs: &Acc2SRS<E>,
        instance_1: &Acc2Instance<E>,
        instance_2: &Acc2Instance<E>,
        Q: E::G1Affine,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> Acc2Instance<E> {
        // compute the fiat-shamir challenge
        let beta = Accumulator2::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, Q);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        let new_error_term: E::G1Affine = {
            let mut res = instance_1.E.mul(one_minus_beta);
            res = res.add(instance_2.E.mul(beta));
            res.add(Q.mul(one_minus_beta * beta)).into()
        };

        Acc2Instance {
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
                |e1, e2| e1 * one_minus_beta + e2 * beta
            ),
            y: generic_linear_combination(
                &instance_1.y,
                &instance_2.y,
                |e1, e2| e1 * one_minus_beta + e2 * beta
            ),
            z: instance_1.z * one_minus_beta + instance_2.z * beta,
            E: new_error_term,
        }
    }

    pub fn helper_function_Q(srs: &Acc2SRS<E>, acc_1: &Accumulator2<E>, acc_2: &Accumulator2<E>) -> E::G1Affine {
        // unwrap the instances/witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        assert_eq!(witness_1.f_star_poly.num_variables, witness_2.f_star_poly.num_variables);
        assert_eq!(witness_1.tree_x.depth, witness_2.tree_x.depth);
        assert_eq!(witness_1.tree_y.depth, witness_2.tree_y.depth);
        assert_eq!(witness_1.D_x.len(), witness_2.D_x.len());

        // value of 2 in the field
        let two = E::ScalarField::from(2u128);

        // build the accumulator from linear combination to run helper_function_V on it
        let temp_acc = Accumulator2 {
            witness: Acc2Witness {
                D_x: generic_linear_combination(
                    &witness_2.D_x,
                    &witness_1.D_x,
                    |d2, d1| d2.mul(two).sub(d1).into_affine()
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
            },
            instance: Acc2Instance {
                // not used by helper function, so we simply pass them as zero or any other random element
                C: E::G1Affine::zero(),
                T: E::G1Affine::zero(),
                E: E::G1Affine::zero(),

                // used by the helper function
                x: generic_linear_combination(
                    &instance_2.x,
                    &instance_1.x,
                    |e2, e1| e2 * two - e1
                ),
                y: generic_linear_combination(
                    &instance_2.y,
                    &instance_1.y,
                    |e2, e1| e2 * two - e1
                ),
                z: two * instance_2.z - instance_1.z,
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

    pub fn decide(srs: &Acc2SRS<E>, acc: &Accumulator2<E>) -> bool {
        let instance = &acc.instance;
        let witness = &acc.witness;

        // first condition
        let pairing_lhs = E::multi_pairing(&witness.D_x, &srs.pc_srs.V_x);
        let pairing_rhs = E::pairing(instance.C, srs.pc_srs.V_prime);

        // second condition
        let ip_rhs = instance.T;
        let ip_lhs = {
            // Concatenate bases and scalars
            let mut combined_bases = Vec::with_capacity(srs.k_x.len() + srs.k_y.len());
            let mut combined_scalars = Vec::with_capacity(witness.tree_x.nodes.len() + witness.tree_y.nodes.len());

            combined_bases.extend_from_slice(srs.k_x.as_slice());
            combined_bases.extend_from_slice(srs.k_y.as_slice());

            combined_scalars.extend_from_slice(witness.tree_x.nodes.as_slice());
            combined_scalars.extend_from_slice(witness.tree_y.nodes.as_slice());

            // Perform a single MSM
            E::G1::msm_unchecked(combined_bases.as_slice(), combined_scalars.as_slice())
        };

        // third condition
        let verify_lhs = Self::helper_function_decide(srs, acc);
        let verify_rhs = instance.E;

        // return and of all three conditions
        (verify_rhs == verify_lhs.into()) && (ip_lhs == ip_rhs.into()) && (pairing_lhs == pairing_rhs)
    }

    pub fn helper_function_decide(srs: &Acc2SRS<E>, acc: &Accumulator2<E>) -> E::G1Affine {
        let instance = &acc.instance;
        let witness = &acc.witness;

        let e_prime: E::ScalarField = inner_product(
            &witness.f_star_poly.evaluation_over_boolean_hypercube,
            &acc.witness.tree_y.get_leaves(),
        ) - instance.z;

        // Optimize E_G computation by doing one big MSM
        let E_G = {
            // Concatenate the scalar vectors
            let scalars: Vec<E::ScalarField> = {
                // Prepare the right-hand-side scalars to be negated
                let tree_x_leaves = acc.witness.tree_x.get_leaves();
                let negated_tree_x_leaves: Vec<E::ScalarField> = tree_x_leaves.iter().map(|&x| -x).collect();

                let mut res = Vec::new();
                res.extend_from_slice(
                    witness
                        .f_star_poly
                        .evaluation_over_boolean_hypercube
                        .as_slice(),
                );
                res.extend_from_slice(&negated_tree_x_leaves);

                res
            };


            // Concatenate the base point vectors
            let bases: Vec<E::G1Affine> = {
                let mut res = Vec::new();
                res.extend_from_slice(srs.pc_srs.H_y.as_slice());
                res.extend_from_slice(witness.D_x.as_slice());

                res
            };

            E::G1::msm_unchecked(&bases, &scalars)
        };

        let error_tree_x = acc.witness.tree_x.difference(acc.instance.x.as_slice());
        let error_tree_y = acc.witness.tree_y.difference(acc.instance.y.as_slice());

        let mut res: E::G1 = E_G.clone();
        res = res.add(E::G1::msm_unchecked(srs.k_x.as_slice(), error_tree_x.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_y.as_slice(), error_tree_y.nodes.as_slice()));
        res.add(srs.k_prime.mul(e_prime)).into()
    }
}

impl<E: Pairing> Accumulator2<E> {
    /// this function returns a random satisfying accumulator by generating two random frseh accumualtors (KZH openings)
    /// and then accumulating them, so that the error vector wouldn't be zero
    pub fn rand<R: RngCore>(srs: &Acc2SRS<E>, rng: &mut R) -> Accumulator2<E>
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        // random bivariate polynomial
        let polynomial1 = MultilinearPolynomial::rand(
            srs.pc_srs.degree_x.log_2() + srs.pc_srs.degree_y.log_2(),
            rng,
        );
        let polynomial2 = MultilinearPolynomial::rand(
            srs.pc_srs.degree_x.log_2() + srs.pc_srs.degree_y.log_2(),
            rng,
        );

        // random x points
        let mut x1: Vec<E::ScalarField> = Vec::new();
        let mut x2: Vec<E::ScalarField> = Vec::new();
        for _ in 0..srs.pc_srs.degree_x.log_2() {
            x1.push(E::ScalarField::rand(rng));
            x2.push(E::ScalarField::rand(rng));
        }
        // random y points
        let mut y1: Vec<E::ScalarField> = Vec::new();
        let mut y2: Vec<E::ScalarField> = Vec::new();
        for _ in 0..srs.pc_srs.degree_y.log_2() {
            y1.push(E::ScalarField::rand(rng));
            y2.push(E::ScalarField::rand(rng));
        }

        // Get vector: (x1, y1)
        let input_1: Vec<E::ScalarField> = x1.clone().into_iter().chain(y1.clone()).collect();
        // Get vector: (x2, y2)
        let input_2: Vec<E::ScalarField> = x2.clone().into_iter().chain(y2.clone()).collect();

        // get evaluations z1 and z2
        let z1 = polynomial1.evaluate(&input_1);
        let z2 = polynomial2.evaluate(&input_2);

        // commit to the polynomial
        let com1 = KZH2::commit(&srs.pc_srs, &polynomial1);
        let com2 = KZH2::commit(&srs.pc_srs, &polynomial2);

        // open the commitment
        let open1 = KZH2::open(&srs.pc_srs, input_1.as_slice(), &com1, &polynomial1);
        let open2 = KZH2::open(&srs.pc_srs, input_2.as_slice(), &com2, &polynomial2);

        // verify the proof
        KZH2::verify(&srs.pc_srs, &input_1, &z1, &com1, &open1);
        KZH2::verify(&srs.pc_srs, &input_2, &z2, &com2, &open2);

        let instance1 = Accumulator2::proof_to_accumulator_instance(&srs, &com1.C, &x1, &y1, &z1);
        let witness1 = Accumulator2::proof_to_accumulator_witness(&srs, open1, &x1, &y1);
        let instance2 = Accumulator2::proof_to_accumulator_instance(&srs, &com2.C, &x2, &y2, &z2);
        let witness2 = Accumulator2::proof_to_accumulator_witness(&srs, open2, &x2, &y2);

        let acc1 = Accumulator2::new(&instance1, &witness1);
        let acc2 = Accumulator2::new(&instance2, &witness2);

        // verify that the fresh accumulators are satisfied
        debug_assert!(Accumulator2::decide(&srs, &acc1));
        debug_assert!(Accumulator2::decide(&srs, &acc2));

        let mut prover_transcript = Transcript::new(b"new_transcript");

        let (accumulated_instance, accumulated_witness, Q) = Accumulator2::prove(&srs, &acc1, &acc2, &mut prover_transcript);

        let accumulated_acc = Accumulator2::new(&accumulated_instance, &accumulated_witness);

        // verify the accumulated instance is satisfied
        debug_assert!(Accumulator2::decide(&srs, &accumulated_acc));

        accumulated_acc
    }
}

pub mod test {
    use ark_ec::pairing::Pairing;
    use rand::thread_rng;

    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use crate::kzh::kzh2::{KZH2, KZH2SRS};

    #[test]
    fn test_accumulator_end_to_end() {
        let (degree_x, degree_y) = (8usize, 8usize);
        let srs_pcs: KZH2SRS<E> = KZH2::setup((degree_x * degree_y).log_2(), &mut thread_rng());
        let srs = Accumulator2::setup(srs_pcs.clone(), &mut thread_rng());

        let acc1 = Accumulator2::rand(&srs, &mut thread_rng());
        let acc2 = Accumulator2::rand(&srs, &mut thread_rng());

        let mut prover_transcript = Transcript::new(b"new_transcript");
        let mut verifier_transcript = prover_transcript.clone();

        let (instance, witness, Q) = Accumulator2::prove(&srs, &acc1, &acc2, &mut prover_transcript);

        let instance_expected = Accumulator2::verify(&srs, &acc1.instance, &acc2.instance, Q, &mut verifier_transcript);

        assert_eq!(instance, instance_expected);

        let decide_timer = start_timer!(|| "decide");
        assert!(Accumulator2::decide(&srs, &Accumulator2 { witness, instance }));
        end_timer!(decide_timer);
    }

    #[test]
    fn test_accumulator_sizes() {
        // change the degrees later, it takes too long
        let degrees = vec![(2, 2), (4, 4), (8, 8), (16, 16), (32, 32)];
        let rng = &mut thread_rng();

        for (degree_x, degree_y) in degrees {
            // set accumulator sts
            let srs_pcs: KZH2SRS<E> = KZH2::setup((degree_x * degree_y).log_2(), rng);
            let srs_acc = Accumulator2::setup(srs_pcs.clone(), rng);
            // get random accumulator
            let acc = Accumulator2::rand(&srs_acc, rng);

            let witness_len = degree_x * degree_y;
            let witness_polynomial: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::rand(
                degree_x.log_2() + degree_y.log_2(),
                rng,
            );

            println!("witness length: {} ({} bytes compressed):\n\taccumulator size: {} bytes (compressed: {} bytes)\n\t\tinstance compressed: {} bytes\n\t\twitness compressed: {} bytes",
                     witness_len, witness_polynomial.compressed_size(),
                     acc.uncompressed_size(), acc.compressed_size(),
                     acc.instance.compressed_size(),
                     acc.witness.compressed_size()
            );
        }
    }
}


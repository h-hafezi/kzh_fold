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

use crate::accumulation::eq_tree::EqTree;
use crate::accumulation::generate_random_elements;
use crate::gadgets::non_native::util::convert_affine_to_scalars;
use crate::math::Math;
use crate::pcs::multilinear_pcs::{PCSEngine, PCSOpeningProof, PolynomialCommitmentSRS};
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::{AppendToTranscript, Transcript};
use crate::utils::inner_product;

#[derive(Clone, Debug)]
pub struct AccSRS<E: Pairing> {
    // vector of size 2 * degree_x - 1
    pub k_x: Vec<E::G1Affine>,

    // vector of size 2 * degree_y - 1
    pub k_y: Vec<E::G1Affine>,

    pub k_prime: E::G1Affine,
    pub pc_srs: PolynomialCommitmentSRS<E>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct AccInstance<E: Pairing> {
    pub C: E::G1Affine,
    pub T: E::G1Affine,
    pub E: E::G1Affine,

    // vector of length log2(degree_x)
    pub x: Vec<E::ScalarField>,

    // vector of length log2(degree_y)
    pub y: Vec<E::ScalarField>,

    pub z: E::ScalarField,
}

impl<E: Pairing> AccInstance<E>
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
pub struct AccWitness<E: Pairing> {
    /// size of degree_x
    pub vec_D: Vec<E::G1Affine>,

    /// MLP of degree_y
    pub f_star_poly: MultilinearPolynomial<E::ScalarField>,

    pub tree_x: EqTree<E::ScalarField>,

    pub tree_y: EqTree<E::ScalarField>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct Accumulator<E: Pairing> {
    pub witness: AccWitness<E>,
    pub instance: AccInstance<E>,
}

impl<E: Pairing> Accumulator<E>
where
    <E as Pairing>::ScalarField: Absorb,
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
{
    pub fn setup<R: RngCore>(pc_srs: PolynomialCommitmentSRS<E>, rng: &mut R) -> AccSRS<E> {
        AccSRS {
            pc_srs: pc_srs.clone(),
            k_x: generate_random_elements::<E, R>(2 * pc_srs.degree_x - 1, rng),
            k_y: generate_random_elements::<E, R>(2 * pc_srs.degree_y - 1, rng),
            k_prime: E::G1Affine::rand(rng),
        }
    }

    pub fn new(instance: &AccInstance<E>, witness: &AccWitness<E>) -> Accumulator<E> {
        Accumulator {
            witness: witness.clone(),
            instance: instance.clone(),
        }
    }

    /// the fiat-shamir challenge is computed as part the transcript operations via hashing two accumulator instances and proof Q
    pub fn compute_fiat_shamir_challenge(transcript: &mut Transcript<E::ScalarField>, instance_1: &AccInstance<E>, instance_2: &AccInstance<E>, Q: E::G1Affine) -> E::ScalarField {
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
        srs: &AccSRS<E>,
        C: &E::G1Affine,
        x: &[E::ScalarField],
        y: &[E::ScalarField],
        z: &E::ScalarField,
    ) -> AccInstance<E> {
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


        AccInstance {
            C: *C,
            T: T.into(),
            E: E::G1Affine::zero(),
            x: x.to_vec(),
            y: y.to_vec(),
            z: z.clone(),
        }
    }

    pub fn new_accumulator_witness_from_fresh_kzh_witness(
        srs: &AccSRS<E>,
        proof: PCSOpeningProof<E>,
        x: &[E::ScalarField],
        y: &[E::ScalarField],
    ) -> AccWitness<E> {
        // asserting the sizes are correct
        assert_eq!(1 << x.len(), srs.pc_srs.degree_x, "invalid size of vector x");
        assert_eq!(1 << y.len(), srs.pc_srs.degree_y, "invalid size of vector y");
        assert_eq!(proof.vec_D.len(), srs.pc_srs.degree_x, "invalid proof size");

        // return a fresh instance by simply computing the two EqTrees
        AccWitness {
            vec_D: proof.vec_D,
            f_star_poly: proof.f_star_poly,
            tree_x: EqTree::new(x),
            tree_y: EqTree::new(y),
        }
    }

    pub fn prove(
        srs: &AccSRS<E>,
        acc_1: &Accumulator<E>,
        acc_2: &Accumulator<E>,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> (AccInstance<E>, AccWitness<E>, E::G1Affine)
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
        let beta = Accumulator::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, Q);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        // get the accumulated new_instance
        let new_instance = Self::verify(srs, instance_1, instance_2, Q, &mut transcript_clone);

        // get the accumulated witness
        let new_witness = AccWitness {
            vec_D: witness_1.vec_D.par_iter().zip(witness_2.vec_D.par_iter())
                .map(|(&d_1, &d_2)| d_1.mul(one_minus_beta).add(d_2.mul(beta)).into_affine())
                .collect(),
            f_star_poly: MultilinearPolynomial {
                num_variables: witness_1.f_star_poly.num_variables,
                evaluation_over_boolean_hypercube: witness_1.f_star_poly.evaluation_over_boolean_hypercube.par_iter()
                    .zip(witness_2.f_star_poly.evaluation_over_boolean_hypercube.par_iter())
                    .map(
                        |(&a, &b)|
                            a * (one_minus_beta) + (b * beta)
                    )
                    .collect(),
                len: witness_1.f_star_poly.len(),
            },
            tree_x: EqTree {
                nodes: witness_1.tree_x.nodes.par_iter()
                    .zip(witness_2.tree_x.nodes.par_iter())
                    .map(
                        |(&a, &b)|
                            a * (one_minus_beta) + (b * beta)
                    )
                    .collect(),
                depth: witness_1.tree_x.depth,
            },
            tree_y: EqTree {
                nodes: witness_1.tree_y.nodes.par_iter()
                    .zip(witness_2.tree_y.nodes.par_iter())
                    .map(
                        |(&a, &b)|
                            a * (one_minus_beta) + (b * beta)
                    )
                    .collect(),
                depth: witness_1.tree_y.depth,
            },
        };

        (new_instance, new_witness, Q)
    }

    pub fn verify(
        srs: &AccSRS<E>,
        instance_1: &AccInstance<E>,
        instance_2: &AccInstance<E>,
        Q: E::G1Affine,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> AccInstance<E> {
        // compute the fiat-shamir challenge
        let beta = Accumulator::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, Q);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        let new_error_term: E::G1Affine = {
            let mut res = instance_1.E.mul(one_minus_beta);
            res = res.add(instance_2.E.mul(beta));
            res.add(Q.mul(one_minus_beta * beta)).into()
        };

        AccInstance {
            C: {
                let res = instance_1.C.mul(one_minus_beta);
                res.add(instance_2.C.mul(beta)).into()
            },
            T: {
                let res = instance_1.T.mul(one_minus_beta);
                res.add(instance_2.T.mul(beta)).into()
            },
            x: instance_1.x.iter()
                .zip(instance_2.x.iter())
                .map(|(&e1, &e2)| e1 * one_minus_beta + e2 * beta)
                .collect(),
            y: instance_1.y.iter()
                .zip(instance_2.y.iter())
                .map(|(&e1, &e2)| e1 * one_minus_beta + e2 * beta)
                .collect(),
            z: instance_1.z * one_minus_beta + instance_2.z * beta,
            E: new_error_term,
        }
    }

    pub fn helper_function_Q(srs: &AccSRS<E>, acc_1: &Accumulator<E>, acc_2: &Accumulator<E>) -> E::G1Affine {
        // unwrap the instances/witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        assert_eq!(witness_1.f_star_poly.num_variables, witness_2.f_star_poly.num_variables);
        assert_eq!(witness_1.tree_x.depth, witness_2.tree_x.depth);
        assert_eq!(witness_1.tree_y.depth, witness_2.tree_y.depth);
        assert_eq!(witness_1.vec_D.len(), witness_2.vec_D.len());

        // value of 2 in the field
        let two = E::ScalarField::from(2u128);

        // build the accumulator from linear combination to run helper_function_V on it
        let temp_acc = Accumulator {
            witness: AccWitness {
                vec_D: witness_2.vec_D.par_iter()
                    .zip(witness_1.vec_D.par_iter())
                    .map(|(&d2, &d1)| d2.mul(two).sub(d1).into_affine())
                    .collect(),
                f_star_poly: MultilinearPolynomial {
                    num_variables: witness_1.f_star_poly.num_variables,
                    evaluation_over_boolean_hypercube: witness_2.f_star_poly.evaluation_over_boolean_hypercube.par_iter()
                        .zip(witness_1.f_star_poly.evaluation_over_boolean_hypercube.par_iter())
                        .map(|(&e2, &e1)| e2 * two - e1)
                        .collect(),
                    len: witness_1.f_star_poly.len(),
                },
                tree_x: EqTree {
                    nodes: witness_2.tree_x.nodes.par_iter()
                        .zip(witness_1.tree_x.nodes.par_iter())
                        .map(|(&w2, &w1)| w2 * two - w1)
                        .collect(),
                    depth: witness_1.tree_x.depth,
                },
                tree_y: EqTree {
                    nodes: witness_2.tree_y.nodes.par_iter()
                        .zip(witness_1.tree_y.nodes.par_iter())
                        .map(|(&w2, &w1)| w2 * two - w1)
                        .collect(),
                    depth: witness_1.tree_y.depth,
                },
            },
            instance: AccInstance {
                // not used by helper function, so we simply pass them as zero or any other random element
                C: E::G1Affine::zero(),
                T: E::G1Affine::zero(),
                E: E::G1Affine::zero(),
                // used by the helper function
                x: instance_2.x.par_iter()
                    .zip(instance_1.x.par_iter())
                    .map(|(&e2, &e1)| e2 * two - e1)
                    .collect(),
                y: instance_2.y.par_iter()
                    .zip(instance_1.y.par_iter())
                    .map(|(&e2, &e1)| e2 * two - e1)
                    .collect(),
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

    pub fn decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> bool {
        let instance = &acc.instance;
        let witness = &acc.witness;

        // first condition
        let pairing_lhs = E::multi_pairing(&witness.vec_D, &srs.pc_srs.vec_V);
        let pairing_rhs = E::pairing(instance.C, srs.pc_srs.V_prime);

        // second condition
        let ip_rhs = instance.T;
        let ip_lhs = {
            let res = E::G1::msm_unchecked(srs.k_x.as_slice(), witness.tree_x.nodes.as_slice());
            res.add(E::G1::msm_unchecked(srs.k_y.as_slice(), witness.tree_y.nodes.as_slice()))
        };

        // third condition
        let verify_lhs = Self::helper_function_decide(srs, acc);
        let verify_rhs = instance.E;

        // return and of all three conditions
        (verify_rhs == verify_lhs.into()) && (ip_lhs == ip_rhs.into()) && (pairing_lhs == pairing_rhs)
    }

    pub fn helper_function_decide(srs: &AccSRS<E>, acc: &Accumulator<E>) -> E::G1Affine {
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
                res.extend_from_slice(srs.pc_srs.vec_H.as_slice());
                res.extend_from_slice(witness.vec_D.as_slice());

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

impl<E: Pairing> Accumulator<E> {
    /// this function returns a random satisfying accumulator by generating two random frseh accumualtors (KZH openings)
    /// and then accumulating them, so that the error vector wouldn't be zero
    pub fn rand<R: RngCore>(srs: &AccSRS<E>, rng: &mut R) -> Accumulator<E>
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
        let com1 = PCSEngine::commit(&srs.pc_srs, &polynomial1);
        let com2 = PCSEngine::commit(&srs.pc_srs, &polynomial2);

        // open the commitment
        let open1 = PCSEngine::open(&polynomial1, com1.clone(), &x1);
        let open2 = PCSEngine::open(&polynomial2, com2.clone(), &x2);

        // verify the proof
        PCSEngine::verify(&srs.pc_srs, &com1, &open1, &x1, &y1, &z1);
        PCSEngine::verify(&srs.pc_srs, &com2, &open2, &x2, &y2, &z2);

        let instance1 = Accumulator::new_accumulator_instance_from_fresh_kzh_instance(&srs, &com1.C, &x1, &y1, &z1);
        let witness1 = Accumulator::new_accumulator_witness_from_fresh_kzh_witness(&srs, open1, &x1, &y1);
        let instance2 = Accumulator::new_accumulator_instance_from_fresh_kzh_instance(&srs, &com2.C, &x2, &y2, &z2);
        let witness2 = Accumulator::new_accumulator_witness_from_fresh_kzh_witness(&srs, open2, &x2, &y2);

        let acc1 = Accumulator::new(&instance1, &witness1);
        let acc2 = Accumulator::new(&instance2, &witness2);

        // verify that the fresh accumulators are satisfied
        debug_assert!(Accumulator::decide(&srs, &acc1));
        debug_assert!(Accumulator::decide(&srs, &acc2));

        let mut prover_transcript = Transcript::new(b"new_transcript");

        let (accumulated_instance, accumulated_witness, Q) = Accumulator::prove(&srs, &acc1, &acc2, &mut prover_transcript);

        let accumulated_acc = Accumulator::new(&accumulated_instance, &accumulated_witness);

        // verify the accumulated instance is satisfied
        debug_assert!(Accumulator::decide(&srs, &accumulated_acc));

        accumulated_acc
    }
}

pub mod test {
    use ark_ec::pairing::Pairing;
    use rand::thread_rng;

    use super::*;
    use crate::constant_for_curves::{ScalarField, E};
    use crate::pcs::multilinear_pcs::{PCSEngine, PolynomialCommitmentSRS};

    #[test]
    fn test_accumulator_end_to_end() {
        let (degree_x, degree_y) = (128usize, 128usize);
        let srs_pcs: PolynomialCommitmentSRS<E> = PCSEngine::setup(degree_x, degree_y, &mut thread_rng());
        let srs = Accumulator::setup(srs_pcs.clone(), &mut thread_rng());

        let acc1 = Accumulator::rand(&srs, &mut thread_rng());
        let acc2 = Accumulator::rand(&srs, &mut thread_rng());

        let mut prover_transcript = Transcript::new(b"new_transcript");
        let mut verifier_transcript = prover_transcript.clone();

        let (instance, witness, Q) = Accumulator::prove(&srs, &acc1, &acc2, &mut prover_transcript);

        let instance_expected = Accumulator::verify(&srs, &acc1.instance, &acc2.instance, Q, &mut verifier_transcript);

        assert_eq!(instance, instance_expected);

        let decide_timer = start_timer!(|| "decide");
        assert!(Accumulator::decide(&srs, &Accumulator { witness, instance }));
        end_timer!(decide_timer);
    }

    #[test]
    fn test_accumulator_sizes() {
        // Hossein: change the degrees later, it takes too long
        let degrees = vec![(2, 2), (4, 4), (8, 8), (16, 16), (32, 32)];
        let rng = &mut thread_rng();

        for (degree_x, degree_y) in degrees {
            // set accumulator sts
            let srs_pcs: PolynomialCommitmentSRS<E> = PCSEngine::setup(degree_x, degree_y, rng);
            let srs_acc = Accumulator::setup(srs_pcs.clone(), rng);
            // get random accumulator
            let acc = Accumulator::rand(&srs_acc, rng);

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


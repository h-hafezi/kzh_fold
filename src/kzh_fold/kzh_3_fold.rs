use crate::gadgets::non_native::util::convert_affine_to_scalars;
use crate::kzh::kzh3::{KZH3Commitment, KZH3Opening, KZH3, KZH3SRS};
use crate::kzh::KZH;
use crate::kzh_fold::eq_tree::EqTree;
use crate::kzh_fold::{generate_random_elements, generic_linear_combination};
use crate::math::Math;
use crate::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use crate::transcript::transcript::Transcript;
use crate::utils::inner_product;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ec::{AffineRepr, VariableBaseMSM};
use ark_ff::{AdditiveGroup, Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use derivative::Derivative;
use rand::{thread_rng, RngCore};
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Acc3Error<E: Pairing> {
    E: E::G1Affine,
}

impl<E: Pairing> Acc3Error<E> {
    /// Returns a new `Acc3Proof` where all elements are `E::G1Affine::zero()`.
    pub fn zero() -> Self {
        Self {
            E: E::G1Affine::zero(),
        }
    }

    /// Returns the elements as a vector.
    pub fn to_vec(&self) -> Vec<E::G1Affine> {
        vec![self.E]
    }
}

impl<E: Pairing> Acc3Error<E> {
    /// Updates the `Acc3Proof` instance based on the provided error instances, `pf`, and `beta`.
    pub fn update(
        error_instance_1: &Self,
        error_instance_2: &Self,
        pf: &Self,
        beta: E::ScalarField,
    ) -> Self {
        let one_minus_beta = E::ScalarField::ONE - beta;
        let beta_one_minus_beta = beta * one_minus_beta;

        let E = error_instance_1.E.mul(one_minus_beta)
            + error_instance_2.E.mul(beta)
            + pf.E.mul(beta_one_minus_beta);


        Self { E: E.into() }
    }
}

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

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Acc3Instance<E: Pairing> {
    pub C: E::G1Affine,
    pub C_y: E::G1Affine,
    pub T: E::G1Affine,
    pub E: Acc3Error<E>,

    // vector of length log2(degree_x)
    pub x: Vec<E::ScalarField>,
    // vector of length log2(degree_y)
    pub y: Vec<E::ScalarField>,
    // vector of length log2(degree_z)
    pub z: Vec<E::ScalarField>,
    // result of polynomial evaluation
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

        for E in self.E.to_vec() {
            let (e_x, e_y) = convert_affine_to_scalars::<E>(E);
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
        instance_2: &Acc3Instance<E>,
        proof: &Acc3Error<E>,
    ) -> E::ScalarField {
        // add the instances to the transcript
        transcript.append_scalars(b"instance 1", instance_1.to_sponge_field_elements().as_slice());
        transcript.append_scalars(b"instance 2", instance_2.to_sponge_field_elements().as_slice());

        // convert the proof Q into scalar field elements and add to the transcript
        for E in proof.to_vec() {
            let (p1, p2) = convert_affine_to_scalars::<E>(E);
            transcript.append_scalars(b"proof", &[p1, p2]);
        }

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
            E: Acc3Error::zero(),
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
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
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

// verify and prove function
impl<E: Pairing> Accumulator3<E> {
    pub fn prove(
        srs: &Acc3SRS<E>,
        acc_1: &Accumulator3<E>,
        acc_2: &Accumulator3<E>,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> (Acc3Instance<E>, Acc3Witness<E>, Acc3Error<E>)
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        // unwrap the instances and witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        // compute the quotient variable Q
        let proof: Acc3Error<E> = Self::compute_error_term(srs, acc_1, acc_2);

        // the transcript is cloned because is used twice, once to compute the fiat_shamir challenge beta directly
        // once to call Self::verify() to get the new accumulated instance, despite that they generate the same randomness
        // and one can do fewer hashes by computing accumulated instance directly, it's for the sake of cleaner code.
        let mut transcript_clone = transcript.clone();

        // get challenge beta
        let beta = Accumulator3::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, &proof);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        // get the accumulated new_instance
        let new_instance = Self::verify(srs, instance_1, instance_2, &proof, &mut transcript_clone);

        // get the accumulated witness
        let new_witness = Acc3Witness {
            D_x: generic_linear_combination(
                &witness_1.D_x,
                &witness_2.D_x,
                |d_1, d_2| d_1.mul(one_minus_beta).add(d_2.mul(beta)),
            ),
            D_y: generic_linear_combination(
                &witness_1.D_y,
                &witness_2.D_y,
                |d_1, d_2| d_1.mul(one_minus_beta).add(d_2.mul(beta)),
            ),
            f_star: MultilinearPolynomial::linear_combination(
                &witness_1.f_star,
                &witness_2.f_star,
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
        };

        (new_instance, new_witness, proof)
    }

    pub fn verify(
        _srs: &Acc3SRS<E>,
        instance_1: &Acc3Instance<E>,
        instance_2: &Acc3Instance<E>,
        proof: &Acc3Error<E>,
        transcript: &mut Transcript<E::ScalarField>,
    ) -> Acc3Instance<E>
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        let beta = Accumulator3::compute_fiat_shamir_challenge(transcript, instance_1, instance_2, &proof);
        let one_minus_beta: E::ScalarField = E::ScalarField::ONE - beta;

        let new_error_term = Acc3Error::update(
            &instance_1.E,
            &instance_2.E,
            &proof,
            beta,
        );

        Acc3Instance {
            C: (instance_1.C.mul(one_minus_beta) + instance_2.C.mul(beta)).into(),
            C_y: (instance_1.C_y.mul(one_minus_beta) + instance_2.C_y.mul(beta)).into(),
            T: (instance_1.T.mul(one_minus_beta) + instance_2.T.mul(beta)).into(),
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
            E: new_error_term,
            output: instance_1.output * one_minus_beta + instance_2.output * beta,
        }
    }

    pub fn compute_error_term(
        srs: &Acc3SRS<E>,
        acc_1: &Accumulator3<E>,
        acc_2: &Accumulator3<E>,
    ) -> Acc3Error<E>
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        // unwrap the instances/witnesses
        let instance_1 = &acc_1.instance;
        let instance_2 = &acc_2.instance;
        let witness_1 = &acc_1.witness;
        let witness_2 = &acc_2.witness;

        assert_eq!(witness_1.f_star.num_variables, witness_2.f_star.num_variables);
        assert_eq!(witness_1.tree_x.depth, witness_2.tree_x.depth);
        assert_eq!(witness_1.tree_y.depth, witness_2.tree_y.depth);
        assert_eq!(witness_1.D_x.len(), witness_2.D_x.len());

        // value of 2 in the field
        let two = E::ScalarField::from(2u128);

        let witness = Acc3Witness {
            D_x: generic_linear_combination(
                &witness_2.D_x,
                &witness_1.D_x,
                |d2, d1| d2.mul(two).sub(d1),
            ),
            D_y: generic_linear_combination(
                &witness_2.D_y,
                &witness_1.D_y,
                |d2, d1| d2.mul(two).sub(d1),
            ),
            f_star: MultilinearPolynomial::linear_combination(
                &witness_2.f_star,
                &witness_1.f_star,
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
        };

        let instance = Acc3Instance {
            C_y: (instance_2.C_y + instance_2.C_y - instance_1.C_y).into(),
            // not used by helper function, so we simply pass them as zero or any other random element
            C: E::G1Affine::zero(),
            T: E::G1Affine::zero(),
            E: Acc3Error::zero(),
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
            output: two * instance_2.output - instance_1.output,
        };

        let acc = Accumulator3::new(&instance, &witness);

        let mut E = {
            let mut res = Self::dec(&srs, &acc)
                + &instance_1.E.E
                - &instance_2.E.E
                - &instance_2.E.E;

            // -1/2 in the scalar field
            let minus_one_over_two: E::ScalarField = two.neg().inverse().unwrap();
            res.mul(minus_one_over_two).into()
        };


        Acc3Error {
            E
        }
    }
}

// deciding functions
impl<E: Pairing> Accumulator3<E> {
    pub fn dec(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
        // Extract required components
        let error_tree_x = acc.witness.tree_x.difference(acc.instance.x.as_slice());
        let error_tree_y = acc.witness.tree_y.difference(acc.instance.y.as_slice());
        let error_tree_z = acc.witness.tree_z.difference(acc.instance.z.as_slice());

        // Compute e_prime for dec_2
        let e_prime: E::ScalarField = inner_product(
            &acc.witness.f_star.evaluation_over_boolean_hypercube,
            &acc.witness.tree_z.get_leaves(),
        ) - acc.instance.output;

        // Prepare scalars and bases for the combined MSM
        let mut scalars = Vec::new();
        let mut bases = Vec::new();

        // Add dec_1 components
        scalars.extend_from_slice(error_tree_x.nodes.as_slice());
        bases.extend_from_slice(srs.k_x.as_slice());

        scalars.extend_from_slice(error_tree_y.nodes.as_slice());
        bases.extend_from_slice(srs.k_y.as_slice());

        scalars.extend_from_slice(error_tree_z.nodes.as_slice());
        bases.extend_from_slice(srs.k_z.as_slice());

        // Add dec_2 component
        scalars.push(e_prime);
        bases.push(srs.k_prime);

        // Add dec_3 components
        scalars.extend_from_slice(
            acc.witness.f_star.evaluation_over_boolean_hypercube.as_slice(),
        );
        bases.extend_from_slice(srs.pc_srs.H_z.as_slice());

        scalars.extend_from_slice(
            acc.witness
                .tree_y
                .get_leaves()
                .iter()
                .map(|g| g.neg())
                .collect::<Vec<_>>()
                .as_slice(),
        );
        bases.extend(
            acc.witness
                .D_y
                .iter()
                .map(|g| g.clone().into())
                .collect::<Vec<_>>(),
        );

        // Add dec_4 components
        scalars.extend_from_slice(
            acc.witness
                .tree_x
                .get_leaves()
        );
        bases.extend(
            acc.witness
                .D_x
                .iter()
                .map(|g| g.clone().into())
                .collect::<Vec<_>>(),
        );

        // Perform the combined MSM
        let msm_result = E::G1::msm_unchecked(bases.as_slice(), scalars.as_slice());

        // Add final adjustments from dec_3 and dec_4
        let final_result = msm_result - acc.instance.C_y;

        // Return as affine point
        final_result.into()
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
            Self::dec(srs, acc),
            acc.instance.E.E.into(),
            "third condition fails"
        );

        // forth condition
        let pairing_lhs = E::multi_pairing(&witness.D_y, &srs.pc_srs.V_y);
        let pairing_rhs = E::pairing(instance.C_y, srs.pc_srs.v);
        assert_eq!(pairing_lhs, pairing_rhs, "forth condition fails");
    }
}


// get fresh satisfying accumulator (zero error terms)
impl<E: Pairing<ScalarField=F>, F: PrimeField + Absorb> Accumulator3<E> {
    // Helper function to create an accumulator from a random polynomial
    fn rand_fresh_accumulator(
        srs: &Acc3SRS<E>,
    ) -> Accumulator3<E>
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        let num_vars = srs.pc_srs.degree_x.log_2() + srs.pc_srs.degree_y.log_2() + srs.pc_srs.degree_z.log_2();

        let input: Vec<F> = (0..num_vars)
            .map(|_| F::rand(&mut thread_rng()))
            .collect();

        let polynomial = MultilinearPolynomial::rand(num_vars, &mut thread_rng());

        let output = polynomial.evaluate(input.as_slice());

        let commitment = KZH3::commit(&srs.pc_srs, &polynomial);

        let opening = KZH3::open(&srs.pc_srs, input.as_slice(), &commitment, &polynomial);

        // Convert proof to instance and witness
        let acc_instance = Accumulator3::proof_to_accumulator_instance(
            &srs,
            input.as_slice(),
            &output,
            &commitment,
            &opening,
        );

        let acc_witness = Accumulator3::proof_to_accumulator_witness(
            &srs,
            commitment,
            opening,
            input.as_slice(),
        );

        Accumulator3::new(&acc_instance, &acc_witness)
    }

    pub fn rand(srs: &Acc3SRS<E>) -> Accumulator3<E>
    where
        <E as Pairing>::ScalarField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb,
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
    {
        // Create two accumulators
        let acc1 = Accumulator3::rand_fresh_accumulator(srs);
        let acc2 = Accumulator3::rand_fresh_accumulator(srs);

        // Prove and decide final accumulator
        let (instance, witness, _proof) = Accumulator3::prove(
            srs,
            &acc1,
            &acc2,
            &mut Transcript::new(b"hi"),
        );

        let final_acc = Accumulator3::new(&instance, &witness);

        Accumulator3::decide(srs, &final_acc);

        final_acc
    }
}

#[cfg(test)]
mod test {
    use crate::constant_for_curves::E;
    use crate::kzh::kzh3::{KZH3, KZH3SRS};
    use crate::kzh::KZH;
    use crate::kzh_fold::kzh_3_fold::{Acc3SRS, Accumulator3};
    use crate::math::Math;
    use crate::transcript::transcript::Transcript;
    use crate::utils::inner_product;
    use ark_ec::pairing::Pairing;
    use ark_ec::CurveGroup;
    use ark_ec::VariableBaseMSM;
    use ark_ff::AdditiveGroup;
    use rand::thread_rng;
    use std::ops::{Add, Mul, Neg};

    fn dec_1<E: Pairing>(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> <E as Pairing>::G1Affine {
        let error_tree_x = acc.witness.tree_x.difference(acc.instance.x.as_slice());
        let error_tree_y = acc.witness.tree_y.difference(acc.instance.y.as_slice());
        let error_tree_z = acc.witness.tree_z.difference(acc.instance.z.as_slice());

        let mut res: E::G1 = E::G1::ZERO;
        res = res.add(E::G1::msm_unchecked(srs.k_x.as_slice(), error_tree_x.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_y.as_slice(), error_tree_y.nodes.as_slice()));
        res = res.add(E::G1::msm_unchecked(srs.k_z.as_slice(), error_tree_z.nodes.as_slice()));

        res.into()
    }

    fn dec_2<E: Pairing>(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
        let instance = &acc.instance;
        let witness = &acc.witness;

        let e_prime: E::ScalarField = inner_product(
            &witness.f_star.evaluation_over_boolean_hypercube,
            &acc.witness.tree_z.get_leaves(),
        ) - instance.output;

        srs.k_prime.mul(e_prime).into()
    }

    fn dec_3<E: Pairing>(srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
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

    fn dec_4<E: Pairing>(_srs: &Acc3SRS<E>, acc: &Accumulator3<E>) -> E::G1Affine {
        let lhs = E::G1::msm_unchecked(
            acc.witness.D_x.iter().map(|g| g.clone().into()).collect::<Vec<_>>().as_slice(),
            acc.witness.tree_x.get_leaves(),
        );

        acc.instance.C_y.add(lhs.neg()).neg().into()
    }

    #[test]
    fn test_dec() {
        let (degree_x, degree_y, degree_z) = (4usize, 2usize, 8usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        // build the srs
        let pcs_srs: KZH3SRS<E> = KZH3::setup((degree_x * degree_y * degree_z).log_2(), &mut thread_rng());
        let acc_srs = Accumulator3::setup(pcs_srs, &mut thread_rng());

        let acc = Accumulator3::rand(&acc_srs);

        assert_eq!(
            Accumulator3::dec(&acc_srs, &acc),
            (dec_1(&acc_srs, &acc)
                + dec_2(&acc_srs, &acc)
                + dec_3(&acc_srs, &acc)
                + dec_4(&acc_srs, &acc)).into_affine(),
        )
    }

    #[test]
    fn kzh3_fold_test() {
        let (degree_x, degree_y, degree_z) = (4usize, 2usize, 8usize);
        let num_vars = degree_x.log_2() + degree_y.log_2() + degree_z.log_2();

        // build the srs
        let pcs_srs: KZH3SRS<E> = KZH3::setup((degree_x * degree_y * degree_z).log_2(), &mut thread_rng());
        let acc_srs = Accumulator3::setup(pcs_srs, &mut thread_rng());

        let acc_1 = Accumulator3::rand(&acc_srs);
        Accumulator3::decide(&acc_srs, &acc_1);
        let acc_2 = Accumulator3::rand(&acc_srs);
        Accumulator3::decide(&acc_srs, &acc_2);

        let (instance, witness, proof) = Accumulator3::prove(&acc_srs, &acc_1, &acc_2, &mut Transcript::new(b"hi"));

        Accumulator3::decide(&acc_srs, &Accumulator3::new(&instance, &witness));
    }
}

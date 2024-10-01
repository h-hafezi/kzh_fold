use std::iter;

use ark_crypto_primitives::sponge::Absorb;
use ark_ec::AffineRepr;
use ark_ff::{PrimeField, UniformRand};
use rand::RngCore;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use transcript::IOPTranscript;

use crate::accumulation::accumulator::{AccInstance, AccWitness, Accumulator};
use ark_ff::Zero;
use crate::polynomial::multilinear_poly::MultilinearPolynomial;
use crate::spartan::sumcheck::SumcheckInstanceProof;
use crate::{accumulation, pcs};
use crate::pcs::multilinear_pcs::{OpeningProof, PolyCommit, PolyCommitTrait, Commitment, SRS as PcsSRS};

// XXX move to mod.rs or somewhere neutral
#[derive(Clone, Debug)]
pub struct SRS<E: Pairing> {
    pub acc_srs: accumulation::accumulator::AccSRS<E>,
}

impl<E: Pairing> SRS<E>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    <E as Pairing>::ScalarField: Absorb,
{
    fn new<T: RngCore>(degree_x: usize, degree_y: usize, rng: &mut T) -> Self {
        let pcs_srs =  PolyCommit::setup(degree_x, degree_y, rng);

        SRS {
            acc_srs: Accumulator::setup(pcs_srs, rng),
        }
    }
}

pub struct SignatureAggrData<E: Pairing> {
    // TODO comment this out for now. we will figure out the BLS stuff later.
    //pk: E::G1Affine,
    //sig: E::G2Affine,
    bitfield_poly: MultilinearPolynomial<E::ScalarField>,
    bitfield_commitment: Commitment<E>,
    sumcheck_proof: Option<SumcheckInstanceProof<E::ScalarField>>,
    // TODO Hossein: For now, instead of a proof, let's just put the R1CS circuit here
    // ivc_proof: IVCProof<E>
}

impl<E: Pairing> SignatureAggrData<E> {
    pub fn new(bitfield_poly: MultilinearPolynomial<E::ScalarField>, _sumcheck_proof: Option<SumcheckInstanceProof<E::ScalarField>>, srs: &SRS<E>) -> Self {
        // XXX this PolyCommit is not very ergonomic
        let poly_commit = PolyCommit { srs: srs.acc_srs.pc_srs.clone() }; // XXX no clone
        let bitfield_commitment = poly_commit.commit(&bitfield_poly);
        SignatureAggrData {
            bitfield_poly,
            bitfield_commitment,
            sumcheck_proof: None
        }
    }
}

/// This struct represents an aggregator on the network that receives data from two parties and needs to aggregate them
/// into one. After aggregation, the aggregator forwards the aggregated data.
pub struct Aggregator<E: Pairing> {
    srs: SRS<E>,
    A_1: SignatureAggrData<E>,
    A_2: SignatureAggrData<E>,
}

impl<E: Pairing> Aggregator<E>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    <E as Pairing>::ScalarField: Absorb,
{
    fn get_accumulator_from_evaluation(&self,
                                       bitfield_poly: &MultilinearPolynomial<E::ScalarField>,
                                       bitfield_commitment: &Commitment<E>,
                                       alpha: &Vec<E::ScalarField>,
                                       beta: &Vec<E::ScalarField>) -> Accumulator<E> {
        let poly_commit = PolyCommit { srs: self.srs.acc_srs.pc_srs.clone() }; // XXX no clone. bad ergonomics

        // Get (alpha, beta) as the evaluation point
        let eval_point: Vec<E::ScalarField> = alpha.clone().into_iter().chain(beta.clone()).collect(); // XXX bad clone

        let y = bitfield_poly.evaluate(&eval_point);
        let opening_proof = poly_commit.open(bitfield_poly, bitfield_commitment.clone(), &alpha); // XXX needless clone

        // XXX bad name for function
        let acc_instance = Accumulator::new_accumulator_instance_from_proof(
            &self.srs.acc_srs,
            &bitfield_commitment.C,
            alpha,
            beta,
            &y
        );
        let acc_witness = Accumulator::new_accumulator_witness_from_proof(
            &self.srs.acc_srs,
            opening_proof,
            alpha,
            beta,
        );
        Accumulator {
            witness: acc_witness,
            instance: acc_instance,
        }
    }

    pub fn aggregate(&self, _transcript: &mut IOPTranscript<E::ScalarField>) -> SignatureAggrData<E> {
        let poly_commit = PolyCommit { srs: self.srs.acc_srs.pc_srs.clone() }; // XXX no clone. bad ergonomics
        // let pk = self.A_1.pk + self.A_2.pk;
        // let sk = self.A_1.sig + self.A_2.sig;

        let c_poly = self.A_1.bitfield_poly.get_bitfield_union_poly(&self.A_2.bitfield_poly);
        let C_commitment = poly_commit.commit(&c_poly);

        // We do sumcheck for the following polynomial:
        // eq(r,x) * (b_1 + b_2 - b_1 * b_2 - c)
        let union_comb_func =
            |poly_eq: &E::ScalarField, poly_b_1: &E::ScalarField, poly_b_2: &E::ScalarField, poly_c: &E::ScalarField|
                                              -> E::ScalarField { *poly_eq * (*poly_b_1 + *poly_b_2 - *poly_b_1 * *poly_b_2 - *poly_c) };

        let (sumcheck_proof, (alpha, beta)) = SumcheckInstanceProof::prove_cubic(&E::ScalarField::zero());
                                                          

        // XXX remove
        let alpha: Vec<E::ScalarField> = iter::repeat_with(|| E::ScalarField::zero()).take(3).collect();
        let beta: Vec<E::ScalarField> = iter::repeat_with(|| E::ScalarField::zero()).take(3).collect();

        // Now the verifier will need:
        // y_1 = b_1(alpha, beta)
        // y_2 = b_2(alpha, beta), and
        // y_3 = c(alpha, beta)
        // to verify the sumcheck
        // Compute the evaluations and its accumulations
        let y_1_accumulator = self.get_accumulator_from_evaluation(
            &self.A_1.bitfield_poly,
            &self.A_1.bitfield_commitment,
            &alpha,
            &beta,
        );
        let y_2_accumulator = self.get_accumulator_from_evaluation(
            &self.A_2.bitfield_poly,
            &self.A_2.bitfield_commitment,
            &alpha,
            &beta,
        );
        let _y_3_accumulator = self.get_accumulator_from_evaluation(
            &self.A_2.bitfield_poly,
            &self.A_2.bitfield_commitment,
            &alpha,
            &beta,
        );

        // Here we need to accumulate y_1 acc, y_2 acc, and y_3 acc into one.
        // TODO Hossein: let's just do y_1 with y_2 for now. but we will need a tree for later.
        let _acc_prime = Accumulator::prove(&self.srs.acc_srs, &y_1_accumulator, &y_2_accumulator);

        // TODO Hossein: Now we want an IVC proof of the accumulation
        // let ivc_proof = accumulation_circuit::prove_accumulation(&acc_prime, &y_1_accumulator, &y_2_accumulator, &self.srs.acc_srs);

        SignatureAggrData {
            bitfield_poly: c_poly,
            bitfield_commitment: C_commitment,
            sumcheck_proof: None,
            // ivc_proof: ivc_proof
        }
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use ark_poly::{EvaluationDomain,GeneralEvaluationDomain};
    use ark_std::test_rng;
    use ark_std::UniformRand;
    use crate::constant_for_curves::{E, ScalarField};

    #[test]
    fn test_aggregate() {
        let rng = &mut rand::thread_rng();
        let mut transcript = IOPTranscript::<ScalarField>::new(b"aggr");

        // num_vars = log(degree_x) + log(degree_y)
        let degree_x = 8usize;
        let degree_y = 8usize;
        let num_vars = 6usize;

        let srs = SRS::<E>::new(degree_x, degree_y, rng);

        let b_1 = MultilinearPolynomial::random_binary(num_vars, rng);
        let sig_aggr_data_1 = SignatureAggrData::new(b_1, None, &srs);

        let b_2 = MultilinearPolynomial::random_binary(num_vars, rng);
        let sig_aggr_data_2 = SignatureAggrData::new(b_2, None, &srs);

        let aggregator = Aggregator {
            srs: srs.clone(),
            A_1: sig_aggr_data_1,
            A_2: sig_aggr_data_2,
        };

        let _agg_data = aggregator.aggregate(&mut transcript);
        // TODO Hossein: Print the constraint count of the R1CS circuit

        // TODO Hossein: Check that the witness satisfies the witness and examine the witness for 1s and 0s

        // Now let's do verification
        // let verifier = Verifier {
        //     srs,
        //     A: agg_data
        // };

        // assert_eq!(true, verifier.verify())
    }
}



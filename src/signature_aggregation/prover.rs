use ark_ec::AffineRepr;
use ark_ff::UniformRand;
use rand::RngCore;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use transcript::IOPTranscript;

use crate::accumulation::accumulator::get_new_acc_srs;
use crate::{accumulation, polynomial_commitment};
use crate::signature_aggregation::bivariate_sumcheck;
use crate::signature_aggregation::bivariate_sumcheck::SumcheckProof;
use crate::{polynomial::bivariate_poly::BivariatePolynomial};

use crate::polynomial_commitment::pcs::{Commitment, OpeningProof, PolyCommit};
use crate::polynomial_commitment::pcs;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SRS<E: Pairing> {
    pub pcs_srs: polynomial_commitment::pcs::SRS<E>,
    pub acc_srs: accumulation::accumulator::AccSRS<E>,
}

impl<E: Pairing> SRS<E> {
    fn new<T: RngCore>(degree_x: usize, degree_y: usize, rng: &mut T) -> Self {
        SRS {
            acc_srs: get_new_acc_srs(degree_x, degree_y),
            pcs_srs: PolyCommit::setup(degree_x, degree_y, rng)
        }
    }
}

pub struct SignatureAggrData<E: Pairing> {
    pk: E::G1Affine,
    sig: E::G2Affine,
    bitfield_commitment: Commitment<E>,
    bitfield_poly: BivariatePolynomial<E::ScalarField>,
    sumcheck_proof: SumcheckProof<E>,
    // Also need SNARK proof for IVC verifier
}

/// This struct represents an aggregator on the network that receives data from two parties and needs to aggregate them
/// into one. After aggregation, the aggregator forwards the aggregated data.
pub struct Aggregator<E: Pairing> {
    srs: SRS<E>,
    A_1: SignatureAggrData<E>,
    A_2: SignatureAggrData<E>,
}

impl<E: Pairing> Aggregator<E> {
    #[allow(unused_variables)] // XXX remove
    pub fn aggregate(&self, transcript: &mut IOPTranscript<E::ScalarField>) -> () { // SignatureAggrData<E> {
        let pk = self.A_1.pk + self.A_2.pk;
        let sk = self.A_1.sig + self.A_2.sig;

        let c_poly = self.A_1.bitfield_poly.bitfield_union(&self.A_2.bitfield_poly);
        // let C_commitment = CoeffFormPCS::commit(&c_poly, &self.srs);

        // Now aggregate all three polys into one
        // XXX
        // let f_poly = b_1 + b_2 - b_1*b_2 - c;
        // for now let's pretend it's c_poly

        let f_poly = c_poly.clone();
        let (sumcheck_proof, (alpha, beta)) = bivariate_sumcheck::prove::<E>(&f_poly, transcript);

        

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
        let mut _transcript = IOPTranscript::<ScalarField>::new(b"aggr");

        let degree_x = 16usize;
        let degree_y = 4usize;

        let srs = SRS::new(degree_x, degree_y, rng);

        // XXX Create two valid signature aggr data and aggregate
        let domain_x = GeneralEvaluationDomain::<ScalarField>::new(degree_x).unwrap();
        let domain_y = GeneralEvaluationDomain::<ScalarField>::new(degree_y).unwrap();
        let _b_1 = BivariatePolynomial::random_binary(rng, domain_x, domain_y, degree_x, degree_y);
        let _b_2 = BivariatePolynomial::random_binary(rng, domain_x, domain_y, degree_x, degree_y);
    }
}


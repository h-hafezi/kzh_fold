use ark_ec::AffineRepr;
use ark_ff::UniformRand;
use rand::RngCore;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use transcript::IOPTranscript;

use crate::signature_aggregation::bivariate_sumcheck::{bivariate_sumcheck, SumcheckProof};
use crate::{polynomial::bivariate_poly::BivariatePolynomial};

use crate::pcs::{Commitment, OpeningProof, CoeffFormPCS, SRS};

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
        let C_commitment = CoeffFormPCS::commit(&c_poly, &self.srs);

        // Now aggregate all three polys into one
        // let f_poly = b_1 + b_2 - b_1*b_2 - c

        // for now let's pretend it's c_poly
        let f_poly = c_poly.clone();
        let sumcheck_proof: SumcheckProof<E> = bivariate_sumcheck(transcript, &f_poly);

        unimplemented!()
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use ark_std::test_rng;
    use ark_std::UniformRand;
    use ark_bn254::{Fr, G1Projective};

    #[test]
    fn test_aggregate() {
        let rng = &mut rand::thread_rng();
        let mut _transcript = IOPTranscript::<Fr>::new(b"aggr");

        // XXX Create two valid signature aggr data and aggregate
        let _b_1 = BivariatePolynomial::<Fr>::random_binary(rng, 4);
        let _b_2 = BivariatePolynomial::<Fr>::random_binary(rng, 4);
    }
}


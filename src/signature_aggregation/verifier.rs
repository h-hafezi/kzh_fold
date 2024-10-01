/*use ark_crypto_primitives::sponge::Absorb;
use ark_ec::AffineRepr;
use ark_ff::{PrimeField, UniformRand};
use rand::RngCore;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use transcript::IOPTranscript;

use crate::accumulation::accumulator::{get_srs, AccInstance, AccWitness, Accumulator};
use crate::{accumulation, pcs};
use crate::signature_aggregation::bivariate_sumcheck;
use crate::signature_aggregation::bivariate_sumcheck::SumcheckProof;
use crate::signature_aggregation::prover::{SRS, SignatureAggrData};

use crate::pcs::bivariate_pcs::{Commitment, OpeningProof, PolyCommit, PolyCommitTrait};
use crate::pcs::bivariate_pcs;

/// This struct represents a network node that just received an aggregate signature. The verifier needs to verify the
/// aggregate signature (and later aggregate it with more signatures herself).
/// For the purposes of this module, we will only do the verification.
pub struct Verifier<E: Pairing> {
    pub srs: SRS<E>,
    pub A: SignatureAggrData<E>,
}

impl<E: Pairing> Verifier<E>
where
    <<E as Pairing>::G1Affine as AffineRepr>::BaseField: Absorb + PrimeField,
    <E as Pairing>::ScalarField: Absorb,
{
    pub fn verify(&self) -> bool {
        // TODO Verify the sumcheck proof (need to build the prover first)

        // Verify the IVC proof

        // Verify the BLS signature

        // At some point, run the decider

        true
    }
}



 */

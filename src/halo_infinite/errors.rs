use thiserror::Error;

/// Represents an error in proof creation, verification, or parsing.
#[derive(Error, Clone, Debug, Eq, PartialEq)]
pub enum ProofError {
    /// This error occurs when a proof failed to verify.
    #[error("Proof verification failed.")]
    VerificationError,
}
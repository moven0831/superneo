//! Error types for the compression layer.

use thiserror::Error;

/// Errors in the final SNARK or its polynomial commitment scheme.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum SnarkError {
    /// The SNARK verifier rejected the proof.
    #[error("SNARK verification failed: {0}")]
    Verify(String),

    /// A polynomial-commitment opening (FRI proximity / Merkle path) was invalid.
    #[error("PCS opening invalid: {0}")]
    PcsOpening(String),

    /// A structural / dimension inconsistency in the final relation.
    #[error("malformed final relation: {0}")]
    Malformed(String),
}

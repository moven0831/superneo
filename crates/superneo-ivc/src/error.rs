//! Error types for the IVC layer.

use thiserror::Error;

/// Errors in the native IVC loop or the recursive verifier circuit.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum IvcError {
    /// A fold step's verifier check failed.
    #[error("fold step {step} failed verification")]
    StepVerify { step: usize },

    /// The Construction-2 IO digest binding two iterations did not match.
    #[error("IVC digest mismatch at step {step}")]
    DigestMismatch { step: usize },

    /// The recursive verifier circuit was unsatisfied (P2).
    #[error("recursive verifier circuit unsatisfied: {0}")]
    CircuitUnsat(String),
}

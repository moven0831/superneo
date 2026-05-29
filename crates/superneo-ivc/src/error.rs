//! Error types for the IVC layer.

use superneo_fold::FoldError;
use thiserror::Error;

/// Errors in the native IVC loop or the recursive verifier circuit.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum IvcError {
    /// A fold step's reduction failed verification.
    #[error("fold step {step} failed: {source}")]
    Step { step: usize, source: FoldError },

    /// The Construction-2 IO digest binding the iterations did not match.
    #[error("IVC digest mismatch")]
    DigestMismatch,

    /// The recomputed accumulator did not match the proof's claimed accumulator.
    #[error("final accumulator mismatch")]
    AccumulatorMismatch,

    /// The recursive verifier circuit was unsatisfied (P2).
    #[error("recursive verifier circuit unsatisfied: {0}")]
    CircuitUnsat(String),
}

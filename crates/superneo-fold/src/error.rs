//! Error types for the folding layer.

use thiserror::Error;

/// Errors in sum-check, the interactive reductions, or verification.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum FoldError {
    /// The sum-check verifier's final evaluation check failed (Definition 6).
    #[error("sum-check verification failed at round {round}")]
    SumCheck { round: usize },

    /// A reduction's verifier check failed (e.g. Π_DEC commitment additivity
    /// `c = Σ bⁱ⁻¹ cᵢ`, or a folded evaluation mismatch).
    #[error("reduction verification failed: {0}")]
    Reduction(String),

    /// A structural / dimension inconsistency in the relation instance.
    #[error("malformed relation instance: {0}")]
    Malformed(String),
}

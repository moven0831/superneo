//! Error types for the commitment layer.

use thiserror::Error;

/// Errors in commitment, parameter validation, or challenge sampling.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum CommitError {
    /// The parameter guard `(K+k)·T·(b−1) < B` (Definition 14 / Appendix D.8) failed.
    #[error("parameter guard violated: (K+k)·T·(b−1) = {lhs} >= B = {rhs}")]
    GuardFailed { lhs: u64, rhs: u64 },

    /// A structural dimension mismatch (e.g. `B != b^k`, or a vector whose length
    /// disagrees with the public matrix `A`).
    #[error("dimension mismatch: {0}")]
    DimMismatch(String),

    /// A committed vector exceeded the ℓ∞ norm bound the commitment binds at.
    #[error("commitment norm bound violated: ‖z‖∞ not < {bound}")]
    NormBound { bound: u64 },
}

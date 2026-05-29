//! Error types for the ring layer.

use thiserror::Error;

/// Errors that can arise in ring arithmetic, embedding, or decomposition.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum RingError {
    /// A vector length was not a multiple of `d = 54` where the SuperNeo embedding
    /// requires whole `d`-blocks (Definition 7).
    #[error("vector length {0} is not a multiple of the ring degree d=54")]
    Misaligned(usize),

    /// A decomposition or reconstruction violated the norm bound (Definition 3).
    #[error("norm bound violated: ‖·‖∞ = {got} not < {bound}")]
    NormBound { got: u64, bound: u64 },
}

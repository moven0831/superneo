//! Error types for the field layer.

use thiserror::Error;

/// Errors that can arise in `F_q` / `F_{q^2}` arithmetic.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum FieldError {
    /// A raw `u64` was supplied that is not a canonical residue (`>= q`) where a
    /// canonical value was required.
    #[error("value {0:#x} is not canonical (>= q)")]
    NonCanonical(u64),

    /// Attempted to invert the zero element (Definition 1: only nonzero elements
    /// are invertible).
    #[error("attempted to invert zero")]
    InverseOfZero,
}

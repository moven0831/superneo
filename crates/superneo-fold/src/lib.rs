//! `superneo-fold` — the folding protocol core.
//!
//! Implements (clean-room):
//!   * a Blake3-backed Fiat–Shamir [`Transcript`] (the protocol is transcript-backend
//!     agnostic),
//!   * the norm-bounded CCS and CCS-evaluation relations (Definitions 11–13),
//!   * the sum-check protocol over the extension field `K` (Definition 6),
//!   * the three interactive reductions `Π_CCS` (strong, Lemma 3), `Π_RLC` (weak,
//!     Lemma 4), `Π_DEC` (reduction of knowledge, Theorem 9), and their composition
//!     `Π_SuperNeo = Π_DEC ∘ Π_RLC ∘ Π_CCS` (Theorem 8) exposed as `fold` / `verify_fold`.
//!
//! Modules are added in milestone M4.

pub mod error;

pub use error::FoldError;

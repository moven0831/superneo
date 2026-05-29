//! `superneo-commit` — the Ajtai / Module-SIS commitment layer.
//!
//! Implements (clean-room):
//!   * the homomorphic Ajtai commitment `Commit(A, z) = A·z` (Definitions 4 & 18),
//!     with an efficient rot-and-add (pay-per-bit) inner product,
//!   * the Goldilocks parameter set of Appendix B.2 with a `validate()` guard
//!     mirroring the paper's Sage relations (`(K+k)·T·(b−1) < B`),
//!   * the strong sampling set `C` (Definition 17) and a sampler.
//!
//! Modules are added in milestone M3.

pub mod error;

pub use error::CommitError;

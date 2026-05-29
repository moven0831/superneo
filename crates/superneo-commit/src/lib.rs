//! `superneo-commit` — the Ajtai / Module-SIS commitment layer.
//!
//! Implements (clean-room):
//!   * the homomorphic Ajtai commitment `Commit(A, z) = A·z` ([`pp`], [`commit`];
//!     Definitions 4 & 18), an `R_F`-module homomorphism,
//!   * the Goldilocks parameter set of Appendix B.2 with a `validate()` guard
//!     mirroring the Sage relation `(K+k)·T·(b−1) < B` ([`params`]),
//!   * the strong sampling set `C` and its sampler ([`challenge`]; Definition 17).
//!
//! Security note: binding rests on the hardness of Module-SIS for these parameters;
//! that assumption is not (and cannot be) unit-tested — only the linear/structural
//! properties are.

pub mod challenge;
pub mod commit;
pub mod error;
pub mod params;
pub mod pp;

pub use commit::Commitment;
pub use error::CommitError;
pub use params::Params;
pub use pp::PublicParams;

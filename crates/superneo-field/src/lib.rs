//! `superneo-field` — the arithmetic foundation of the SuperNeo PoC.
//!
//! Implements, from scratch (clean-room; see the project plan):
//!   * [`Fp`](fp::Fp): the Goldilocks prime field `F_q` with `q = 2^64 − 2^32 + 1`,
//!     using Solinas reduction (`2^64 ≡ 2^32 − 1 (mod q)`).
//!   * [`Ext2`](ext2::Ext2): the degree-2 extension `K = F_{q^2} = F_q[u]/(u² − 7)`,
//!     over which the sum-check protocol attains negligible soundness error.
//!
//! Paper references: Definition 1 (Fields), Appendix B.2 (Goldilocks parameters).
//!
//! This crate has zero internal dependencies — it sits at the bottom of the
//! workspace DAG.

pub mod error;
pub mod ext2;
pub mod fp;

pub use error::FieldError;

/// The Goldilocks prime modulus `q = 2^64 − 2^32 + 1`.
pub const Q: u64 = 0xFFFF_FFFF_0000_0001;

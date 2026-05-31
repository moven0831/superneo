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
//! Built in milestone M4.

pub mod error;
pub mod fold;
pub mod multilinear;
pub mod pi_ccs;
pub mod pi_dec;
pub mod pi_rlc;
pub mod sumcheck;
pub mod transcript;
pub mod types;

pub use error::FoldError;
pub use fold::{fold, verify_fold};
pub use pi_ccs::{pi_ccs_prove, pi_ccs_verify, pi_ccs_verify_traced, PiCcsTrace};
pub use pi_dec::{pi_dec_prove, pi_dec_verify};
pub use pi_rlc::{pi_rlc_prove, pi_rlc_verify};
pub use sumcheck::{sumcheck_prove, sumcheck_verify, RoundPoly, SumCheckProof};
pub use transcript::Transcript;
pub use types::{
    CcsInstance, CcsStructure, CcsWitness, CeInstance, CeWitness, FoldProof, GlobalParams,
    PiCcsProof, PiDecProof, PiRlcProof, SparsePoly,
};

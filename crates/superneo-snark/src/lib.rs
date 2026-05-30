//! `superneo-snark` — final compression of the SuperNeo accumulator.
//!
//! Implements (clean-room) a Spartan-style sum-check SNARK that compresses the final
//! CCS-evaluation accumulator — its constant-term evaluation claims `ct(y_j) = ~(M_j z)(r)`
//! (Theorem 6) and the ℓ∞ norm bound `‖z‖∞ < b` — into a proof whose size is independent
//! of the IVC length (constant in the number of folds), with the witness committed by a
//! BaseFold/FRI multilinear polynomial commitment over Goldilocks ([79], [84]). The PCS
//! is post-quantum (hash-based, via Blake3 Merkle trees). At these PoC parameters the
//! proof is dominated by the FRI query phase and is not smaller than the witness.
//!
//! Layering:
//!   * [`code`] — the foldable Reed–Solomon code (roots of unity, NTT encoding, the
//!     multilinear basis transforms, and the codeword/coefficient fold commutation).
//!   * [`merkle`] — Blake3 Merkle commitments for the codeword layers.
//!   * [`mle`] — little-endian multilinear helpers shared by the PCS and the reduction.
//!   * [`basefold`] — `commit` / `open` / `verify` for the multilinear PCS.
//!   * [`spartan`] — the linear-claim batching + norm sum-check reducing the accumulator
//!     relation to witness-MLE openings.
//!   * [`compress`] / [`verify`] — the top-level API over an IVC accumulator.
//!
//! Built in milestone M6.

pub mod basefold;
pub mod code;
pub mod error;
pub mod merkle;
pub mod mle;
pub mod pipeline;
pub mod spartan;

pub use error::SnarkError;
pub use pipeline::{run_pipeline, PipelineReport};
pub use spartan::{compress, verify, SnarkProof};

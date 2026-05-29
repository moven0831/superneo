//! `superneo-snark` — final compression of the SuperNeo accumulator.
//!
//! Implements (clean-room) a Spartan-style sum-check SNARK that proves
//! satisfiability of the final CCS-evaluation accumulator — the Ajtai opening
//! `c = A·z`, the ℓ∞ norm bound, and the multilinear evaluations `yⱼ = ~(M̄ⱼz)(r)` —
//! compressed with a BaseFold/FRI multilinear polynomial commitment over Goldilocks
//! ([79], [84]). The `MlPcs` trait keeps the lattice-native Hachi PCS [69] as a
//! drop-in alternative.
//!
//! Built in milestone M6.

pub mod error;

pub use error::SnarkError;

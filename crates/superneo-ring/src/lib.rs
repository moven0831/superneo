//! `superneo-ring` — the cyclotomic ring layer.
//!
//! Implements (clean-room) the ring `R_F = F[X]/(Φ)` with `Φ = X^54 + X^27 + 1`
//! (the 81st cyclotomic polynomial, `d = 54`) and its extension `R_K = K[X]/(Φ)`,
//! plus the SuperNeo embedding machinery:
//!   * coefficient maps `cf`, `cf_inv`, `ct` (Definition 2),
//!   * the rotation / S-action matrices,
//!   * the **bar / inner-product transform** (Theorem 5) and its lifts
//!     (Definition 8, Theorems 6 & 7) — the technical heart of the scheme,
//!   * the balanced ℓ∞ norm (Definition 3) and base-`b` decomposition `split_b`.
//!
//! Modules are added milestone-by-milestone (M2). For now this crate exposes its
//! error type so the workspace builds end-to-end.

pub mod error;

pub use error::RingError;

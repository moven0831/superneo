//! `superneo-ring` — the cyclotomic ring layer.
//!
//! Implements (clean-room) the ring `R_F = F[X]/(Φ)` with `Φ = X^54 + X^27 + 1`
//! (the 81st cyclotomic polynomial, `d = 54`) plus the SuperNeo embedding machinery:
//!   * coefficient maps `cf`, `cf_inv`, `ct` ([`maps`], Definition 2),
//!   * the rotation / S-action matrix ([`s_action`]),
//!   * the **bar / inner-product transform** ([`bar`], Theorem 5) and its lifts
//!     (Definition 8, Theorem 6) — the technical heart of the scheme,
//!   * the balanced ℓ∞ norm ([`norm`], Definition 3) and base-`b` decomposition
//!     `split_b` ([`decomp`]).

pub mod bar;
pub mod decomp;
pub mod error;
pub mod maps;
pub mod norm;
pub mod ring;
pub mod ring_k;
pub mod s_action;

pub use error::RingError;
pub use ring::{RingElem, D};
pub use ring_k::RingK;

//! `superneo-ivc` — incrementally verifiable computation on top of the folding scheme.
//!
//! Implements (clean-room), in two phases:
//!   * **P1 (native, M5):** an IVC state machine that folds a sequence of step
//!     instances natively (`prove_native` / `verify_native`), with the
//!     Construction-2 IO encoding `enc_inst(hash(...))` binding each iteration to
//!     the next.
//!   * **P2 (recursive, M7):** the SuperNeo folding *verifier* expressed as CCS
//!     gadgets (K-field, sum-check, RLC, DEC, transcript), closing a true recursive
//!     IVC loop (`prove_ivc` / `verify_ivc`).
//!
//! Compiler blueprint: HyperNova / Nova / PCD ([56], [58], [86]).

pub mod error;

pub use error::IvcError;

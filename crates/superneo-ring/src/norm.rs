//! The balanced ℓ∞ norm (Definition 3).
//!
//! For a field element, `‖a‖∞ = |balanced(a)|`. For a vector / ring element /
//! ring vector, it is the maximum over coordinates / coefficients.

use superneo_field::fp::Fp;

use crate::ring::RingElem;

/// `‖·‖∞` of a field vector.
pub fn norm_inf_fp(v: &[Fp]) -> u64 {
    v.iter().map(|x| x.abs_balanced()).max().unwrap_or(0)
}

/// `‖·‖∞` of a single ring element (max over its coefficients).
pub fn norm_inf_ring(a: &RingElem) -> u64 {
    norm_inf_fp(a.coeffs())
}

/// `‖·‖∞` of a ring vector (max over all elements).
pub fn norm_inf_ring_vec(v: &[RingElem]) -> u64 {
    v.iter().map(norm_inf_ring).max().unwrap_or(0)
}

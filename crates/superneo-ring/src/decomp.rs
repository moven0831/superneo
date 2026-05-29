//! Base-`b` balanced decomposition `split_b` (Definition 3).
//!
//! Given a field vector `z` whose balanced coordinates satisfy `|v| < b^k`, produce
//! `k` vectors `z_1, …, z_k` with
//!
//! ```text
//!   z = Σ_{i=1}^k b^{i−1} · z_i      and      ‖z_i‖∞ < b.
//! ```
//!
//! Each coordinate is decomposed by sign-magnitude base-`b` expansion: the digits of
//! `|v|` lie in `{0, …, b−1}` and inherit the sign of `v`, so every digit has
//! magnitude `≤ b−1 < b`, and the (signed) reconstruction is exact over `F`.

use superneo_field::fp::Fp;

use crate::error::RingError;

/// Decompose `z` into `k` low-norm vectors in base `b` (Definition 3).
///
/// Returns `[z_1, …, z_k]` where `z_1` carries weight `b^0`. Errors with
/// [`RingError::NormBound`] if some coordinate's magnitude does not fit in `k`
/// base-`b` digits (i.e. `|balanced| ≥ b^k`).
pub fn split_b(z: &[Fp], b: u64, k: usize) -> Result<Vec<Vec<Fp>>, RingError> {
    assert!(b >= 2, "decomposition base must be ≥ 2");
    let bound = b
        .checked_pow(k as u32)
        .expect("b^k overflowed u64 — k too large for this base");
    let mut parts = vec![vec![Fp::ZERO; z.len()]; k];
    for (col, x) in z.iter().enumerate() {
        let v = x.balanced();
        let sign: i64 = if v < 0 { -1 } else { 1 };
        let mut mag = v.unsigned_abs();
        if mag >= bound {
            return Err(RingError::NormBound { got: mag, bound });
        }
        for part in parts.iter_mut() {
            let digit = (mag % b) as i64;
            part[col] = Fp::from_i64(sign * digit);
            mag /= b;
        }
    }
    Ok(parts)
}

/// Reconstruct `z = Σ_{i=1}^k b^{i−1} · z_i` from its base-`b` parts.
pub fn recompose_b(parts: &[Vec<Fp>], b: u64) -> Vec<Fp> {
    if parts.is_empty() {
        return Vec::new();
    }
    let len = parts[0].len();
    let mut out = vec![Fp::ZERO; len];
    // Horner from the most significant part: acc = acc·b + part.
    let b_fp = Fp::new(b);
    for part in parts.iter().rev() {
        for col in 0..len {
            out[col] = out[col] * b_fp + part[col];
        }
    }
    out
}

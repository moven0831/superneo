//! Coefficient maps and the SuperNeo embedding (Definitions 2, 7, 8).
//!
//! * `cf` / `cf_inv` / `ct`: coefficient vector, its inverse, and constant term.
//! * `embed_vector`: the SuperNeo embedding of a length-`d·n_R` field vector into a
//!   length-`n_R` ring vector (Definition 7) — split into `d`-blocks, each block the
//!   coefficient vector of one ring element.
//! * `bar_matrix_row`: the bar-lift `M̄` of a single matrix row (Definition 8).
//! * `mbar_row_dot_z`: the ring inner product `Σ_i M̄_row[i]·z[i]` whose constant
//!   term recovers `(Mz)_row` (Theorem 6).

use superneo_field::fp::Fp;

use crate::bar::bar_ring;
use crate::error::RingError;
use crate::ring::{RingElem, D};

/// `cf(a)`: the coefficient vector of a ring element (Definition 2).
#[inline]
pub fn cf(a: &RingElem) -> &[Fp; D] {
    a.coeffs()
}

/// `cf⁻¹`: the ring element whose coefficient vector is `c`.
#[inline]
pub fn cf_inv(c: [Fp; D]) -> RingElem {
    RingElem::from_coeffs(c)
}

/// `ct(a) = c_0`: the constant term (Definition 2).
#[inline]
pub fn ct(a: &RingElem) -> Fp {
    a.ct()
}

/// The SuperNeo embedding of a field vector into a ring vector (Definition 7).
///
/// `z ∈ F^{d·n_R}` is split into `n_R` consecutive `d`-blocks; block `i` becomes the
/// coefficient vector of ring element `z_i`. Errors if `len` is not a multiple of `d`.
pub fn embed_vector(z: &[Fp]) -> Result<Vec<RingElem>, RingError> {
    if !z.len().is_multiple_of(D) {
        return Err(RingError::Misaligned(z.len()));
    }
    Ok(z.chunks_exact(D)
        .map(|blk| {
            let mut c = [Fp::ZERO; D];
            c.copy_from_slice(blk);
            RingElem::from_coeffs(c)
        })
        .collect())
}

/// The bar-lift `M̄` of one matrix row (Definition 8): split the length-`d·n_R` row
/// into `d`-blocks and apply the bar transform to each, yielding `n_R` ring elements.
pub fn bar_matrix_row(row: &[Fp]) -> Result<Vec<RingElem>, RingError> {
    if !row.len().is_multiple_of(D) {
        return Err(RingError::Misaligned(row.len()));
    }
    Ok(row
        .chunks_exact(D)
        .map(|blk| {
            let mut c = [Fp::ZERO; D];
            c.copy_from_slice(blk);
            bar_ring(&c)
        })
        .collect())
}

/// The ring inner product `Σ_i bar_row[i] · z[i]`. Its constant term equals
/// `(Mz)_row` (Theorem 6) when `bar_row = M̄_row` and `z` is the embedding of the
/// witness. Panics if the two vectors differ in length.
pub fn mbar_row_dot_z(bar_row: &[RingElem], z: &[RingElem]) -> RingElem {
    assert_eq!(bar_row.len(), z.len(), "row/witness length mismatch");
    let mut acc = RingElem::ZERO;
    for (m, zi) in bar_row.iter().zip(z.iter()) {
        acc = acc + (*m * *zi);
    }
    acc
}

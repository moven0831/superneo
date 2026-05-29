//! Multilinear-extension utilities over the extension field `K = Ext2`.
//!
//! Convention: a length-`2^ℓ` value table is indexed so that variable `0` is the
//! **most significant** bit, i.e. `table[Σ_i x_i · 2^{ℓ−1−i}]`. The sum-check folds
//! variable `0` first. `mle_eval` and `eq_table` share this convention, so the
//! identity `Σ_x eq_table(r)[x]·vals[x] = mle_eval(vals, r)` holds.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_ring::{RingElem, RingK, D};

/// The smallest power of two `≥ n` (and `≥ 1`).
pub fn next_pow2(n: usize) -> usize {
    n.max(1).next_power_of_two()
}

/// `log2` of a power of two.
pub fn log2_pow2(m: usize) -> usize {
    debug_assert!(m.is_power_of_two());
    m.trailing_zeros() as usize
}

/// `eq(a, b) = Π_i (a_i·b_i + (1−a_i)(1−b_i))`.
pub fn eq(a: &[Ext2], b: &[Ext2]) -> Ext2 {
    assert_eq!(a.len(), b.len());
    let mut acc = Ext2::ONE;
    for (&ai, &bi) in a.iter().zip(b.iter()) {
        acc *= ai * bi + (Ext2::ONE - ai) * (Ext2::ONE - bi);
    }
    acc
}

/// The table `[eq(x, point)]_{x ∈ {0,1}^ℓ}` (length `2^ℓ`), variable 0 = MSB.
pub fn eq_table(point: &[Ext2]) -> Vec<Ext2> {
    let ell = point.len();
    let n = 1usize << ell;
    let mut tab = vec![Ext2::ZERO; n];
    for (idx, slot) in tab.iter_mut().enumerate() {
        let mut acc = Ext2::ONE;
        for (i, &pi) in point.iter().enumerate() {
            let bit = (idx >> (ell - 1 - i)) & 1;
            acc *= if bit == 1 { pi } else { Ext2::ONE - pi };
        }
        *slot = acc;
    }
    tab
}

/// Evaluate the multilinear extension of an extension-valued table at `point`.
///
/// `vals` is zero-padded to `2^point.len()` if shorter. Folds variable 0 (MSB) first.
pub fn mle_eval_ext(vals: &[Ext2], point: &[Ext2]) -> Ext2 {
    let ell = point.len();
    let mut size = 1usize << ell;
    let mut tab = vec![Ext2::ZERO; size];
    tab[..vals.len().min(size)].copy_from_slice(&vals[..vals.len().min(size)]);
    for &p in point.iter() {
        let half = size / 2;
        for idx in 0..half {
            tab[idx] = tab[idx] * (Ext2::ONE - p) + tab[half + idx] * p;
        }
        size = half;
    }
    tab[0]
}

/// Evaluate the multilinear extension of a field-valued table at `point`.
pub fn mle_eval_fp(vals: &[Fp], point: &[Ext2]) -> Ext2 {
    let ext: Vec<Ext2> = vals.iter().map(|&x| Ext2::from_base(x)).collect();
    mle_eval_ext(&ext, point)
}

/// Evaluate the multilinear extension of a **ring** vector at an extension point,
/// yielding an `R_K` element: `~(v)(point)` where `v ∈ R_F^len` (coefficient-wise).
///
/// Coefficient `ℓ` of the result is the MLE of `[cf(v_row)_ℓ]_row` at `point`.
pub fn mle_eval_ring(vals: &[RingElem], point: &[Ext2]) -> RingK {
    let mut coeffs = [Ext2::ZERO; D];
    let mut column = vec![Fp::ZERO; vals.len()];
    for (l, slot) in coeffs.iter_mut().enumerate() {
        for (row, v) in vals.iter().enumerate() {
            column[row] = v.coeffs()[l];
        }
        *slot = mle_eval_fp(&column, point);
    }
    RingK::from_coeffs(coeffs)
}

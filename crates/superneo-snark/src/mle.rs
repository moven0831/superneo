//! Little-endian multilinear helpers shared by the PCS and the Spartan reduction.
//!
//! These mirror `superneo_fold::multilinear` but bind variable `x_k` to bit `k` (LSB
//! first) so a sum-check's reduced point is directly usable as a BaseFold evaluation
//! point. (The fold crate's versions are MSB-first; the two conventions are not
//! interchangeable, which is why this small parallel set exists.)

use superneo_field::ext2::Ext2;

/// `eq(point, x)` evaluation table over `{0,1}^ν`, little-endian (`x_k` = bit `k`).
pub fn eq_table(point: &[Ext2]) -> Vec<Ext2> {
    let nu = point.len();
    let mut tab = vec![Ext2::ONE; 1 << nu];
    for (k, &pk) in point.iter().enumerate() {
        let one_minus = Ext2::ONE - pk;
        for (i, slot) in tab.iter_mut().enumerate() {
            *slot *= if (i >> k) & 1 == 1 { pk } else { one_minus };
        }
    }
    tab
}

/// `eq(a, b) = Π_k (a_k·b_k + (1−a_k)(1−b_k))`.
pub fn eq_eval(a: &[Ext2], b: &[Ext2]) -> Ext2 {
    a.iter()
        .zip(b.iter())
        .map(|(&ak, &bk)| ak * bk + (Ext2::ONE - ak) * (Ext2::ONE - bk))
        .fold(Ext2::ONE, |x, y| x * y)
}

/// Bind the LSB of an evaluation table to `alpha`: `t'[j] = (1−α)·t[2j] + α·t[2j+1]`.
pub fn fold_evals(table: &[Ext2], alpha: Ext2) -> Vec<Ext2> {
    table
        .chunks_exact(2)
        .map(|p| p[0] + alpha * (p[1] - p[0]))
        .collect()
}

/// Evaluate the multilinear extension of `evals` (little-endian) at `point`, LSB first.
pub fn eval(evals: &[Ext2], point: &[Ext2]) -> Ext2 {
    let mut cur = evals.to_vec();
    for &pk in point {
        cur = fold_evals(&cur, pk);
    }
    cur[0]
}

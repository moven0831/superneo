//! The SuperNeo inner-product (bar) transform (Theorem 5).
//!
//! Let `G ∈ F^{d×d}` be the symmetric Gram matrix `G[i][j] = ct(X^{i+j} mod Φ)`.
//! The transform is `bar(a) := G⁻¹·a`; then for ring elements with coefficient
//! vectors `bar(a)` and `b`,
//!
//! ```text
//!   ct( cf⁻¹(bar(a)) · cf⁻¹(b) ) = ⟨a, b⟩.
//! ```
//!
//! The paper proves `det(G) = ±1`, so `G` is invertible over `F`. We build `G⁻¹`
//! once (lazily) and assert `G⁻¹·G = I` as an always-on post-condition — a failure
//! there signals a ring-reduction bug, not a linear-algebra one (plan, Risk R2).

use std::sync::OnceLock;

use superneo_field::fp::Fp;

use crate::ring::{RingElem, D};

/// The Gram matrix `G[i][j] = ct(X^{i+j} mod Φ)`.
pub fn gram() -> [[Fp; D]; D] {
    // ct(X^k mod Φ) for k = 0..=2D−2.
    let mut ct_pow = [Fp::ZERO; 2 * D - 1];
    for (k, slot) in ct_pow.iter_mut().enumerate() {
        *slot = RingElem::x_pow_mod_phi(k).ct();
    }
    // Row i of the Gram matrix is exactly ct_pow[i ..= i+D-1] (since G[i][j] = ct_pow[i+j]).
    let mut g = [[Fp::ZERO; D]; D];
    for i in 0..D {
        g[i].copy_from_slice(&ct_pow[i..i + D]);
    }
    g
}

/// Invert a `d×d` matrix over `F` by Gauss–Jordan elimination; `None` if singular.
fn invert(mat: &[[Fp; D]; D]) -> Option<[[Fp; D]; D]> {
    let mut a = *mat;
    let mut inv = identity();
    for col in 0..D {
        // Find a pivot at or below the diagonal.
        let piv = (col..D).find(|&r| !a[r][col].is_zero())?;
        a.swap(col, piv);
        inv.swap(col, piv);
        // Normalize the pivot row.
        let inv_p = a[col][col].inv().ok()?;
        for j in 0..D {
            a[col][j] *= inv_p;
            inv[col][j] *= inv_p;
        }
        // Eliminate the column from all other rows.
        for r in 0..D {
            if r == col {
                continue;
            }
            let f = a[r][col];
            if f.is_zero() {
                continue;
            }
            for j in 0..D {
                a[r][j] -= f * a[col][j];
                inv[r][j] -= f * inv[col][j];
            }
        }
    }
    Some(inv)
}

/// The determinant of a `d×d` matrix over `F` (Gaussian elimination, tracking sign).
pub fn determinant(mat: &[[Fp; D]; D]) -> Fp {
    let mut a = *mat;
    let mut det = Fp::ONE;
    for col in 0..D {
        let piv = match (col..D).find(|&r| !a[r][col].is_zero()) {
            Some(p) => p,
            None => return Fp::ZERO,
        };
        if piv != col {
            a.swap(col, piv);
            det = -det;
        }
        det *= a[col][col];
        let inv_p = a[col][col].inv().expect("nonzero pivot");
        for r in (col + 1)..D {
            let f = a[r][col] * inv_p;
            if f.is_zero() {
                continue;
            }
            for j in col..D {
                a[r][j] -= f * a[col][j];
            }
        }
    }
    det
}

fn identity() -> [[Fp; D]; D] {
    let mut m = [[Fp::ZERO; D]; D];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = Fp::ONE;
    }
    m
}

/// The cached bar-transform matrix `G⁻¹`.
pub fn bar_matrix() -> &'static [[Fp; D]; D] {
    static M: OnceLock<[[Fp; D]; D]> = OnceLock::new();
    M.get_or_init(|| {
        let g = gram();
        let inv = invert(&g).expect("Gram matrix must be invertible (Theorem 5: det = ±1)");
        // Always-on post-condition: G⁻¹·G = I.
        for i in 0..D {
            for j in 0..D {
                let mut acc = Fp::ZERO;
                for k in 0..D {
                    acc += inv[i][k] * g[k][j];
                }
                let want = if i == j { Fp::ONE } else { Fp::ZERO };
                assert_eq!(
                    acc, want,
                    "bar_matrix · Gram != I at ({i}, {j}) — ring bug?"
                );
            }
        }
        inv
    })
}

/// Apply the bar transform to one length-`d` block: `bar(v) = G⁻¹·v`.
pub fn bar_block(v: &[Fp; D]) -> [Fp; D] {
    let m = bar_matrix();
    let mut out = [Fp::ZERO; D];
    for i in 0..D {
        let mut acc = Fp::ZERO;
        for j in 0..D {
            acc += m[i][j] * v[j];
        }
        out[i] = acc;
    }
    out
}

/// Apply the bar transform to a length-`d` block, returning a [`RingElem`].
pub fn bar_ring(v: &[Fp; D]) -> RingElem {
    RingElem::from_coeffs(bar_block(v))
}

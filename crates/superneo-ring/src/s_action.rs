//! The rotation / S-action matrix `rot(a) ∈ F^{d×d}` satisfying
//! `cf(a·b) = rot(a)·cf(b)` (the matrix form of ring multiplication).
//!
//! Column `j` of `rot(a)` is `cf(a · X^j mod Φ)`, since
//! `a·b = Σ_j b_j · (a·X^j)` and `cf` is linear.

use superneo_field::fp::Fp;

use crate::ring::{RingElem, D};

/// Compute the `d×d` rotation matrix `rot(a)` with `rot(a)·cf(b) = cf(a·b)`.
pub fn rot(a: &RingElem) -> [[Fp; D]; D] {
    let mut m = [[Fp::ZERO; D]; D];
    for col in 0..D {
        let shifted = a.mul_x_pow(col); // a · X^col mod Φ
        let c = shifted.coeffs();
        for row in 0..D {
            m[row][col] = c[row];
        }
    }
    m
}

/// Apply a `d×d` matrix to a length-`d` coefficient vector (`M·v`).
pub fn matvec(m: &[[Fp; D]; D], v: &[Fp; D]) -> [Fp; D] {
    let mut out = [Fp::ZERO; D];
    for (row, mrow) in m.iter().enumerate() {
        let mut acc = Fp::ZERO;
        for col in 0..D {
            acc += mrow[col] * v[col];
        }
        out[row] = acc;
    }
    out
}

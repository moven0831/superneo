//! The extension ring `R_K = K[X]/(Φ)`, `Φ = X^54 + X^27 + 1`.
//!
//! Evaluation claims `y_j = ~(M̄_j z)(r) ∈ R_K` live here (Definition 13): a ring
//! element whose `d` coefficients are extension-field (`K = Ext2`) values. We need
//! addition, the constant term `ct ∈ K`, scaling by `F`/`K`, and the mixed product
//! `R_K × R_F → R_K` (used to apply the Π_RLC challenges `ρ ∈ R_F` and the Π_DEC
//! weights). `F_q ⊂ K` gives the embedding `R_F ↪ R_K`.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;

use crate::ring::{RingElem, D};

const MID: usize = D / 2; // 27
const PROD_LEN: usize = 2 * D - 1; // 107

/// An element of `R_K`, stored as its coefficient vector `[c_0, …, c_{53}] ∈ K^d`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct RingK([Ext2; D]);

impl RingK {
    /// The zero element.
    pub const ZERO: RingK = RingK([Ext2::ZERO; D]);

    /// Build from a full coefficient array.
    #[inline]
    pub const fn from_coeffs(c: [Ext2; D]) -> RingK {
        RingK(c)
    }

    /// Embed an `R_F` element into `R_K` (coefficient-wise `F ↪ K`).
    pub fn from_rf(a: &RingElem) -> RingK {
        let mut c = [Ext2::ZERO; D];
        for (ci, &ai) in c.iter_mut().zip(a.coeffs().iter()) {
            *ci = Ext2::from_base(ai);
        }
        RingK(c)
    }

    /// A constant ring element `k ∈ K` (degree-0).
    pub fn from_const(k: Ext2) -> RingK {
        let mut c = [Ext2::ZERO; D];
        c[0] = k;
        RingK(c)
    }

    /// Borrow the coefficient array.
    #[inline]
    pub fn coeffs(&self) -> &[Ext2; D] {
        &self.0
    }

    /// The `ℓ`-th coefficient `cf(·)_ℓ`.
    #[inline]
    pub fn coeff(&self, l: usize) -> Ext2 {
        self.0[l]
    }

    /// The constant term `ct(·) = c_0 ∈ K` (Definition 2).
    #[inline]
    pub fn ct(&self) -> Ext2 {
        self.0[0]
    }

    /// Component-wise sum.
    pub fn add(&self, other: &RingK) -> RingK {
        let mut c = [Ext2::ZERO; D];
        for i in 0..D {
            c[i] = self.0[i] + other.0[i];
        }
        RingK(c)
    }

    /// Scale by a base-field element `s ∈ F`.
    pub fn scale_fp(&self, s: Fp) -> RingK {
        let se = Ext2::from_base(s);
        let mut c = [Ext2::ZERO; D];
        for i in 0..D {
            c[i] = self.0[i] * se;
        }
        RingK(c)
    }

    /// Scale by an extension element `k ∈ K`.
    pub fn scale_ext(&self, k: Ext2) -> RingK {
        let mut c = [Ext2::ZERO; D];
        for i in 0..D {
            c[i] = self.0[i] * k;
        }
        RingK(c)
    }

    /// The mixed ring product `self · rf` for `rf ∈ R_F`, reduced modulo `Φ`.
    ///
    /// Used to apply the Π_RLC challenges `ρ ∈ R_F` to the evaluation claims
    /// `y_j ∈ R_K`. Same trinomial reduction as `R_F`, over `K`-coefficients.
    pub fn mul_rf(&self, rf: &RingElem) -> RingK {
        let mut prod = [Ext2::ZERO; PROD_LEN];
        let rfc = rf.coeffs();
        for i in 0..D {
            let ai = self.0[i];
            if ai == Ext2::ZERO {
                continue;
            }
            for j in 0..D {
                if rfc[j].is_zero() {
                    continue;
                }
                prod[i + j] += ai * Ext2::from_base(rfc[j]);
            }
        }
        // Reduce: X^i = −X^{i−27} − X^{i−54}, reverse pass (see ring::reduce).
        for i in (D..PROD_LEN).rev() {
            let t = prod[i];
            if t != Ext2::ZERO {
                prod[i] = Ext2::ZERO;
                prod[i - MID] -= t;
                prod[i - D] -= t;
            }
        }
        let mut out = [Ext2::ZERO; D];
        out.copy_from_slice(&prod[0..D]);
        RingK(out)
    }
}

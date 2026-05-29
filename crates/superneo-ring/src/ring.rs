//! The cyclotomic ring `R_F = F[X]/(Φ)` with `Φ = X^54 + X^27 + 1` (the 81st
//! cyclotomic polynomial, degree `d = 54`).
//!
//! Reduction uses the identity `X^54 = −X^27 − 1`, applied to product degrees
//! `[54, 106]` in **reverse** (high to low) so that the cascade `X^i → −X^{i−27} − X^{i−54}`
//! — where `i−27` may itself still be `≥ 54` for `i ≥ 81` — is resolved by the time
//! the loop reaches that lower degree (see plan, Risk R1).

use core::ops::{Add, Mul, Neg, Sub};

use superneo_field::fp::Fp;

/// Degree of the cyclotomic ring (`φ(81) = 54`).
pub const D: usize = 54;

/// Exponent of the middle term of `Φ = X^54 + X^27 + 1`.
const MID: usize = D / 2; // 27

/// Length of an unreduced product of two degree-`<D` polynomials (`2D − 1 = 107`).
const PROD_LEN: usize = 2 * D - 1;

/// An element of `R_F`, stored as its coefficient vector `[c_0, …, c_{53}]`.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct RingElem([Fp; D]);

impl RingElem {
    /// The zero ring element.
    pub const ZERO: RingElem = RingElem([Fp::ZERO; D]);

    /// The multiplicative identity `1`.
    pub fn one() -> RingElem {
        let mut c = [Fp::ZERO; D];
        c[0] = Fp::ONE;
        RingElem(c)
    }

    /// The generator `X`.
    pub fn x() -> RingElem {
        let mut c = [Fp::ZERO; D];
        c[1] = Fp::ONE;
        RingElem(c)
    }

    /// Build from a full coefficient array.
    #[inline]
    pub const fn from_coeffs(c: [Fp; D]) -> RingElem {
        RingElem(c)
    }

    /// Borrow the coefficient array.
    #[inline]
    pub fn coeffs(&self) -> &[Fp; D] {
        &self.0
    }

    /// Consume into the coefficient array.
    #[inline]
    pub fn into_coeffs(self) -> [Fp; D] {
        self.0
    }

    /// The constant term `ct(·) = c_0` (Definition 2).
    #[inline]
    pub fn ct(&self) -> Fp {
        self.0[0]
    }

    /// Whether this is the zero element.
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|c| c.is_zero())
    }

    /// Reduce an unreduced coefficient buffer (`len = 2D − 1`) modulo `Φ`.
    ///
    /// Mutates `prod` in place and returns the reduced element. Uses the reverse
    /// pass described in the module docs.
    fn reduce(prod: &mut [Fp; PROD_LEN]) -> RingElem {
        for i in (D..PROD_LEN).rev() {
            let t = prod[i];
            if !t.is_zero() {
                prod[i] = Fp::ZERO;
                // X^i = X^{i−54}·X^54 = X^{i−54}·(−X^27 − 1) = −X^{i−27} − X^{i−54}.
                prod[i - MID] -= t;
                prod[i - D] -= t;
            }
        }
        let mut out = [Fp::ZERO; D];
        out.copy_from_slice(&prod[0..D]);
        RingElem(out)
    }

    /// `X^k mod Φ` for `0 ≤ k ≤ 2D − 2`.
    pub fn x_pow_mod_phi(k: usize) -> RingElem {
        assert!(k < PROD_LEN, "x_pow_mod_phi: k={k} out of range");
        let mut prod = [Fp::ZERO; PROD_LEN];
        prod[k] = Fp::ONE;
        RingElem::reduce(&mut prod)
    }

    /// `self · X^k mod Φ` for `0 ≤ k < D` (the column generator of the rotation map).
    pub fn mul_x_pow(&self, k: usize) -> RingElem {
        assert!(k < D, "mul_x_pow: k={k} out of range");
        let mut prod = [Fp::ZERO; PROD_LEN];
        prod[k..k + D].copy_from_slice(&self.0);
        RingElem::reduce(&mut prod)
    }

    /// Scale by a base-field scalar.
    pub fn scale(&self, s: Fp) -> RingElem {
        let mut out = [Fp::ZERO; D];
        for i in 0..D {
            out[i] = self.0[i] * s;
        }
        RingElem(out)
    }
}

impl Add for RingElem {
    type Output = RingElem;
    fn add(self, rhs: RingElem) -> RingElem {
        let mut out = [Fp::ZERO; D];
        for i in 0..D {
            out[i] = self.0[i] + rhs.0[i];
        }
        RingElem(out)
    }
}

impl Sub for RingElem {
    type Output = RingElem;
    fn sub(self, rhs: RingElem) -> RingElem {
        let mut out = [Fp::ZERO; D];
        for i in 0..D {
            out[i] = self.0[i] - rhs.0[i];
        }
        RingElem(out)
    }
}

impl Neg for RingElem {
    type Output = RingElem;
    fn neg(self) -> RingElem {
        let mut out = [Fp::ZERO; D];
        for i in 0..D {
            out[i] = -self.0[i];
        }
        RingElem(out)
    }
}

impl Mul for RingElem {
    type Output = RingElem;
    fn mul(self, rhs: RingElem) -> RingElem {
        let mut prod = [Fp::ZERO; PROD_LEN];
        for i in 0..D {
            let ai = self.0[i];
            if ai.is_zero() {
                continue;
            }
            for j in 0..D {
                prod[i + j] += ai * rhs.0[j];
            }
        }
        RingElem::reduce(&mut prod)
    }
}

impl core::fmt::Debug for RingElem {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "RingElem(")?;
        let mut first = true;
        for (i, c) in self.0.iter().enumerate() {
            if !c.is_zero() {
                if !first {
                    write!(f, " + ")?;
                }
                write!(f, "{}·X^{}", c, i)?;
                first = false;
            }
        }
        if first {
            write!(f, "0")?;
        }
        write!(f, ")")
    }
}

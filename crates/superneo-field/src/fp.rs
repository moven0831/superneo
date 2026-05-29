//! The Goldilocks prime field `F_q`, `q = 2^64 − 2^32 + 1`.
//!
//! Elements are stored as **canonical** residues in `[0, q)`. Multiplication uses
//! Solinas reduction built on the identities
//!
//! ```text
//!   2^64 ≡ 2^32 − 1   (mod q)        // = EPSILON
//!   2^96 ≡ −1         (mod q)
//! ```
//!
//! See Definition 1 and Appendix B.2 of the paper. The reduction is cross-checked
//! against the reference `x mod q` (computed in `u128`) in the test suite.

use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::error::FieldError;
use crate::Q;

/// `2^32 − 1`, the residue of `2^64` modulo `q`.
const EPSILON: u64 = (1 << 32) - 1;

/// An element of the Goldilocks field `F_q`, stored canonically in `[0, q)`.
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct Fp(u64);

impl Fp {
    /// The additive identity.
    pub const ZERO: Fp = Fp(0);
    /// The multiplicative identity.
    pub const ONE: Fp = Fp(1);
    /// The modulus, as an `Fp`-adjacent constant for callers.
    pub const MODULUS: u64 = Q;

    /// Construct from a `u64`, reducing into `[0, q)`.
    ///
    /// One conditional subtraction suffices because any `u64` is `< 2q`.
    #[inline]
    pub const fn new(x: u64) -> Self {
        Fp(if x >= Q { x - Q } else { x })
    }

    /// Construct from a value already known to be canonical (`< q`).
    ///
    /// Returns [`FieldError::NonCanonical`] otherwise.
    #[inline]
    pub const fn from_canonical(x: u64) -> Result<Self, FieldError> {
        if x < Q {
            Ok(Fp(x))
        } else {
            Err(FieldError::NonCanonical(x))
        }
    }

    /// The canonical `u64` representative in `[0, q)`.
    #[inline]
    pub const fn to_u64(self) -> u64 {
        self.0
    }

    /// Whether this is the zero element.
    #[inline]
    pub const fn is_zero(self) -> bool {
        self.0 == 0
    }

    /// The balanced representative in `(−q/2, q/2]`.
    ///
    /// Magnitude is `< q/2 < 2^63`, so it fits in `i64` (Definition 3).
    #[inline]
    pub fn balanced(self) -> i64 {
        if self.0 > Q / 2 {
            // value − q, computed in i128 to avoid overflow, then narrowed.
            (self.0 as i128 - Q as i128) as i64
        } else {
            self.0 as i64
        }
    }

    /// The ℓ∞ magnitude `|balanced|` (Definition 3).
    #[inline]
    pub fn abs_balanced(self) -> u64 {
        self.balanced().unsigned_abs()
    }

    /// Reduce a 128-bit value modulo `q` (Solinas).
    ///
    /// Splits `x = x_lo + 2^64·x_hi`, and `x_hi = 2^32·hi_hi + hi_lo`, then applies
    /// `2^64 ≡ EPSILON` and `2^96 ≡ −1`.
    #[inline]
    pub const fn reduce128(x: u128) -> Self {
        let lo = x as u64;
        let hi = (x >> 64) as u64;
        let hi_hi = hi >> 32;
        let hi_lo = hi & EPSILON;

        // 2^96 ≡ −1: subtract hi_hi, correcting a borrow by −EPSILON.
        let (mut t0, borrow) = lo.overflowing_sub(hi_hi);
        if borrow {
            t0 = t0.wrapping_sub(EPSILON);
        }
        // 2^64 ≡ EPSILON: add hi_lo·EPSILON (< 2^64), correcting a carry by +EPSILON.
        let t1 = hi_lo * EPSILON;
        let (sum, carry) = t0.overflowing_add(t1);
        let res = if carry { sum.wrapping_add(EPSILON) } else { sum };
        Fp::new(res)
    }

    /// `self²`.
    #[inline]
    pub fn square(self) -> Self {
        self * self
    }

    /// `self^exp` via square-and-multiply.
    pub fn pow(self, mut exp: u64) -> Self {
        let mut base = self;
        let mut acc = Fp::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                acc *= base;
            }
            base = base.square();
            exp >>= 1;
        }
        acc
    }

    /// Multiplicative inverse via Fermat's little theorem (`a^{q−2}`).
    ///
    /// Returns [`FieldError::InverseOfZero`] for zero.
    pub fn inv(self) -> Result<Self, FieldError> {
        if self.is_zero() {
            return Err(FieldError::InverseOfZero);
        }
        Ok(self.pow(Q - 2))
    }
}

// ---- Arithmetic operator impls ------------------------------------------------

impl Add for Fp {
    type Output = Fp;
    #[inline]
    fn add(self, rhs: Fp) -> Fp {
        // a + b < 2q < 2^65, so compute in u128 and subtract q at most once.
        let s = self.0 as u128 + rhs.0 as u128;
        let q = Q as u128;
        Fp(if s >= q { (s - q) as u64 } else { s as u64 })
    }
}

impl Sub for Fp {
    type Output = Fp;
    #[inline]
    fn sub(self, rhs: Fp) -> Fp {
        let (r, borrow) = self.0.overflowing_sub(rhs.0);
        // On borrow, r = a − b + 2^64; adding q (mod 2^64) yields a − b + q ∈ (0, q).
        Fp(if borrow { r.wrapping_add(Q) } else { r })
    }
}

impl Neg for Fp {
    type Output = Fp;
    #[inline]
    fn neg(self) -> Fp {
        if self.0 == 0 {
            Fp::ZERO
        } else {
            Fp(Q - self.0)
        }
    }
}

impl Mul for Fp {
    type Output = Fp;
    #[inline]
    fn mul(self, rhs: Fp) -> Fp {
        Fp::reduce128(self.0 as u128 * rhs.0 as u128)
    }
}

impl AddAssign for Fp {
    #[inline]
    fn add_assign(&mut self, rhs: Fp) {
        *self = *self + rhs;
    }
}
impl SubAssign for Fp {
    #[inline]
    fn sub_assign(&mut self, rhs: Fp) {
        *self = *self - rhs;
    }
}
impl MulAssign for Fp {
    #[inline]
    fn mul_assign(&mut self, rhs: Fp) {
        *self = *self * rhs;
    }
}

impl From<u64> for Fp {
    #[inline]
    fn from(x: u64) -> Fp {
        Fp::new(x)
    }
}

impl fmt::Debug for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Fp({})", self.0)
    }
}

impl fmt::Display for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(feature = "serde")]
impl serde::Serialize for Fp {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_u64(self.0)
    }
}

#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Fp {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let x = u64::deserialize(d)?;
        Fp::from_canonical(x).map_err(serde::de::Error::custom)
    }
}

//! The degree-2 extension `K = F_{q^2} = F_q[u]/(u² − 7)`.
//!
//! `7` is a quadratic non-residue modulo the Goldilocks prime (verified in tests),
//! so `u² − 7` is irreducible and the quotient is a field. The sum-check protocol
//! runs its challenges over `K`, where the soundness error `ℓ·d/|K|` is negligible
//! (Definition 6); `F_q ⊂ K` as the constant sub-field (Definition 1).

use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::error::FieldError;
use crate::fp::Fp;

/// The non-residue `W = 7` defining the extension `u² = W`.
const W: Fp = Fp::new(7);

/// An element `c0 + c1·u` of `K = F_{q^2}`.
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct Ext2 {
    /// Constant coefficient.
    pub c0: Fp,
    /// Coefficient of `u`.
    pub c1: Fp,
}

impl Ext2 {
    /// The additive identity.
    pub const ZERO: Ext2 = Ext2 {
        c0: Fp::ZERO,
        c1: Fp::ZERO,
    };
    /// The multiplicative identity.
    pub const ONE: Ext2 = Ext2 {
        c0: Fp::ONE,
        c1: Fp::ZERO,
    };

    /// Construct `c0 + c1·u`.
    #[inline]
    pub const fn new(c0: Fp, c1: Fp) -> Self {
        Ext2 { c0, c1 }
    }

    /// Embed a base-field element as a constant polynomial (`F_q ↪ K`).
    #[inline]
    pub const fn from_base(a: Fp) -> Self {
        Ext2 {
            c0: a,
            c1: Fp::ZERO,
        }
    }

    /// Whether this is the zero element.
    #[inline]
    pub const fn is_zero(self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// The Frobenius conjugate `c0 − c1·u` (the nontrivial automorphism `x ↦ x^q`).
    #[inline]
    pub fn conjugate(self) -> Self {
        Ext2 {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    /// The field norm `N(z) = z·conj(z) = c0² − W·c1² ∈ F_q`.
    #[inline]
    pub fn norm(self) -> Fp {
        self.c0.square() - W * self.c1.square()
    }

    /// Multiply by a base-field scalar.
    #[inline]
    pub fn mul_base(self, s: Fp) -> Self {
        Ext2 {
            c0: self.c0 * s,
            c1: self.c1 * s,
        }
    }

    /// `self²`.
    #[inline]
    pub fn square(self) -> Self {
        self * self
    }

    /// `self^exp` via square-and-multiply.
    pub fn pow(self, mut exp: u64) -> Self {
        let mut base = self;
        let mut acc = Ext2::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                acc *= base;
            }
            base = base.square();
            exp >>= 1;
        }
        acc
    }

    /// Multiplicative inverse: `conj(z) · N(z)^{-1}`.
    ///
    /// Returns [`FieldError::InverseOfZero`] for zero.
    pub fn inv(self) -> Result<Self, FieldError> {
        if self.is_zero() {
            return Err(FieldError::InverseOfZero);
        }
        let n_inv = self.norm().inv()?;
        Ok(self.conjugate().mul_base(n_inv))
    }
}

// ---- Arithmetic operator impls ------------------------------------------------

impl Add for Ext2 {
    type Output = Ext2;
    #[inline]
    fn add(self, rhs: Ext2) -> Ext2 {
        Ext2 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl Sub for Ext2 {
    type Output = Ext2;
    #[inline]
    fn sub(self, rhs: Ext2) -> Ext2 {
        Ext2 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl Neg for Ext2 {
    type Output = Ext2;
    #[inline]
    fn neg(self) -> Ext2 {
        Ext2 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

impl Mul for Ext2 {
    type Output = Ext2;
    #[inline]
    fn mul(self, rhs: Ext2) -> Ext2 {
        // (a + b·u)(c + d·u) = (ac + W·bd) + (ad + bc)·u, since u² = W.
        let ac = self.c0 * rhs.c0;
        let bd = self.c1 * rhs.c1;
        let ad = self.c0 * rhs.c1;
        let bc = self.c1 * rhs.c0;
        Ext2 {
            c0: ac + W * bd,
            c1: ad + bc,
        }
    }
}

impl AddAssign for Ext2 {
    #[inline]
    fn add_assign(&mut self, rhs: Ext2) {
        *self = *self + rhs;
    }
}
impl SubAssign for Ext2 {
    #[inline]
    fn sub_assign(&mut self, rhs: Ext2) {
        *self = *self - rhs;
    }
}
impl MulAssign for Ext2 {
    #[inline]
    fn mul_assign(&mut self, rhs: Ext2) {
        *self = *self * rhs;
    }
}

impl From<Fp> for Ext2 {
    #[inline]
    fn from(a: Fp) -> Ext2 {
        Ext2::from_base(a)
    }
}

impl fmt::Debug for Ext2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Ext2({} + {}·u)", self.c0, self.c1)
    }
}

#[cfg(feature = "serde")]
impl serde::Serialize for Ext2 {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeTuple;
        let mut t = s.serialize_tuple(2)?;
        t.serialize_element(&self.c0)?;
        t.serialize_element(&self.c1)?;
        t.end()
    }
}

#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Ext2 {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let (c0, c1) = <(Fp, Fp)>::deserialize(d)?;
        Ok(Ext2 { c0, c1 })
    }
}

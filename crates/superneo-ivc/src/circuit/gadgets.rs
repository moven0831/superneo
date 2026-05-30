//! Extension-field (`K = F_q[u]/(u²−7)`) gadgets over the R1CS builder (M7).
//!
//! A `K` value is a pair of base-field linear combinations `(c0, c1)` for `c0 + c1·u`.
//! Addition and scaling are free (linear); multiplication uses Karatsuba — three base
//! multiplications, matching the recursion-cost accounting in the plan:
//! `(a0+a1u)(b0+b1u) = (a0b0 + 7·a1b1) + ((a0+a1)(b0+b1) − a0b0 − a1b1)·u`.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;

use super::cs::{ConstraintSystem, Lc};

/// The non-residue `u² = 7`.
const W: u64 = 7;

/// A `K`-valued wire `c0 + c1·u`, each part a base-field linear combination.
#[derive(Clone, Debug)]
pub struct KVar {
    /// Constant part.
    pub c0: Lc,
    /// Coefficient of `u`.
    pub c1: Lc,
}

impl KVar {
    /// A constant `K` value.
    pub fn constant(cs: &ConstraintSystem, v: Ext2) -> KVar {
        KVar {
            c0: cs.constant(v.c0),
            c1: cs.constant(v.c1),
        }
    }

    /// Allocate a `K` value as two witness variables (advice).
    pub fn alloc(cs: &mut ConstraintSystem, v: Ext2) -> KVar {
        KVar {
            c0: Lc::from_var(cs.alloc(v.c0)),
            c1: Lc::from_var(cs.alloc(v.c1)),
        }
    }

    /// `self + other` (free).
    pub fn add(&self, other: &KVar) -> KVar {
        KVar {
            c0: self.c0.add(&other.c0),
            c1: self.c1.add(&other.c1),
        }
    }

    /// `self − other` (free).
    pub fn sub(&self, other: &KVar) -> KVar {
        KVar {
            c0: self.c0.sub(&other.c0),
            c1: self.c1.sub(&other.c1),
        }
    }

    /// `s · self` for a base-field constant `s` (free).
    pub fn scale(&self, s: Fp) -> KVar {
        KVar {
            c0: self.c0.scale(s),
            c1: self.c1.scale(s),
        }
    }

    /// `self · other` via Karatsuba (allocates three product variables).
    pub fn mul(&self, cs: &mut ConstraintSystem, other: &KVar) -> KVar {
        let m0 = Lc::from_var(cs.mul(&self.c0, &other.c0)); // a0·b0
        let m1 = Lc::from_var(cs.mul(&self.c1, &other.c1)); // a1·b1
        let m2 = Lc::from_var(cs.mul(&self.c0.add(&self.c1), &other.c0.add(&other.c1)));
        KVar {
            c0: m0.add(&m1.scale(Fp::new(W))), // a0b0 + 7·a1b1
            c1: m2.sub(&m0).sub(&m1),          // (a0+a1)(b0+b1) − a0b0 − a1b1
        }
    }

    /// Impose `self = other` (two base-field equality constraints).
    pub fn assert_eq(&self, cs: &mut ConstraintSystem, other: &KVar) {
        cs.assert_eq(&self.c0, &other.c0);
        cs.assert_eq(&self.c1, &other.c1);
    }
}

//! In-circuit ring gadgets (M7, Phase P2 support): `R_F` and `R_K` arithmetic.
//!
//! Mirrors `superneo_ring`. An [`RFVar`] is `D = 54` base-field coefficient wires for an
//! element of `R_F = F[X]/(X⁵⁴+X²⁷+1)`; an [`RKVar`] is `D` extension coefficient
//! ([`KVar`]) wires for an element of `R_K`. Multiplication is schoolbook into a `2D−1`
//! buffer followed by the trinomial reduction `X⁵⁴ = −X²⁷ − 1` (the reverse pass in
//! `ring::reduce` / `ring_k::mul_rf`).
//!
//! Cost note: the reduction is linear (only adds/subs — free), but schoolbook
//! multiplication costs `D²` base multiplications per `R_F` product (and `2D²` for
//! `R_K · R_F`). That is why the Π_RLC commitment/eval combine — `κ·(K+k)` ring products —
//! must run at scaled-down `κ, K, n_R` in this PoC; it is the dominant recursion cost.

use superneo_field::ext2::Ext2;
use superneo_ring::{RingElem, RingK, D};

use super::cs::{ConstraintSystem, Lc};
use super::gadgets::KVar;

const MID: usize = D / 2; // 27, the middle exponent of Φ = X^D + X^MID + 1
const PROD_LEN: usize = 2 * D - 1; // 107, the unreduced product length

/// Trinomial reduction of an unreduced length-`PROD_LEN` buffer, in place: for `i` from
/// the top down, `X^i = −X^{i−MID} − X^{i−D}`. Generic over the (linear) coefficient ops
/// so it serves both `R_F` (over `Lc`) and `R_K` (over `KVar`). Returns the low `D` slots.
fn reduce<T: Clone>(prod: &mut [T], sub: impl Fn(&T, &T) -> T) -> Vec<T> {
    for i in (D..PROD_LEN).rev() {
        let t = prod[i].clone();
        prod[i - MID] = sub(&prod[i - MID], &t);
        prod[i - D] = sub(&prod[i - D], &t);
    }
    prod[0..D].to_vec()
}

/// An `R_F` element: `D` base-field coefficient wires.
#[derive(Clone)]
pub struct RFVar {
    /// Coefficients `c_0, …, c_{D-1}` (each a base-field linear combination).
    pub coeffs: Vec<Lc>,
}

impl RFVar {
    /// A constant `R_F` value (all coefficients pinned).
    pub fn constant(cs: &ConstraintSystem, v: &RingElem) -> RFVar {
        RFVar {
            coeffs: v.coeffs().iter().map(|&c| cs.constant(c)).collect(),
        }
    }

    /// Allocate an `R_F` value as `D` advice wires.
    pub fn alloc(cs: &mut ConstraintSystem, v: &RingElem) -> RFVar {
        RFVar {
            coeffs: v.coeffs().iter().map(|&c| Lc::from_var(cs.alloc(c))).collect(),
        }
    }

    /// `self + other` (free).
    pub fn add(&self, other: &RFVar) -> RFVar {
        RFVar {
            coeffs: self.coeffs.iter().zip(&other.coeffs).map(|(a, b)| a.add(b)).collect(),
        }
    }

    /// `self · other mod Φ` (schoolbook `D²` base mults, then linear reduction).
    pub fn mul(&self, cs: &mut ConstraintSystem, other: &RFVar) -> RFVar {
        let mut prod: Vec<Lc> = vec![Lc::zero(); PROD_LEN];
        for i in 0..D {
            for j in 0..D {
                let p = cs.mul(&self.coeffs[i], &other.coeffs[j]);
                prod[i + j] = prod[i + j].add(&Lc::from_var(p));
            }
        }
        RFVar {
            coeffs: reduce(&mut prod, |a, b| a.sub(b)),
        }
    }

    /// Impose `self = other` coefficient-wise.
    pub fn assert_eq(&self, cs: &mut ConstraintSystem, other: &RFVar) {
        for (a, b) in self.coeffs.iter().zip(&other.coeffs) {
            cs.assert_eq(a, b);
        }
    }
}

/// An `R_K` element: `D` extension-field coefficient wires.
#[derive(Clone)]
pub struct RKVar {
    /// Coefficients `c_0, …, c_{D-1}` (each a `K`-valued [`KVar`]).
    pub coeffs: Vec<KVar>,
}

impl RKVar {
    /// A constant `R_K` value (all coefficients pinned).
    pub fn constant(cs: &ConstraintSystem, v: &RingK) -> RKVar {
        RKVar {
            coeffs: v.coeffs().iter().map(|&c| KVar::constant(cs, c)).collect(),
        }
    }

    /// `self + other` (free).
    pub fn add(&self, other: &RKVar) -> RKVar {
        RKVar {
            coeffs: self.coeffs.iter().zip(&other.coeffs).map(|(a, b)| a.add(b)).collect(),
        }
    }

    /// The mixed product `self · rf mod Φ` for `rf ∈ R_F`, mirroring [`RingK::mul_rf`].
    /// Each coefficient product is `K · F` (two base mults).
    pub fn mul_rf(&self, cs: &mut ConstraintSystem, rf: &RFVar) -> RKVar {
        let zero = KVar::constant(cs, Ext2::ZERO);
        let mut prod: Vec<KVar> = vec![zero; PROD_LEN];
        for i in 0..D {
            for j in 0..D {
                let term = kvar_scale_by_var(cs, &self.coeffs[i], &rf.coeffs[j]);
                prod[i + j] = prod[i + j].add(&term);
            }
        }
        RKVar {
            coeffs: reduce(&mut prod, |a, b| a.sub(b)),
        }
    }

    /// Impose `self = other` coefficient-wise.
    pub fn assert_eq(&self, cs: &mut ConstraintSystem, other: &RKVar) {
        for (a, b) in self.coeffs.iter().zip(&other.coeffs) {
            a.assert_eq(cs, b);
        }
    }
}

/// `(field-variable s) · (K value y)` — scale each component of `y` by the variable `s`.
/// Unlike [`KVar::scale`] (constant), `s` here is a wire, so this costs two base mults.
fn kvar_scale_by_var(cs: &mut ConstraintSystem, y: &KVar, s: &Lc) -> KVar {
    KVar {
        c0: Lc::from_var(cs.mul(&y.c0, s)),
        c1: Lc::from_var(cs.mul(&y.c1, s)),
    }
}

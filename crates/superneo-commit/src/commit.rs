//! The Ajtai commitment value and its `R_F`-module operations (Definition 4).
//!
//! A commitment is a vector in `R_F^κ` (the codomain `C` of the homomorphism
//! `L = Commit(A, ·)`). It supports the linear operations the folding reductions
//! need: addition, scaling by a ring element (Π_RLC), and scaling by a field element
//! (Π_DEC's `b^{i-1}` weights).

use superneo_field::fp::Fp;
use superneo_ring::RingElem;

/// An Ajtai commitment `c = A·z ∈ R_F^κ`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Commitment(pub Vec<RingElem>);

impl Commitment {
    /// The all-zero commitment of rank `kappa`.
    pub fn zero(kappa: usize) -> Self {
        Commitment(vec![RingElem::ZERO; kappa])
    }

    /// The module rank `κ`.
    pub fn dim(&self) -> usize {
        self.0.len()
    }

    /// Component-wise sum `self + other` (panics on rank mismatch).
    pub fn add(&self, other: &Commitment) -> Commitment {
        assert_eq!(self.dim(), other.dim(), "commitment rank mismatch");
        Commitment(
            self.0
                .iter()
                .zip(other.0.iter())
                .map(|(a, b)| *a + *b)
                .collect(),
        )
    }

    /// Scale by a ring element `ρ ∈ R_F` (used by Π_RLC).
    pub fn scale_ring(&self, rho: &RingElem) -> Commitment {
        Commitment(self.0.iter().map(|c| *c * *rho).collect())
    }

    /// Scale by a base-field element `s ∈ F` (used by Π_DEC's `b^{i-1}` weights).
    pub fn scale_fp(&self, s: Fp) -> Commitment {
        Commitment(self.0.iter().map(|c| c.scale(s)).collect())
    }
}

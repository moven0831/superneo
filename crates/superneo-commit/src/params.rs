//! SuperNeo parameter sets and their validation guard (Definition 14, Appendix B.2/D.8).

use crate::error::CommitError;

/// A SuperNeo parameter set (the security-relevant constants).
///
/// Logical witness dimensions (`m`, `n_R`, instance counts per call) are supplied
/// separately; these are the fixed cryptographic constants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Params {
    /// Ring degree `d = φ(η)`.
    pub d: usize,
    /// Module rank `κ` of the Ajtai matrix (number of commitment slots).
    pub kappa: usize,
    /// Decomposition base `b`.
    pub b: u64,
    /// Decomposition depth `k` (also the number of carried CE instances; `B = b^k`).
    pub k: usize,
    /// High norm bound `B = b^k`.
    pub big_b: u64,
    /// Expansion factor `T` of the strong sampling set (Theorem 11).
    pub t_exp: u64,
    /// Maximum number of fresh CCS instances folded per step (`K`).
    pub max_fresh: u64,
    /// Soundness parameter `λ` (bits).
    pub lambda: u32,
}

impl Params {
    /// The Goldilocks parameter set of Appendix B.2:
    /// `η = 81`, `Φ = X^54+X^27+1`, `d = 54`, `κ = 18`, `b = 2`, `k = 14`,
    /// `B = 2^14`, `T = 216`, `K ≤ 61`, `λ ≈ 125`.
    pub const GOLDILOCKS_B2: Params = Params {
        d: 54,
        kappa: 18,
        b: 2,
        k: 14,
        big_b: 1 << 14,
        t_exp: 216,
        max_fresh: 61,
        lambda: 125,
    };

    /// Validate the parameter set against the paper's relations.
    ///
    /// Checks `B = b^k` and the folding-admissibility guard
    /// `(K + k)·T·(b − 1) < B` (Definition 14; printed as `True` by the Appendix D.8
    /// Sage script). Returns [`CommitError::GuardFailed`] / [`CommitError::DimMismatch`].
    pub fn validate(&self) -> Result<(), CommitError> {
        // B must equal b^k.
        let pow = self
            .b
            .checked_pow(self.k as u32)
            .ok_or_else(|| CommitError::DimMismatch("b^k overflowed u64".into()))?;
        if pow != self.big_b {
            return Err(CommitError::DimMismatch(format!(
                "B = {} != b^k = {}",
                self.big_b, pow
            )));
        }
        // Folding guard: (K + k)·T·(b − 1) < B.
        let lhs = (self.max_fresh + self.k as u64)
            .saturating_mul(self.t_exp)
            .saturating_mul(self.b - 1);
        if lhs >= self.big_b {
            return Err(CommitError::GuardFailed {
                lhs,
                rhs: self.big_b,
            });
        }
        Ok(())
    }

    /// The maximum ℓ∞ coefficient magnitude allowed in a strong-sampling-set element
    /// `C` (Definition 17). For the Goldilocks set, coefficients lie in `{−2,…,2}`.
    pub const fn challenge_coeff_bound(&self) -> i64 {
        2
    }
}

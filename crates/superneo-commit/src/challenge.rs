//! The strong sampling set `C` (Definition 17) and its sampler.
//!
//! For the Goldilocks parameter set, `C` is the set of ring elements whose
//! coefficients lie in `{−2, −1, 0, 1, 2}`. Differences of distinct elements then
//! have ℓ∞ norm `≤ 4 ≪ b_inv ≈ 2.5·10⁹`, so they are invertible (Theorem 10) — the
//! property the extractor in Π_RLC's soundness proof relies on.

use rand::RngCore;
use superneo_field::fp::Fp;
use superneo_ring::{RingElem, D};

use crate::params::Params;

/// Sample one challenge from `C`: each coefficient uniform in `{−bound, …, bound}`.
pub fn sample_challenge<R: RngCore>(rng: &mut R, params: &Params) -> RingElem {
    let bound = params.challenge_coeff_bound();
    let span = (2 * bound + 1) as u64;
    let mut c = [Fp::ZERO; D];
    for ci in c.iter_mut() {
        // Uniform in {−bound, …, bound}; modulo bias over a 5-way split is negligible.
        let digit = (rng.next_u64() % span) as i64 - bound;
        *ci = Fp::from_i64(digit);
    }
    RingElem::from_coeffs(c)
}

/// Sample `n` independent challenges from `C`.
pub fn sample_challenges<R: RngCore>(rng: &mut R, params: &Params, n: usize) -> Vec<RingElem> {
    (0..n).map(|_| sample_challenge(rng, params)).collect()
}

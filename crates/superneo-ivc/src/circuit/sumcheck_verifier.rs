//! The in-circuit sum-check verifier (M7) — the dominant component of `verify_fold`.
//!
//! Given a sum-check transcript (round polynomials and challenges, supplied as advice),
//! the gadget enforces the verifier's chain over `K`: each round's `g(0) + g(1)` equals
//! the running claim, and the next claim is `g` interpolated at the round challenge. The
//! returned [`KVar`] is the final reduced claim, which a full recursive verifier ties to
//! `Q(r′)` (reconstructed from the claimed evaluations — the documented residual).
//!
//! The round challenges are advice here: deriving them in-circuit is the Fiat–Shamir
//! (Blake3) hashing that is left native in this PoC.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;

use super::cs::ConstraintSystem;
use super::gadgets::KVar;

/// In-circuit Lagrange interpolation of `g` (evaluations at `0..g.len()`) at `r`.
fn k_lagrange(cs: &mut ConstraintSystem, g: &[KVar], r: &KVar) -> KVar {
    let n = g.len();
    let mut acc = KVar::constant(cs, Ext2::ZERO);
    for j in 0..n {
        // numerator Π_{l≠j} (r − l) and denominator Π_{l≠j} (j − l) (a constant).
        let factors: Vec<KVar> = (0..n)
            .filter(|&l| l != j)
            .map(|l| r.sub(&KVar::constant(cs, Ext2::from_base(Fp::new(l as u64)))))
            .collect();
        let mut num = factors[0].clone();
        for f in &factors[1..] {
            num = num.mul(cs, f);
        }
        let mut den = Fp::ONE;
        for l in 0..n {
            if l != j {
                den *= Fp::from_i64(j as i64 - l as i64);
            }
        }
        let inv_den = den.inv().expect("distinct interpolation points");
        let term = g[j].mul(cs, &num).scale(inv_den);
        acc = acc.add(&term);
    }
    acc
}

/// Synthesize the sum-check verifier into `cs`, returning the final reduced claim.
///
/// `init_claim` is the claimed sum `T`; `round_evals[i]` are the degree-`d` round
/// polynomial's evaluations at `0..=d`; `challenges[i]` is round `i`'s challenge.
pub fn synthesize_sumcheck_verifier(
    cs: &mut ConstraintSystem,
    init_claim: Ext2,
    round_evals: &[Vec<Ext2>],
    challenges: &[Ext2],
) -> KVar {
    assert_eq!(round_evals.len(), challenges.len());
    let mut claim = KVar::alloc(cs, init_claim);
    for (g_evals, &r) in round_evals.iter().zip(challenges.iter()) {
        let g: Vec<KVar> = g_evals.iter().map(|&e| KVar::alloc(cs, e)).collect();
        let rk = KVar::alloc(cs, r);
        // g(0) + g(1) == claim.
        g[0].add(&g[1]).assert_eq(cs, &claim);
        // claim ← g(r).
        claim = k_lagrange(cs, &g, &rk);
    }
    claim
}

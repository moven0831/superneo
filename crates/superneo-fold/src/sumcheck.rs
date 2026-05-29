//! The sum-check protocol over the extension field `K` (Definition 6).
//!
//! Generic and oracle-free: the prover is given a set of multilinear tables (each of
//! length `2^ℓ`) and a `combine` closure that forms the summand from one value per
//! table at a point. The combined polynomial `Q(X) = combine(table_1(X), …)` has a
//! known per-variable `degree`; each round sends its evaluations at `0..=degree`.
//! The verifier is stateless given the proof, the claimed sum, and the degree.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;

use crate::error::FoldError;
use crate::transcript::Transcript;

/// A round polynomial, represented by its evaluations at `0, 1, …, degree`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RoundPoly(pub Vec<Ext2>);

impl RoundPoly {
    /// `g(0) + g(1)`, the sum-check round invariant.
    pub fn sum_over_bit(&self) -> Ext2 {
        self.0[0] + self.0[1]
    }

    /// Evaluate at an arbitrary `r ∈ K` by Lagrange interpolation through `0..=deg`.
    pub fn eval_at(&self, r: Ext2) -> Ext2 {
        lagrange_eval(&self.0, r)
    }
}

/// A sum-check proof: one round polynomial per variable.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SumCheckProof {
    /// One [`RoundPoly`] per variable (`ell` of them).
    pub rounds: Vec<RoundPoly>,
}

/// Lagrange interpolation: evaluate at `r` given values at integer points `0..len`.
pub fn lagrange_eval(evals: &[Ext2], r: Ext2) -> Ext2 {
    let n = evals.len();
    let pts: Vec<Ext2> = (0..n).map(|i| Ext2::from_base(Fp::new(i as u64))).collect();
    let mut acc = Ext2::ZERO;
    for i in 0..n {
        let mut num = Ext2::ONE;
        let mut den = Ext2::ONE;
        for j in 0..n {
            if i == j {
                continue;
            }
            num *= r - pts[j];
            den *= pts[i] - pts[j];
        }
        acc += evals[i] * num * den.inv().expect("distinct interpolation points");
    }
    acc
}

/// Prove `Σ_{x ∈ {0,1}^ℓ} combine(tables(x)) = T`.
///
/// Returns the proof and the evaluation point `r' ∈ K^ℓ` fixed by the transcript.
/// Each table is consumed and folded in place. (`T` is not needed by the prover —
/// it is the verifier's input — but equals `proof.rounds[0].sum_over_bit()`.)
pub fn sumcheck_prove(
    mut tables: Vec<Vec<Ext2>>,
    combine: impl Fn(&[Ext2]) -> Ext2,
    degree: usize,
    ell: usize,
    tr: &mut Transcript,
) -> (SumCheckProof, Vec<Ext2>) {
    let n_tables = tables.len();
    let mut size = 1usize << ell;
    debug_assert!(tables.iter().all(|t| t.len() == size));
    let mut rounds = Vec::with_capacity(ell);
    let mut r_prime = Vec::with_capacity(ell);
    let mut vals = vec![Ext2::ZERO; n_tables];

    for _ in 0..ell {
        let half = size / 2;
        let mut g = vec![Ext2::ZERO; degree + 1];
        for (c, gc) in g.iter_mut().enumerate() {
            let pt = Ext2::from_base(Fp::new(c as u64));
            let one_minus = Ext2::ONE - pt;
            let mut acc = Ext2::ZERO;
            for x in 0..half {
                for (ti, t) in tables.iter().enumerate() {
                    vals[ti] = t[x] * one_minus + t[half + x] * pt;
                }
                acc += combine(&vals);
            }
            *gc = acc;
        }
        for &e in &g {
            tr.absorb_ext(b"sc/g", e);
        }
        let r = tr.challenge_ext(b"sc/r");
        for t in tables.iter_mut() {
            for x in 0..half {
                t[x] = t[x] * (Ext2::ONE - r) + t[half + x] * r;
            }
            t.truncate(half);
        }
        size = half;
        rounds.push(RoundPoly(g));
        r_prime.push(r);
    }

    (SumCheckProof { rounds }, r_prime)
}

/// Verify a sum-check proof against the claimed sum `T`.
///
/// Returns the transcript-fixed point `r'` and the final reduced claim `Q(r')`,
/// which the caller must independently match against the protocol's own evaluation.
pub fn sumcheck_verify(
    proof: &SumCheckProof,
    claimed_sum: Ext2,
    degree: usize,
    ell: usize,
    tr: &mut Transcript,
) -> Result<(Vec<Ext2>, Ext2), FoldError> {
    if proof.rounds.len() != ell {
        return Err(FoldError::Malformed(format!(
            "sum-check has {} rounds, expected {ell}",
            proof.rounds.len()
        )));
    }
    let mut claim = claimed_sum;
    let mut r_prime = Vec::with_capacity(ell);
    for (round, poly) in proof.rounds.iter().enumerate() {
        if poly.0.len() != degree + 1 {
            return Err(FoldError::Malformed(format!(
                "round {round} polynomial has {} evals, expected {}",
                poly.0.len(),
                degree + 1
            )));
        }
        if poly.sum_over_bit() != claim {
            return Err(FoldError::SumCheck { round });
        }
        for &e in &poly.0 {
            tr.absorb_ext(b"sc/g", e);
        }
        let r = tr.challenge_ext(b"sc/r");
        claim = poly.eval_at(r);
        r_prime.push(r);
    }
    Ok((r_prime, claim))
}

//! In-circuit Π_CCS verifier (M7, Phase P1): the final `Q(r')` reconstruction.
//!
//! The native verifier ([`superneo_fold::pi_ccs::pi_ccs_verify`], lines 353-381) accepts
//! iff the sum-check's reduced claim equals the combined oracle reconstructed from the
//! claimed evaluations:
//!
//! ```text
//! Q(r') = eq(r',α)·(F + γ^K·NC) + γ^{2K+k}·Eval,   F   = Σ_{i<K}   γ^i  f(ct(y'_{i,·})),
//!                                                  NC  = Σ_{i<K+k} γ^i  Π_{|j|<b}(ct(y'_{i,0}) − j),
//!                                                  Eval= eq(r',r)·Σ_{ci,j,ℓ} γ^{I} cf(y'_{K+ci,j})_ℓ.
//! ```
//!
//! This module synthesizes that reconstruction over [`KVar`] (the `K`-extension gadget)
//! and asserts the reconstructed `Q(r')` equals the wire the [`sumcheck verifier`] reduces
//! to — replacing the advice `expected_final` constant the standalone sum-check circuit
//! used. The sum-check point `r'` is the *same* set of challenge wires the chain consumes
//! (see [`synthesize_sumcheck_verifier`]), so the chain and the tie are bound to one `r'`.
//!
//! Scope (documented residual): the challenges `α, γ, r'` are advice here — the in-circuit
//! Fiat–Shamir that would derive them is a later phase. Feeding a circuit the genuine
//! [`PiCcsTrace`] therefore reproduces the native verifier's decision exactly, which is
//! what the tests assert (tamper an eval or a round polynomial → both reject).
//!
//! [`sumcheck verifier`]: super::sumcheck_verifier
//! [`synthesize_sumcheck_verifier`]: super::sumcheck_verifier::synthesize_sumcheck_verifier

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::{CcsInstance, CcsStructure, CeInstance, GlobalParams, PiCcsProof, PiCcsTrace};
use superneo_fold::SparsePoly;

use super::cs::ConstraintSystem;
use super::gadgets::KVar;
use super::sumcheck_verifier::synthesize_sumcheck_verifier;

/// Highest γ-power index used by `Q(r')`, mirroring `pi_ccs::max_gamma_index`.
fn max_gamma_index(cap_k: usize, k: usize, t: usize, d: usize) -> usize {
    let mut hi = (2 * cap_k + k).max(cap_k + k);
    if k > 0 {
        hi = hi.max(k * t * d - 1);
    }
    hi
}

/// The chain `γ^0, …, γ^hi` as wires, with `γ^0 = 1` constant and `γ^i = γ^{i-1}·γ`
/// enforced (one Karatsuba mul per power). The chain length `k·t·d` dominates the cost.
fn gamma_power_chain(cs: &mut ConstraintSystem, gamma: &KVar, hi: usize) -> Vec<KVar> {
    let mut p = Vec::with_capacity(hi + 1);
    p.push(KVar::constant(cs, Ext2::ONE));
    for i in 1..=hi {
        let prev = p[i - 1].clone();
        p.push(prev.mul(cs, gamma));
    }
    p
}

/// In-circuit `eq(a,b) = Π_i (a_i·b_i + (1−a_i)(1−b_i))`, mirroring `multilinear::eq`.
fn eq_circuit(cs: &mut ConstraintSystem, a: &[KVar], b: &[KVar]) -> KVar {
    assert_eq!(a.len(), b.len());
    let one = KVar::constant(cs, Ext2::ONE);
    let mut acc: Option<KVar> = None;
    for (ai, bi) in a.iter().zip(b.iter()) {
        let ab = ai.mul(cs, bi);
        let lo = one.sub(ai);
        let hi = one.sub(bi);
        let comp = lo.mul(cs, &hi);
        let factor = ab.add(&comp);
        acc = Some(match acc {
            None => factor,
            Some(prev) => prev.mul(cs, &factor),
        });
    }
    acc.unwrap_or(one)
}

/// In-circuit `Π_{j=-(b-1)}^{b-1} (z − j)`, mirroring `pi_ccs::range_product`.
fn range_product_circuit(cs: &mut ConstraintSystem, z: &KVar, b: u64) -> KVar {
    let bi = b as i64;
    let mut acc: Option<KVar> = None;
    for j in -(bi - 1)..=(bi - 1) {
        let factor = z.sub(&KVar::constant(cs, Ext2::from_base(Fp::from_i64(j))));
        acc = Some(match acc {
            None => factor,
            Some(prev) => prev.mul(cs, &factor),
        });
    }
    acc.unwrap_or_else(|| KVar::constant(cs, Ext2::ONE))
}

/// In-circuit `f(args) = Σ_terms coeff·Π_v args[v]^{exp_v}`, mirroring `SparsePoly::eval_ext`.
fn eval_f_circuit(cs: &mut ConstraintSystem, f: &SparsePoly, args: &[KVar]) -> KVar {
    let mut acc = KVar::constant(cs, Ext2::ZERO);
    for (coeff, exps) in &f.terms {
        let mut prod: Option<KVar> = None;
        for (v, &e) in exps.iter().enumerate() {
            for _ in 0..e {
                prod = Some(match prod {
                    None => args[v].clone(),
                    Some(p) => p.mul(cs, &args[v]),
                });
            }
        }
        let term = match prod {
            None => KVar::constant(cs, Ext2::from_base(*coeff)), // empty product = constant term
            Some(p) => p.scale(*coeff),
        };
        acc = acc.add(&term);
    }
    acc
}

/// Build the in-circuit Π_CCS verifier from a real proof and its verification trace.
///
/// Synthesizes (1) the sum-check verifier chain from the public claimed sum `T` to the
/// reduced claim, then (2) the `Q(r')` reconstruction over the claimed evaluations, and
/// asserts the two are equal. The system [`is_satisfied`](ConstraintSystem::is_satisfied)
/// iff [`pi_ccs_verify`](superneo_fold::pi_ccs::pi_ccs_verify) accepts the same proof.
///
/// `trace` must be the [`PiCcsTrace`] produced by
/// [`pi_ccs_verify_traced`](superneo_fold::pi_ccs::pi_ccs_verify_traced) on `(s, fresh,
/// carried, proof)` — its `α, γ, r'` are taken as advice (in-circuit Fiat–Shamir is a
/// later phase). `claimed_sum` and the carried evaluation point are public inputs and are
/// pinned as constants.
pub fn build_pi_ccs_verifier_circuit(
    gp: &GlobalParams,
    s: &CcsStructure,
    fresh: &[CcsInstance],
    carried: &[CeInstance],
    proof: &PiCcsProof,
    trace: &PiCcsTrace,
) -> ConstraintSystem {
    let cap_k = fresh.len();
    let k = carried.len();
    let t = s.t();
    let d = gp.d;

    let mut cs = ConstraintSystem::new();

    // (1) Sum-check chain: public T → reduced claim, exposing the r' wires the tie reuses.
    let round_evals: Vec<Vec<Ext2>> = proof.sumcheck.rounds.iter().map(|p| p.0.clone()).collect();
    let (claim, r_wires) =
        synthesize_sumcheck_verifier(&mut cs, trace.claimed_sum, &round_evals, &trace.r_prime);

    // (2) Challenge advice and the γ-power chain.
    let alpha: Vec<KVar> = trace.alpha.iter().map(|&a| KVar::alloc(&mut cs, a)).collect();
    let gamma = KVar::alloc(&mut cs, trace.gamma);
    let gpow = gamma_power_chain(&mut cs, &gamma, max_gamma_index(cap_k, k, t, d));

    // (3) Claimed-eval advice: constant terms ct[i][j], plus the higher coefficients of the
    //     carried instances (l ≥ 1; l = 0 reuses ct) for the Eval term.
    let ct: Vec<Vec<KVar>> = (0..cap_k + k)
        .map(|i| {
            (0..t)
                .map(|j| KVar::alloc(&mut cs, proof.evals[i][j].ct()))
                .collect()
        })
        .collect();
    let carr: Vec<Vec<Vec<KVar>>> = (0..k)
        .map(|ci| {
            (0..t)
                .map(|j| {
                    let mut v = Vec::with_capacity(d);
                    v.push(ct[cap_k + ci][j].clone()); // l = 0
                    for l in 1..d {
                        v.push(KVar::alloc(&mut cs, proof.evals[cap_k + ci][j].coeff(l)));
                    }
                    v
                })
                .collect()
        })
        .collect();

    // (4) F = Σ_{i<K} γ^i · f(ct(y'_{i,·})).
    let mut f_val = KVar::constant(&cs, Ext2::ZERO);
    for i in 0..cap_k {
        let args: Vec<KVar> = (0..t).map(|j| ct[i][j].clone()).collect();
        let fi = eval_f_circuit(&mut cs, &s.f, &args);
        let scaled = gpow[i].mul(&mut cs, &fi);
        f_val = f_val.add(&scaled);
    }

    // (5) NC = Σ_{i<K+k} γ^i · range_product(ct(y'_{i,0})).
    let mut n_val = KVar::constant(&cs, Ext2::ZERO);
    for i in 0..(cap_k + k) {
        let rp = range_product_circuit(&mut cs, &ct[i][0], gp.b);
        let scaled = gpow[i].mul(&mut cs, &rp);
        n_val = n_val.add(&scaled);
    }

    // (6) Eval = eq(r', r)·Σ_{ci,j,ℓ} γ^{ci + k·j + k·t·ℓ} · cf(y'_{K+ci,j})_ℓ.
    let mut e_val = KVar::constant(&cs, Ext2::ZERO);
    if k > 0 {
        let mut inner = KVar::constant(&cs, Ext2::ZERO);
        for (ci, carr_ci) in carr.iter().enumerate() {
            for (j, carr_cij) in carr_ci.iter().enumerate() {
                for (l, coeff) in carr_cij.iter().enumerate() {
                    let idx = ci + k * j + k * t * l;
                    let scaled = gpow[idx].mul(&mut cs, coeff);
                    inner = inner.add(&scaled);
                }
            }
        }
        let cr: Vec<KVar> = carried[0].r.iter().map(|&x| KVar::constant(&cs, x)).collect();
        let eq_cr = eq_circuit(&mut cs, &r_wires, &cr);
        e_val = eq_cr.mul(&mut cs, &inner);
    }

    // (7) Q(r') = eq(r',α)·(F + γ^K·NC) + γ^{2K+k}·Eval.
    let eq_a = eq_circuit(&mut cs, &r_wires, &alpha);
    let nc_scaled = gpow[cap_k].mul(&mut cs, &n_val);
    let bracket = f_val.add(&nc_scaled);
    let lhs = eq_a.mul(&mut cs, &bracket);
    let eval_scaled = gpow[2 * cap_k + k].mul(&mut cs, &e_val);
    let q_rp = lhs.add(&eval_scaled);

    // (8) Tie the reconstructed Q(r') to the sum-check's reduced claim.
    q_rp.assert_eq(&mut cs, &claim);

    cs
}

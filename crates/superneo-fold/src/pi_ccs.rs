//! Π_CCS: the strong interactive reduction `CCS^K × CE^k → CE^{K+k}` (§7.3, Lemma 3).
//!
//! A single sum-check over the combined polynomial
//! `Q(X) = eq(X,α)·(F(X) + γ^K·NC(X)) + γ^{2K+k}·Eval(X)`, where
//!   * `F`  = `Σ_{i∈[K]} γ^{i-1} f(~(M_1 z_i), …, ~(M_t z_i))`   (CCS satisfaction, fresh),
//!   * `NC` = `Σ_{i∈[K+k]} γ^{i-1} Π_{j=-(b-1)}^{b-1}(~z_i − j)`  (norm/range, all),
//!   * `Eval` = `eq(X,r)·Σ_{carried i,j,ℓ} γ^{I(i,j,ℓ)} ~(cf(M̄_j ẑ_i)_ℓ)`  (prior claims).
//!
//! The disjoint γ-power ranges (`F∈[0,K)`, `NC∈[K,2K+k)`, `Eval∈[2K+k, …)`) separate the
//! three checks under Schwartz–Zippel. The claimed sum is
//! `T = γ^{2K+k}·Σ γ^{I} cf(y_{i,j})_ℓ`, chosen so `T = Σ_x Q(x)` for honest inputs.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_ring::RingK;

use crate::error::FoldError;
use crate::multilinear::{eq, eq_table, mle_eval_ring};
use crate::sumcheck::{sumcheck_prove, sumcheck_verify};
use crate::transcript::Transcript;
use crate::types::{
    bar_ring_matvec, embed_witness, field_matvec, flatten_ring, CcsInstance, CcsStructure,
    CcsWitness, CeInstance, CeWitness, GlobalParams, PiCcsProof,
};

/// `Π_{j=-(b-1)}^{b-1} (z − j)` — vanishes iff `‖z‖∞ < b`.
fn range_product(z: Ext2, b: u64) -> Ext2 {
    let bi = b as i64;
    let mut acc = Ext2::ONE;
    for j in -(bi - 1)..=(bi - 1) {
        acc *= z - Ext2::from_base(Fp::from_i64(j));
    }
    acc
}

/// Powers `γ^0, …, γ^{hi}`.
fn gamma_powers(gamma: Ext2, hi: usize) -> Vec<Ext2> {
    let mut p = Vec::with_capacity(hi + 1);
    p.push(Ext2::ONE);
    for i in 1..=hi {
        p.push(p[i - 1] * gamma);
    }
    p
}

/// Highest γ-power index used by the combined oracle, for `(K, k, t, d)`.
fn max_gamma_index(cap_k: usize, k: usize, t: usize, d: usize) -> usize {
    let mut hi = (2 * cap_k + k).max(cap_k + k);
    if k > 0 {
        hi = hi.max(k * t * d - 1);
    }
    hi
}

/// Per-variable degree of the round polynomials.
fn round_degree(cap_k: usize, b: u64, f_degree: usize) -> usize {
    let mut deg = 2 * b as usize; // NC contributes 2b
    if cap_k > 0 {
        deg = deg.max(1 + f_degree); // eq·F contributes 1 + deg(f)
    }
    deg
}

fn absorb_ringk(tr: &mut Transcript, label: &'static [u8], y: &RingK) {
    for c in y.coeffs() {
        tr.absorb_ext(label, *c);
    }
}

/// Bind the relation `s`, absorb all input instances, then squeeze `(α, γ)`. Binding `s`
/// here (both prover and verifier call this shared helper) ensures even a standalone
/// `fold`/`verify_fold` — which builds its own transcript with no `absorb_vk` — commits
/// the challenges to the relation being proven.
fn absorb_inputs_and_challenges(
    tr: &mut Transcript,
    gp: &GlobalParams,
    s: &CcsStructure,
    fresh: &[CcsInstance],
    carried: &[CeInstance],
) -> (Vec<Ext2>, Ext2) {
    tr.absorb_structure(s);
    for inst in fresh {
        tr.absorb_commitment(b"piccs/ccs", &inst.c);
    }
    for inst in carried {
        tr.absorb_commitment(b"piccs/ce", &inst.c);
        for &ri in &inst.r {
            tr.absorb_ext(b"piccs/r", ri);
        }
        for yj in &inst.y {
            absorb_ringk(tr, b"piccs/y", yj);
        }
    }
    let alpha = tr.challenge_ext_vec(b"piccs/alpha", gp.log_m);
    let gamma = tr.challenge_ext(b"piccs/gamma");
    (alpha, gamma)
}

/// The verifier's claimed sum `T = γ^{2K+k}·Σ γ^{I(i,j,ℓ)} cf(y_{i,j})_ℓ`.
fn claimed_sum(carried: &[CeInstance], gp: &GlobalParams, gpow: &[Ext2], cap_k: usize) -> Ext2 {
    let (k, d) = (carried.len(), gp.d);
    if k == 0 {
        return Ext2::ZERO;
    }
    let t = carried[0].y.len();
    let mut inner = Ext2::ZERO;
    for (ci, inst) in carried.iter().enumerate() {
        for (j, yj) in inst.y.iter().enumerate() {
            for l in 0..d {
                inner += gpow[ci + k * j + k * t * l] * yj.coeff(l);
            }
        }
    }
    gpow[2 * cap_k + k] * inner
}

/// Pad a length-`m` field column to `cap_m`, as extension values.
fn padded_ext(col: &[Fp], cap_m: usize) -> Vec<Ext2> {
    let mut v = vec![Ext2::ZERO; cap_m];
    for (i, &x) in col.iter().enumerate() {
        v[i] = Ext2::from_base(x);
    }
    v
}

/// Build the sum-check tables and record the layout offsets.
struct Layout {
    eq_r: usize,
    u0: usize,
    z0: usize,
    v0: usize,
}

#[allow(clippy::too_many_arguments)]
fn build_tables(
    gp: &GlobalParams,
    s: &CcsStructure,
    alpha: &[Ext2],
    carried_r: Option<&[Ext2]>,
    fresh_z: &[Vec<Fp>],                         // K field witnesses (length m)
    all_z: &[Vec<Fp>],                           // K+k field witnesses (length m), for NC
    carried_zr: &[Vec<superneo_ring::RingElem>], // k ring witnesses, for Eval
) -> (Vec<Vec<Ext2>>, Layout) {
    let cap_m = gp.cap_m();
    let t = s.t();
    let d = gp.d;
    let cap_k = fresh_z.len();
    let k = carried_zr.len();

    let mut tables: Vec<Vec<Ext2>> = Vec::new();
    // 0: eq_α
    tables.push(eq_table(alpha));
    // 1: eq_r
    let eq_r_idx = tables.len();
    match carried_r {
        Some(r) => tables.push(eq_table(r)),
        None => tables.push(vec![Ext2::ZERO; cap_m]),
    }
    // u[i][j] = M_j z_i (fresh)
    let u0 = tables.len();
    for zi in fresh_z {
        for mat in &s.matrices {
            tables.push(padded_ext(&field_matvec(mat, zi), cap_m));
        }
    }
    // zt[i] = z_i (all)
    let z0 = tables.len();
    for zi in all_z {
        tables.push(padded_ext(zi, cap_m));
    }
    // v[ci][j][ℓ] = cf(M̄_j ẑ_i)_ℓ (carried)
    let v0 = tables.len();
    for zr in carried_zr {
        for mat in &s.matrices {
            let ring_vec = bar_ring_matvec(mat, zr); // length m
            for l in 0..d {
                let col: Vec<Fp> = ring_vec.iter().map(|e| e.coeffs()[l]).collect();
                tables.push(padded_ext(&col, cap_m));
            }
        }
    }
    debug_assert_eq!(tables.len(), v0 + k * t * d);
    let _ = cap_k;
    (
        tables,
        Layout {
            eq_r: eq_r_idx,
            u0,
            z0,
            v0,
        },
    )
}

/// The combined-oracle summand `Q` evaluated from one value per table.
#[allow(clippy::too_many_arguments)]
fn combine(
    vals: &[Ext2],
    lay: &Layout,
    f: &crate::types::SparsePoly,
    gpow: &[Ext2],
    cap_k: usize,
    k: usize,
    t: usize,
    d: usize,
    b: u64,
) -> Ext2 {
    let eq_a = vals[0];
    let eq_r = vals[lay.eq_r];

    let mut f_val = Ext2::ZERO;
    for i in 0..cap_k {
        let args: Vec<Ext2> = (0..t).map(|j| vals[lay.u0 + i * t + j]).collect();
        f_val += gpow[i] * f.eval_ext(&args);
    }
    let mut nc = Ext2::ZERO;
    for i in 0..(cap_k + k) {
        nc += gpow[i] * range_product(vals[lay.z0 + i], b);
    }
    let mut ev = Ext2::ZERO;
    for ci in 0..k {
        for j in 0..t {
            for l in 0..d {
                ev += gpow[ci + k * j + k * t * l] * vals[lay.v0 + (ci * t + j) * d + l];
            }
        }
    }
    eq_a * (f_val + gpow[cap_k] * nc) + gpow[2 * cap_k + k] * eq_r * ev
}

/// Assemble the `K+k` output CE instances at point `r'` with claimed evals.
fn output_instances(
    fresh: &[CcsInstance],
    carried: &[CeInstance],
    r_prime: &[Ext2],
    evals: &[Vec<RingK>],
) -> Vec<CeInstance> {
    let cap_k = fresh.len();
    (0..fresh.len() + carried.len())
        .map(|i| {
            let c = if i < cap_k {
                fresh[i].c.clone()
            } else {
                carried[i - cap_k].c.clone()
            };
            CeInstance {
                c,
                r: r_prime.to_vec(),
                y: evals[i].clone(),
            }
        })
        .collect()
}

/// Π_CCS prover. Returns the `K+k` output CE instances, their witnesses, and the proof.
pub fn pi_ccs_prove(
    tr: &mut Transcript,
    gp: &GlobalParams,
    s: &CcsStructure,
    fresh: &[CcsInstance],
    fresh_wit: &[CcsWitness],
    carried: &[CeInstance],
    carried_wit: &[CeWitness],
) -> (Vec<CeInstance>, Vec<CeWitness>, PiCcsProof) {
    let cap_k = fresh.len();
    let k = carried.len();
    let t = s.t();
    let d = gp.d;

    let (alpha, gamma) = absorb_inputs_and_challenges(tr, gp, s, fresh, carried);
    let gpow = gamma_powers(gamma, max_gamma_index(cap_k, k, t, d));
    let degree = round_degree(cap_k, gp.b, s.f.degree());

    // Witnesses: ring vectors ẑ_i and field vectors z_i (length m), for all K+k.
    let fresh_zr: Vec<Vec<_>> = fresh_wit.iter().map(|w| embed_witness(&w.z)).collect();
    let carried_zr: Vec<Vec<_>> = carried_wit.iter().map(|w| w.z_ring.clone()).collect();
    let fresh_z: Vec<Vec<Fp>> = fresh_wit.iter().map(|w| w.z.clone()).collect();
    let mut all_z = fresh_z.clone();
    all_z.extend(carried_zr.iter().map(|zr| flatten_ring(zr)));
    let all_zr: Vec<Vec<_>> = fresh_zr.iter().chain(carried_zr.iter()).cloned().collect();

    let carried_r = if k > 0 { Some(&carried[0].r[..]) } else { None };
    let (tables, lay) = build_tables(gp, s, &alpha, carried_r, &fresh_z, &all_z, &carried_zr);

    let comb = |vals: &[Ext2]| combine(vals, &lay, &s.f, &gpow, cap_k, k, t, d, gp.b);

    let (sc_proof, r_prime) = sumcheck_prove(tables, comb, degree, gp.log_m, tr);

    // Claimed evaluations y'_{i,j} = ~(M̄_j ẑ_i)(r').
    let evals: Vec<Vec<RingK>> = all_zr
        .iter()
        .map(|zr| {
            s.matrices
                .iter()
                .map(|mat| mle_eval_ring(&bar_ring_matvec(mat, zr), &r_prime))
                .collect()
        })
        .collect();
    for row in &evals {
        for y in row {
            absorb_ringk(tr, b"piccs/yprime", y);
        }
    }

    let out_inst = output_instances(fresh, carried, &r_prime, &evals);
    let out_wit: Vec<CeWitness> = all_zr
        .into_iter()
        .map(|zr| CeWitness { z_ring: zr })
        .collect();
    (
        out_inst,
        out_wit,
        PiCcsProof {
            sumcheck: sc_proof,
            evals,
        },
    )
}

/// Π_CCS verifier. Returns the `K+k` output CE instances.
pub fn pi_ccs_verify(
    tr: &mut Transcript,
    gp: &GlobalParams,
    s: &CcsStructure,
    fresh: &[CcsInstance],
    carried: &[CeInstance],
    proof: &PiCcsProof,
) -> Result<Vec<CeInstance>, FoldError> {
    let cap_k = fresh.len();
    let k = carried.len();
    let t = s.t();
    let d = gp.d;

    let (alpha, gamma) = absorb_inputs_and_challenges(tr, gp, fresh, carried);
    let gpow = gamma_powers(gamma, max_gamma_index(cap_k, k, t, d));
    let degree = round_degree(cap_k, gp.b, s.f.degree());
    let t_claimed = claimed_sum(carried, gp, &gpow, cap_k);

    let (r_prime, claim) = sumcheck_verify(&proof.sumcheck, t_claimed, degree, gp.log_m, tr)?;

    // Shape check on the claimed evaluations.
    if proof.evals.len() != cap_k + k || proof.evals.iter().any(|row| row.len() != t) {
        return Err(FoldError::Malformed(
            "Π_CCS eval matrix has wrong shape".into(),
        ));
    }
    for row in &proof.evals {
        for y in row {
            absorb_ringk(tr, b"piccs/yprime", y);
        }
    }

    // Reconstruct Q(r') from the claimed evals.
    let mut f_val = Ext2::ZERO;
    for i in 0..cap_k {
        let args: Vec<Ext2> = (0..t).map(|j| proof.evals[i][j].ct()).collect();
        f_val += gpow[i] * s.f.eval_ext(&args);
    }
    let mut n_val = Ext2::ZERO;
    for i in 0..(cap_k + k) {
        n_val += gpow[i] * range_product(proof.evals[i][0].ct(), gp.b);
    }
    let mut e_val = Ext2::ZERO;
    if k > 0 {
        let mut inner = Ext2::ZERO;
        for ci in 0..k {
            for j in 0..t {
                for l in 0..d {
                    inner += gpow[ci + k * j + k * t * l] * proof.evals[cap_k + ci][j].coeff(l);
                }
            }
        }
        e_val = eq(&r_prime, &carried[0].r) * inner;
    }
    let q_rp = eq(&r_prime, &alpha) * (f_val + gpow[cap_k] * n_val) + gpow[2 * cap_k + k] * e_val;

    if q_rp != claim {
        return Err(FoldError::Reduction(
            "Π_CCS final Q(r') check failed".into(),
        ));
    }

    Ok(output_instances(fresh, carried, &r_prime, &proof.evals))
}

//! The Spartan-style reduction compressing the final SuperNeo accumulator (M6).
//!
//! The accumulator is `k` CE instances `(c_i, r_i, y_i)` with witnesses `ẑ_i`. Their
//! field witnesses are stacked into one multilinear `Z` (instance `i` occupying block
//! `[i·cap_m, (i+1)·cap_m)`, rows little-endian) and committed once with the BaseFold
//! PCS. The CE relation (Definition 13) is then reduced to evaluations of `Z̃`:
//!
//!   * **Ajtai opening** `c_i = A·ẑ_i` — via `cf(c_i[κ]) = Σ_col rot(A[κ][col])·z_i`
//!     (the rotation identity), one linear claim per output coefficient;
//!   * **evaluation claims** `ct(y_{i,j}) = ~(M_j z_i)(r_i) = Σ_col M̃_j(r_i,col)·z_i[col]`
//!     (Theorem 6), one linear claim per matrix;
//!
//! all batched with powers of a challenge `η` into a single inner-product sum-check
//! `Σ_x W(x)·Z(x) = v`, plus a separate **norm** sum-check `Σ_x eq(τ,x)·NC(Z(x)) = 0`
//! enforcing `‖Z‖∞ < b` (`NC(v) = Π_{j=−(b−1)}^{b−1}(v−j)`). Each sum-check reduces to
//! one opening of `Z̃`, verified against the commitment by the PCS.
//!
//! All `D` coefficients of the evaluation claims `y_{i,j} ∈ R_K` are enforced — the
//! constant term (Theorem 6) *and* the higher coefficients (`ℓ > 0`) — via the bar-lifted
//! matrix and the rotation identity, the same linear-claim shape as the Ajtai binding.
//! The Ajtai opening and the ℓ∞ norm bound are likewise enforced.

use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::multilinear::{eq_table, next_pow2};
use superneo_fold::sumcheck::lagrange_eval;
use superneo_fold::types::flatten_ring;
use superneo_fold::{CcsStructure, CeInstance, CeWitness, GlobalParams, Transcript};
use superneo_ring::maps::bar_matrix_row;
use superneo_ring::s_action::rot;
use superneo_ring::D;

use crate::basefold::{self, OpenProof};
use crate::code::log2;
use crate::error::SnarkError;
use crate::mle::{eq_eval, eq_table as eq_table_le, eval as mle_eval_le, fold_evals};

/// A compressed proof of the final accumulator's CE satisfaction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SnarkProof {
    /// BaseFold commitment to the stacked witness `Z`.
    pub z_commit: basefold::Commitment,
    /// Round polynomials of the batched linear-claim (inner-product) sum-check.
    pub lin_rounds: Vec<Vec<Ext2>>,
    /// PCS opening of `Z̃` at the linear sum-check's reduced point.
    pub lin_open: OpenProof,
    /// Round polynomials of the norm sum-check.
    pub norm_rounds: Vec<Vec<Ext2>>,
    /// PCS opening of `Z̃` at the norm sum-check's reduced point.
    pub norm_open: OpenProof,
}

/// Domain-separation label for the compression transcript.
const SNARK_DOMAIN: &[u8] = b"superneo/snark/v1";

// ---- multilinear helpers --------------------------------------------------

/// The range polynomial `NC(v) = Π_{j=−(b−1)}^{b−1}(v − j)`; vanishes iff `‖v‖∞ < b`.
fn range_nc(v: Ext2, b: u64) -> Ext2 {
    let bi = b as i64;
    let mut acc = Ext2::ONE;
    for j in -(bi - 1)..=(bi - 1) {
        acc *= v - Ext2::from_base(Fp::from_i64(j));
    }
    acc
}

// ---- LSB-first sum-check -----------------------------------------------------

/// Prove `Σ_{x∈{0,1}^ell} combine(tables(x)) = T`, folding the LSB first (so the
/// reduced point is directly usable as a PCS evaluation point). Returns the round
/// polynomials (each `degree+1` evaluations) and the reduced point.
fn sumcheck_prove_le(
    mut tables: Vec<Vec<Ext2>>,
    combine: impl Fn(&[Ext2]) -> Ext2,
    degree: usize,
    ell: usize,
    tr: &mut Transcript,
) -> (Vec<Vec<Ext2>>, Vec<Ext2>) {
    let n_tables = tables.len();
    let mut rounds = Vec::with_capacity(ell);
    let mut point = Vec::with_capacity(ell);
    let mut vals = vec![Ext2::ZERO; n_tables];

    for _ in 0..ell {
        let half = tables[0].len() / 2;
        let mut g = vec![Ext2::ZERO; degree + 1];
        for (c, gc) in g.iter_mut().enumerate() {
            let cc = Ext2::from_base(Fp::new(c as u64));
            let mut acc = Ext2::ZERO;
            for j in 0..half {
                for (ti, t) in tables.iter().enumerate() {
                    vals[ti] = t[2 * j] + cc * (t[2 * j + 1] - t[2 * j]);
                }
                acc += combine(&vals);
            }
            *gc = acc;
        }
        for &e in &g {
            tr.absorb_ext(b"snark/sc", e);
        }
        let alpha = tr.challenge_ext(b"snark/scr");
        for t in tables.iter_mut() {
            *t = fold_evals(t, alpha);
        }
        rounds.push(g);
        point.push(alpha);
    }
    (rounds, point)
}

/// Verify an LSB-first sum-check against the claimed sum, returning the reduced point
/// and the final claim `combine(tables(point))`.
fn sumcheck_verify_le(
    rounds: &[Vec<Ext2>],
    claimed_sum: Ext2,
    degree: usize,
    ell: usize,
    tr: &mut Transcript,
) -> Result<(Vec<Ext2>, Ext2), SnarkError> {
    if rounds.len() != ell {
        return Err(SnarkError::Verify("sum-check length".into()));
    }
    let mut claim = claimed_sum;
    let mut point = Vec::with_capacity(ell);
    for g in rounds {
        if g.len() != degree + 1 {
            return Err(SnarkError::Verify("round-poly degree".into()));
        }
        if g[0] + g[1] != claim {
            return Err(SnarkError::Verify("sum-check round mismatch".into()));
        }
        for &e in g {
            tr.absorb_ext(b"snark/sc", e);
        }
        let alpha = tr.challenge_ext(b"snark/scr");
        claim = lagrange_eval(g, alpha);
        point.push(alpha);
    }
    Ok((point, claim))
}

// ---- the public linear-claim weight table (shared by prover & verifier) ------

/// Build the batched linear-claim weight table `W` (length `2^ν`, over `K`) and the
/// target `v = Σ_claims η^c · public_c`, from public data only. Claims, in η-power
/// order per instance: the `t·d` evaluation-claim coefficients (`y_{i,j}.coeff(ℓ)`,
/// `ℓ ∈ [0,d)`), then the `κ·d` Ajtai opening claims (`cf(c_i[κ])_ℓ`).
fn build_linear(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    acc: &[CeInstance],
    eta: Ext2,
    n: usize,
    cap_m: usize,
) -> (Vec<Ext2>, Ext2) {
    let t = s.t();
    let kappa = pp.kappa;

    // Precompute rotation matrices rot(A[κ][col_ring]) once.
    let rots: Vec<Vec<[[Fp; D]; D]>> = (0..kappa)
        .map(|kp| (0..gp.n_r).map(|cr| rot(pp.entry(kp, cr))).collect())
        .collect();

    let mut w = vec![Ext2::ZERO; n];
    let mut v = Ext2::ZERO;
    let mut eta_pow = Ext2::ONE;

    for (i, inst) in acc.iter().enumerate() {
        let base = i * cap_m;
        let eqtab = eq_table(&inst.r); // MSB-first, matching compute_evals

        // Evaluation claims: bind ALL `D` coefficients of `y_{i,j} ∈ R_K` (the constant
        // term of Theorem 6 *and* the higher coefficients), via the bar-lifted matrix and
        // the rotation identity — the same linear-claim shape as the Ajtai block below:
        //   coeff_ℓ(y_{i,j}) = Σ_{cr,blk} [Σ_out eq(r_i,out)·rot(M̄_j[out][cr])[ℓ][blk]] · z[cr·D+blk]
        // (binding only the constant term would let a higher-coefficient forgery pass.)
        for j in 0..t {
            let mat = &s.matrices[j];
            // racc[cr][ℓ][blk] = Σ_out eq(r_i,out)·rot(M̄_j[out][cr])[ℓ][blk]  (over K).
            let mut racc = vec![[[Ext2::ZERO; D]; D]; gp.n_r];
            for (out, eo) in eqtab.iter().enumerate().take(gp.m) {
                let bar_row = bar_matrix_row(&mat[out]).expect("matrix row length is a multiple of d");
                for cr in 0..gp.n_r {
                    let rmat = rot(&bar_row[cr]);
                    for l in 0..D {
                        for blk in 0..D {
                            racc[cr][l][blk] += eo.mul_base(rmat[l][blk]);
                        }
                    }
                }
            }
            for l in 0..D {
                for cr in 0..gp.n_r {
                    for blk in 0..D {
                        w[base + cr * D + blk] += eta_pow * racc[cr][l][blk];
                    }
                }
                // Target uses full K-multiplication: `coeff(ℓ)` is an `Ext2` (unlike the
                // base-field Ajtai coeffs below, which use `mul_base`).
                v += eta_pow * inst.y[j].coeff(l);
                eta_pow *= eta;
            }
        }

        // Ajtai opening claims: weight rot(A[κ][col_ring])[ℓ][col_in_block].
        for kp in 0..kappa {
            for l in 0..D {
                for cr in 0..gp.n_r {
                    let rmat = &rots[kp][cr];
                    for blk in 0..D {
                        let col = cr * D + blk;
                        w[base + col] += eta_pow.mul_base(rmat[l][blk]);
                    }
                }
                v += eta_pow.mul_base(inst.c.0[kp].coeffs()[l]);
                eta_pow *= eta;
            }
        }
    }
    (w, v)
}

/// Stack the accumulator field witnesses into the multilinear `Z` (length `n`).
fn stack_witness(gp: &GlobalParams, acc_wit: &[CeWitness], n: usize, cap_m: usize) -> Vec<Fp> {
    let mut z = vec![Fp::ZERO; n];
    for (i, w) in acc_wit.iter().enumerate() {
        let zi = flatten_ring(&w.z_ring);
        let base = i * cap_m;
        z[base..base + zi.len().min(gp.n_f)].copy_from_slice(&zi[..zi.len().min(gp.n_f)]);
    }
    z
}

/// Absorb the accumulator instances (binding the proof to the public claim).
fn absorb_acc(tr: &mut Transcript, acc: &[CeInstance]) {
    for inst in acc {
        tr.absorb_commitment(b"snark/c", &inst.c);
        for &ri in &inst.r {
            tr.absorb_ext(b"snark/r", ri);
        }
        for yj in &inst.y {
            for c in yj.coeffs() {
                tr.absorb_ext(b"snark/y", *c);
            }
        }
    }
}

// ---- top-level API -----------------------------------------------------------

/// Compress an accumulator (`k` CE instances + their witnesses) into a proof whose size
/// is independent of the IVC length (the number of folds) — constant in the computation,
/// though not smaller than the witness at these PoC parameters (the FRI query phase
/// dominates). The witness itself is never revealed.
pub fn compress(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    acc: &[CeInstance],
    acc_wit: &[CeWitness],
) -> Result<SnarkProof, SnarkError> {
    if acc.len() != acc_wit.len() {
        return Err(SnarkError::Malformed("instance/witness count".into()));
    }
    let cap_m = gp.cap_m();
    let n = next_pow2(next_pow2(acc.len().max(1)) * cap_m);
    let nu = log2(n);

    let z = stack_witness(gp, acc_wit, n, cap_m);
    let (z_commit, z_data) = basefold::commit(&z);

    let mut tr = Transcript::new(SNARK_DOMAIN);
    tr.absorb_vk(pp, s); // bind the commitment key A and the relation s into Fiat–Shamir
    absorb_acc(&mut tr, acc);
    tr.absorb_bytes(b"snark/zroot", &z_commit.root);
    let eta = tr.challenge_ext(b"snark/eta");
    let tau = tr.challenge_ext_vec(b"snark/tau", nu);

    // Batched linear-claim sum-check: Σ_x W(x)·Z(x) = v.
    let (w, _v) = build_linear(gp, pp, s, acc, eta, n, cap_m);
    let z_ext: Vec<Ext2> = z.iter().map(|&x| Ext2::from_base(x)).collect();
    let (lin_rounds, lin_pt) = sumcheck_prove_le(
        vec![w, z_ext.clone()],
        |vals| vals[0] * vals[1],
        2,
        nu,
        &mut tr,
    );
    let lin_open = basefold::open(&mut tr, &z_commit, &z_data, &lin_pt);

    // Norm sum-check: Σ_x eq(τ,x)·NC(Z(x)) = 0.
    let norm_degree = 2 * gp.b as usize;
    let (norm_rounds, norm_pt) = sumcheck_prove_le(
        vec![eq_table_le(&tau), z_ext],
        |vals| vals[0] * range_nc(vals[1], gp.b),
        norm_degree,
        nu,
        &mut tr,
    );
    let norm_open = basefold::open(&mut tr, &z_commit, &z_data, &norm_pt);

    Ok(SnarkProof {
        z_commit,
        lin_rounds,
        lin_open,
        norm_rounds,
        norm_open,
    })
}

/// Verify a compressed accumulator proof.
pub fn verify(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    acc: &[CeInstance],
    proof: &SnarkProof,
) -> Result<(), SnarkError> {
    let cap_m = gp.cap_m();
    let n = next_pow2(next_pow2(acc.len().max(1)) * cap_m);
    let nu = log2(n);
    if proof.z_commit.num_vars != nu {
        return Err(SnarkError::Verify("witness arity mismatch".into()));
    }

    let mut tr = Transcript::new(SNARK_DOMAIN);
    tr.absorb_vk(pp, s); // bind the commitment key A and the relation s into Fiat–Shamir
    absorb_acc(&mut tr, acc);
    tr.absorb_bytes(b"snark/zroot", &proof.z_commit.root);
    let eta = tr.challenge_ext(b"snark/eta");
    let tau = tr.challenge_ext_vec(b"snark/tau", nu);

    // Linear-claim sum-check.
    let (w, v) = build_linear(gp, pp, s, acc, eta, n, cap_m);
    let (lin_pt, lin_claim) = sumcheck_verify_le(&proof.lin_rounds, v, 2, nu, &mut tr)?;
    basefold::verify(&mut tr, &proof.z_commit, &lin_pt, &proof.lin_open)?;
    // final claim must equal W̃(ρ)·Z̃(ρ).
    let w_at = mle_eval_le(&w, &lin_pt);
    if lin_claim != w_at * proof.lin_open.value {
        return Err(SnarkError::Verify("linear-claim final tie failed".into()));
    }

    // Norm sum-check (claimed sum 0).
    let norm_degree = 2 * gp.b as usize;
    let (norm_pt, norm_claim) =
        sumcheck_verify_le(&proof.norm_rounds, Ext2::ZERO, norm_degree, nu, &mut tr)?;
    basefold::verify(&mut tr, &proof.z_commit, &norm_pt, &proof.norm_open)?;
    if norm_claim != eq_eval(&tau, &norm_pt) * range_nc(proof.norm_open.value, gp.b) {
        return Err(SnarkError::Verify("norm final tie failed".into()));
    }

    Ok(())
}

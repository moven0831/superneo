//! The BaseFold multilinear polynomial commitment scheme ([84], [79]) — the
//! post-quantum, hash-based PCS at the heart of the final SNARK (M6).
//!
//! A multilinear polynomial given by its evaluation table over `{0,1}^ν` is committed
//! by Reed–Solomon-encoding its coefficient vector and Merkle-committing the codeword
//! ([`commit`]). An evaluation `f(r) = v` is proved by running a sum-check on
//! `Σ_x eq(r,x)·f(x) = v` whose per-round challenge simultaneously folds the codeword
//! FRI-style ([`open`]); the verifier replays the sum-check, checks that the
//! fully-folded codeword constant equals the reduced claim, and spot-checks the folding
//! with Merkle queries ([`verify`]). Corrupting the codeword, a layer, or a sum-check
//! message makes one of these checks fail.
//!
//! Conventions are little-endian throughout: index `i` has bit `k` as variable `x_k`,
//! the sum-check binds `x_0` (the LSB) first, and the codeword folds low/high halves —
//! the two are tied by the [`code`](crate::code) commutation invariant.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;

use crate::code::{self, fold_codeword, root_of_unity};
use crate::error::SnarkError;
use crate::merkle::{verify_path, MerklePath, MerkleTree};
use crate::mle::{self, eq_eval, fold_evals};
use superneo_fold::Transcript;

/// Reed–Solomon blow-up exponent (rate `ρ = 2^{-LOG_BLOWUP}`).
const LOG_BLOWUP: usize = 2;
/// The final (degree-0) codeword length `2^{LOG_BLOWUP}`.
const FINAL_LEN: usize = 1 << LOG_BLOWUP;
/// Target soundness of the FRI query phase, in bits.
const SECURITY_BITS: usize = 100;
/// Number of FRI consistency queries, derived from the target soundness and the rate.
///
/// Each query rejects a word `δ`-far from the code with probability `≥ δ`. Under the
/// proximity-gap regime (`δ → 1 − ρ`, the standard deployment assumption) the per-query
/// error is `ρ = 2^{-LOG_BLOWUP}`, so `q` queries give error `2^{-LOG_BLOWUP·q}` and
/// reaching `SECURITY_BITS` needs `q = ⌈SECURITY_BITS / LOG_BLOWUP⌉`. (The *provable*
/// unique-decoding bound `δ = (1−ρ)/2` is weaker — error `((1+ρ)/2)^q` — and needs ~3×
/// more queries for the same target; raise this if a provable guarantee is required.
/// Adding proof-of-work grinding of `g` bits would let `q` drop by `g/LOG_BLOWUP`.)
/// Proof size scales ~linearly in `NUM_QUERIES`.
const NUM_QUERIES: usize = SECURITY_BITS.div_ceil(LOG_BLOWUP);

/// A BaseFold commitment: the Merkle root of the base codeword and the variable count.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment {
    /// Root of the base-layer (`C_0`) Merkle tree.
    pub root: [u8; 32],
    /// Number of variables `ν` of the committed multilinear.
    pub num_vars: usize,
}

/// Prover-retained data: every codeword layer needed to answer queries plus its tree.
#[derive(Clone, Debug)]
pub struct ProverData {
    num_vars: usize,
    evals: Vec<Fp>,
    codeword0: Vec<Ext2>,
    tree0: MerkleTree,
}

/// One layer's opened pair `(C_i[p], C_i[p + half_i])` with authentication paths.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LayerOpening {
    /// The codeword value at the low index `p`.
    pub lo: Ext2,
    /// The codeword value at the high index `p + half_i`.
    pub hi: Ext2,
    /// Authentication path for `lo`.
    pub lo_path: MerklePath,
    /// Authentication path for `hi`.
    pub hi_path: MerklePath,
}

/// One FRI query: the base position and the opened pair at each folding layer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryProof {
    /// The sampled base-layer position in `[0, n0/2)`.
    pub pos: usize,
    /// Opened pairs, one per folding layer `C_0..C_{ν−1}`.
    pub layers: Vec<LayerOpening>,
}

/// A BaseFold evaluation proof for `f(r) = value`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProof {
    /// The claimed value `f(r)`.
    pub value: Ext2,
    /// Sum-check round polynomials, each as evaluations at `0, 1, 2` (degree 2).
    pub round_polys: Vec<[Ext2; 3]>,
    /// Merkle roots of the intermediate folded codewords `C_1, …, C_{ν−1}`.
    pub layer_roots: Vec<[u8; 32]>,
    /// The fully-folded (degree-0) codeword constant.
    pub final_constant: Ext2,
    /// FRI consistency queries.
    pub queries: Vec<QueryProof>,
}

/// Lagrange-interpolate degree-2 evals at `0,1,2` and evaluate at `r`:
/// `g(r) = g0·(r−1)(r−2)/2 − g1·r(r−2) + g2·r(r−1)/2`.
fn lagrange3(g: &[Ext2; 3], r: Ext2) -> Ext2 {
    let half = Ext2::from_base(Fp::new(2).inv().unwrap());
    let r1 = r - Ext2::ONE;
    let r2 = r - Ext2::from_base(Fp::new(2));
    g[0] * (r1 * r2) * half - g[1] * (r * r2) + g[2] * (r * r1) * half
}

/// Absorb the commitment context (root, ν, point) so the proof is bound to the claim.
fn absorb_context(tr: &mut Transcript, comm: &Commitment, point: &[Ext2]) {
    tr.absorb_bytes(b"basefold/root", &comm.root);
    tr.absorb_bytes(b"basefold/nv", &(comm.num_vars as u64).to_le_bytes());
    for &p in point {
        tr.absorb_ext(b"basefold/pt", p);
    }
}

/// Squeeze `NUM_QUERIES` base-layer positions in `[0, half0)`.
fn sample_queries(tr: &mut Transcript, half0: usize) -> Vec<usize> {
    (0..NUM_QUERIES)
        .map(|_| (tr.challenge_fp(b"basefold/query").to_u64() as usize) % half0)
        .collect()
}

/// Commit to a multilinear given by its evaluation table over `{0,1}^ν` (`n = 2^ν ≥ 2`).
pub fn commit(evals: &[Fp]) -> (Commitment, ProverData) {
    let n = evals.len();
    assert!(n.is_power_of_two() && n >= 2, "need ≥ 2 evaluations");
    let num_vars = code::log2(n);
    let coeffs = code::evals_to_coeffs(evals);
    let codeword0 = code::lift(&code::encode(&coeffs, LOG_BLOWUP));
    let tree0 = MerkleTree::commit(&codeword0);
    let root = tree0.root();
    (
        Commitment { root, num_vars },
        ProverData {
            num_vars,
            evals: evals.to_vec(),
            codeword0,
            tree0,
        },
    )
}

/// Prove `f(point) = value` for the committed `f`.
pub fn open(
    tr: &mut Transcript,
    comm: &Commitment,
    data: &ProverData,
    point: &[Ext2],
) -> OpenProof {
    assert_eq!(point.len(), data.num_vars);
    absorb_context(tr, comm, point);

    let nu = data.num_vars;
    let f0 = code::lift(&data.evals);
    let value = mle::eval(&f0, point);

    let mut f_tab = f0;
    let mut eq_tab = mle::eq_table(point);

    // Codeword layers C_0..C_{ν−1} kept for the query phase; trees aligned. `codewords`
    // grows by folding its own last layer, so no separate accumulator is needed.
    let mut codewords = vec![data.codeword0.clone()];
    let mut trees = vec![data.tree0.clone()];

    let mut round_polys = Vec::with_capacity(nu);
    let mut layer_roots = Vec::new();
    let two = Ext2::from_base(Fp::new(2));

    for round in 0..nu {
        // Degree-2 round polynomial g(c) = Σ_j eq_lin(c)·f_lin(c), c ∈ {0,1,2}.
        let mut g = [Ext2::ZERO; 3];
        for j in 0..f_tab.len() / 2 {
            let (f0, f1) = (f_tab[2 * j], f_tab[2 * j + 1]);
            let (e0, e1) = (eq_tab[2 * j], eq_tab[2 * j + 1]);
            g[0] += e0 * f0;
            g[1] += e1 * f1;
            let f2 = f0 + two * (f1 - f0);
            let e2 = e0 + two * (e1 - e0);
            g[2] += e2 * f2;
        }
        for &e in &g {
            tr.absorb_ext(b"basefold/sc", e);
        }
        let alpha = tr.challenge_ext(b"basefold/alpha");
        round_polys.push(g);

        f_tab = fold_evals(&f_tab, alpha);
        eq_tab = fold_evals(&eq_tab, alpha);
        let folded = fold_codeword(codewords.last().unwrap(), alpha);

        if round < nu - 1 {
            // C_{round+1} is an intermediate layer: commit and bind its root.
            let tree = MerkleTree::commit(&folded);
            tr.absorb_bytes(b"basefold/layer", &tree.root());
            layer_roots.push(tree.root());
            trees.push(tree);
        }
        codewords.push(folded);
    }
    debug_assert_eq!(codewords.last().unwrap().len(), FINAL_LEN);

    let final_constant = codewords.last().unwrap()[0];
    tr.absorb_ext(b"basefold/final", final_constant);

    let half0 = data.codeword0.len() / 2;
    let queries = sample_queries(tr, half0)
        .into_iter()
        .map(|pos| {
            let layers = (0..nu)
                .map(|i| {
                    let half_i = (data.codeword0.len() >> i) / 2;
                    let p = pos % half_i;
                    LayerOpening {
                        lo: codewords[i][p],
                        hi: codewords[i][p + half_i],
                        lo_path: trees[i].open(p),
                        hi_path: trees[i].open(p + half_i),
                    }
                })
                .collect();
            QueryProof { pos, layers }
        })
        .collect();

    OpenProof {
        value,
        round_polys,
        layer_roots,
        final_constant,
        queries,
    }
}

/// FRI-fold a single point: `(lo + hi)/2 + α·(lo − hi)/(2·g^p)`, the per-position image
/// of [`fold_codeword`](crate::code::fold_codeword). `two_inv` (= 2⁻¹ in `K`) is passed
/// in so the constant inversion is hoisted out of the query loop.
fn fold_point(lo: Ext2, hi: Ext2, alpha: Ext2, g_pow_p: Fp, two_inv: Ext2) -> Ext2 {
    let s_inv = Ext2::from_base(g_pow_p.inv().expect("domain points are nonzero"));
    (lo + hi) * two_inv + alpha * ((lo - hi) * two_inv * s_inv)
}

/// Verify a BaseFold evaluation proof against `comm` and the claim `f(point) = value`.
pub fn verify(
    tr: &mut Transcript,
    comm: &Commitment,
    point: &[Ext2],
    proof: &OpenProof,
) -> Result<(), SnarkError> {
    let nu = comm.num_vars;
    if point.len() != nu {
        return Err(SnarkError::Malformed("point arity mismatch".into()));
    }
    if proof.round_polys.len() != nu || proof.layer_roots.len() != nu.saturating_sub(1) {
        return Err(SnarkError::Verify("proof shape mismatch".into()));
    }
    absorb_context(tr, comm, point);

    // Replay the sum-check, deriving folding challenges and interleaving layer roots
    // in exactly the order `open` bound them.
    let mut claim = proof.value;
    let mut alphas = Vec::with_capacity(nu);
    for (round, g) in proof.round_polys.iter().enumerate() {
        if g[0] + g[1] != claim {
            return Err(SnarkError::Verify(format!(
                "sum-check round {round} mismatch"
            )));
        }
        for &e in g {
            tr.absorb_ext(b"basefold/sc", e);
        }
        let alpha = tr.challenge_ext(b"basefold/alpha");
        claim = lagrange3(g, alpha);
        alphas.push(alpha);
        if round < nu - 1 {
            tr.absorb_bytes(b"basefold/layer", &proof.layer_roots[round]);
        }
    }
    tr.absorb_ext(b"basefold/final", proof.final_constant);

    // The sum-check reduces Σ_x eq(r,x)f(x) to eq(r,α)·f(α); f(α) is the folded constant.
    if claim != eq_eval(point, &alphas) * proof.final_constant {
        return Err(SnarkError::Verify("final folding/eval tie failed".into()));
    }

    // FRI query phase. Hoist the per-layer domain generators and 2⁻¹ (loop-invariant
    // across all queries) out of the inner loops.
    let n0 = 1usize << (nu + LOG_BLOWUP);
    let half0 = n0 / 2;
    let layer_gen: Vec<Fp> = (0..nu)
        .map(|i| root_of_unity(nu + LOG_BLOWUP - i))
        .collect();
    let two_inv = Ext2::from_base(Fp::new(2).inv().expect("2 ≠ 0 in F_q"));
    let positions = sample_queries(tr, half0);
    if proof.queries.len() != positions.len() {
        return Err(SnarkError::Verify("query count mismatch".into()));
    }
    for (q, &pos) in proof.queries.iter().zip(positions.iter()) {
        if q.pos != pos || q.layers.len() != nu {
            return Err(SnarkError::Verify("query position mismatch".into()));
        }
        for i in 0..nu {
            let n_i = n0 >> i;
            let half_i = n_i / 2;
            let p = pos % half_i;
            let root_i = if i == 0 {
                &comm.root
            } else {
                &proof.layer_roots[i - 1]
            };
            let lyr = &q.layers[i];
            if !verify_path(root_i, n_i, p, lyr.lo, &lyr.lo_path)
                || !verify_path(root_i, n_i, p + half_i, lyr.hi, &lyr.hi_path)
            {
                return Err(SnarkError::PcsOpening(format!("Merkle path layer {i}")));
            }
            let g_pow_p = layer_gen[i].pow(p as u64);
            let folded = fold_point(lyr.lo, lyr.hi, alphas[i], g_pow_p, two_inv);
            let expected = if i + 1 < nu {
                // C_{i+1}[p]: lo or hi of the next layer's opened pair.
                let half_next = half_i / 2;
                if p < half_next {
                    q.layers[i + 1].lo
                } else {
                    q.layers[i + 1].hi
                }
            } else {
                proof.final_constant
            };
            if folded != expected {
                return Err(SnarkError::PcsOpening(format!(
                    "fold consistency layer {i}"
                )));
            }
        }
    }
    Ok(())
}

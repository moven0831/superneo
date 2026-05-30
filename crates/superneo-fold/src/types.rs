//! Relation types for the folding scheme (Definitions 11–14).
//!
//! Conventions used throughout (see the M4 design notes in the plan):
//!   * `M_1 = I_m` (Remark 3), so `m = n_f` and `~(M_1 z) = ~z` is directly available.
//!   * `n_f = d · n_R`; the witness field vector `z ∈ F^{n_f}` corresponds to the ring
//!     witness `ẑ ∈ R_F^{n_R}` via the coefficient embedding.
//!   * the row dimension `m` is zero-padded to `cap_m = 2^{log_m}` for the sum-check.
//!   * no public input in the fold relations (`n_in = 0`); public IO is an IVC-layer
//!     concern, handled separately.

use superneo_commit::{Commitment, Params, PublicParams};
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_ring::maps::{bar_matrix_row, embed_vector, mbar_row_dot_z};
use superneo_ring::{RingElem, RingK, D};

use crate::multilinear::{log2_pow2, next_pow2};

/// A CCS constraint polynomial `f ∈ F^{<u}[X_1, …, X_t]`, as a sum of monomials.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparsePoly {
    /// Number of variables `t`.
    pub t: usize,
    /// Terms `(coefficient, exponent-vector of length t)`.
    pub terms: Vec<(Fp, Vec<u32>)>,
}

impl SparsePoly {
    /// The zero polynomial in `t` variables.
    pub fn zero(t: usize) -> Self {
        SparsePoly {
            t,
            terms: Vec::new(),
        }
    }

    /// Construct from terms, checking each exponent vector has length `t`.
    pub fn new(t: usize, terms: Vec<(Fp, Vec<u32>)>) -> Self {
        for (_, e) in &terms {
            assert_eq!(e.len(), t, "exponent vector must have length t");
        }
        SparsePoly { t, terms }
    }

    /// Total degree `max_term Σ exponents` (0 for the zero polynomial).
    pub fn degree(&self) -> usize {
        self.terms
            .iter()
            .map(|(_, e)| e.iter().map(|&x| x as usize).sum())
            .max()
            .unwrap_or(0)
    }

    /// Evaluate at extension-field arguments.
    pub fn eval_ext(&self, args: &[Ext2]) -> Ext2 {
        debug_assert_eq!(args.len(), self.t);
        let mut acc = Ext2::ZERO;
        for (coeff, exps) in &self.terms {
            let mut term = Ext2::from_base(*coeff);
            for (v, &e) in exps.iter().enumerate() {
                for _ in 0..e {
                    term *= args[v];
                }
            }
            acc += term;
        }
        acc
    }
}

/// A CCS structure `s = ({M_j}_{j∈[t]}, f)` (Definition 11), with `M_1 = I_m`.
#[derive(Clone, Debug)]
pub struct CcsStructure {
    /// The `t` matrices, each dense `m × n_f` (= `m × m`); `matrices[0] = I_m`.
    pub matrices: Vec<Vec<Vec<Fp>>>,
    /// The constraint polynomial.
    pub f: SparsePoly,
    /// Number of rows / constraints (`= n_f`).
    pub m: usize,
}

impl CcsStructure {
    /// Number of matrices `t`.
    pub fn t(&self) -> usize {
        self.matrices.len()
    }

    /// Build and validate: square matrices, `M_1 = I_m`, matching `f.t`.
    pub fn new(matrices: Vec<Vec<Vec<Fp>>>, f: SparsePoly) -> Self {
        assert!(!matrices.is_empty(), "need at least M_1");
        let m = matrices[0].len();
        assert_eq!(f.t, matrices.len(), "f arity must equal t");
        for mat in &matrices {
            assert_eq!(mat.len(), m, "all matrices share row count m");
            for row in mat {
                assert_eq!(row.len(), m, "matrices are square (m = n_f)");
            }
        }
        // M_1 = I_m.
        for (i, row) in matrices[0].iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                let want = if i == j { Fp::ONE } else { Fp::ZERO };
                assert_eq!(v, want, "M_1 must be the identity (Remark 3)");
            }
        }
        CcsStructure { matrices, f, m }
    }

    /// The identity-only structure (`t = 1`, `f = 0`): satisfied by every witness,
    /// useful for exercising the norm / evaluation / RLC / DEC machinery in isolation.
    pub fn identity(m: usize) -> Self {
        let mut id = vec![vec![Fp::ZERO; m]; m];
        for (i, row) in id.iter_mut().enumerate() {
            row[i] = Fp::ONE;
        }
        CcsStructure::new(vec![id], SparsePoly::zero(1))
    }
}

/// Global reduction parameters (Definition 14).
#[derive(Clone)]
pub struct GlobalParams {
    /// Ring degree `d = 54`.
    pub d: usize,
    /// Ring-vector length `n_R`.
    pub n_r: usize,
    /// Witness field length `n_f = d · n_R = m`.
    pub n_f: usize,
    /// Row count `m`.
    pub m: usize,
    /// Padded hypercube dimension `log_m` (so the cube has `2^{log_m} ≥ m` points).
    pub log_m: usize,
    /// Decomposition base `b`.
    pub b: u64,
    /// Decomposition depth / running-instance count `k` (`B = b^k`).
    pub k: usize,
    /// High norm bound `B = b^k`.
    pub cap_b: u64,
    /// ℓ∞ bound on strong-sampling-set coefficients (Definition 17).
    pub chal_bound: i64,
}

impl GlobalParams {
    /// The padded hypercube size `2^{log_m}`.
    pub fn cap_m(&self) -> usize {
        1 << self.log_m
    }

    /// The cryptographic parameter set backing these global parameters.
    pub fn params(&self) -> Params {
        Params::GOLDILOCKS_B2
    }

    /// Goldilocks parameters (Appendix B.2) for a ring-vector length `n_R`.
    pub fn goldilocks(n_r: usize) -> Self {
        let p = Params::GOLDILOCKS_B2;
        // Enforce the folding-admissibility guard (Definition 14) on the production path,
        // so an invalid parameter set fails fast rather than silently breaking binding.
        p.validate()
            .expect("Goldilocks parameter set must satisfy the folding-admissibility guard");
        let n_f = (p.d) * n_r;
        let cap_m = next_pow2(n_f);
        GlobalParams {
            d: p.d,
            n_r,
            n_f,
            m: n_f,
            log_m: log2_pow2(cap_m),
            b: p.b,
            k: p.k,
            cap_b: p.big_b,
            chal_bound: p.challenge_coeff_bound(),
        }
    }
}

// ---- Instances & witnesses ----------------------------------------------------

/// Norm-bounded CCS instance (public part), `n_in = 0` (Definition 12).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CcsInstance {
    /// Commitment `c = L(ẑ)`.
    pub c: Commitment,
}

/// CCS witness (Definition 12): the field vector `z ∈ F^{n_f}`.
#[derive(Clone, Debug)]
pub struct CcsWitness {
    /// Field witness `z`, length `n_f`.
    pub z: Vec<Fp>,
}

/// Norm-bounded CCS evaluation instance (Definition 13).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CeInstance {
    /// Commitment `c = L(ẑ)`.
    pub c: Commitment,
    /// Evaluation point `r ∈ K^{log_m}`.
    pub r: Vec<Ext2>,
    /// Evaluation claims `y_j = ~(M̄_j ẑ)(r) ∈ R_K`, one per matrix.
    pub y: Vec<RingK>,
}

/// CE witness (Definition 13): the ring vector `ẑ ∈ R_F^{n_R}`.
#[derive(Clone, Debug)]
pub struct CeWitness {
    /// Ring witness `ẑ`, length `n_R`.
    pub z_ring: Vec<RingElem>,
}

// ---- Relation helpers ---------------------------------------------------------

/// The SuperNeo embedding `ẑ = embed(z) ∈ R_F^{n_R}` of a field witness.
pub fn embed_witness(z: &[Fp]) -> Vec<RingElem> {
    embed_vector(z).expect("witness length must be a multiple of d")
}

/// Flatten a ring witness `ẑ` to its field coefficient vector `z ∈ F^{d·n_R}`.
pub fn flatten_ring(z_ring: &[RingElem]) -> Vec<Fp> {
    let mut out = Vec::with_capacity(z_ring.len() * D);
    for e in z_ring {
        out.extend_from_slice(e.coeffs());
    }
    out
}

/// The field matrix-vector product `M · z ∈ F^m`.
pub fn field_matvec(mat: &[Vec<Fp>], z: &[Fp]) -> Vec<Fp> {
    mat.iter()
        .map(|row| {
            let mut acc = Fp::ZERO;
            for (a, b) in row.iter().zip(z.iter()) {
                acc += *a * *b;
            }
            acc
        })
        .collect()
}

/// The bar-lifted ring matrix-vector product `M̄ · ẑ ∈ R_F^m`.
///
/// Each row of `mat` (length `n_f`) is bar-lifted to `n_R` ring elements (Definition 8),
/// then dotted with the ring witness (Theorem 6: `ct` of this recovers `M z`).
pub fn bar_ring_matvec(mat: &[Vec<Fp>], z_ring: &[RingElem]) -> Vec<RingElem> {
    mat.iter()
        .map(|row| {
            let bar_row = bar_matrix_row(row).expect("row length must be a multiple of d");
            mbar_row_dot_z(&bar_row, z_ring)
        })
        .collect()
}

/// Evaluation claims `y_j = ~(M̄_j ẑ)(r)` for all matrices (the CE `y` vector).
pub fn compute_evals(s: &CcsStructure, z_ring: &[RingElem], r: &[Ext2]) -> Vec<RingK> {
    use crate::multilinear::mle_eval_ring;
    s.matrices
        .iter()
        .map(|mat| mle_eval_ring(&bar_ring_matvec(mat, z_ring), r))
        .collect()
}

// ---- Proof types --------------------------------------------------------------

/// Π_CCS proof: the sum-check transcript plus the claimed evaluations `y'_{i,j}`.
#[derive(Clone, Debug)]
pub struct PiCcsProof {
    /// Sum-check proof.
    pub sumcheck: crate::sumcheck::SumCheckProof,
    /// `evals[i][j] = y'_{i,j} ∈ R_K`, for `i ∈ [K+k]`, `j ∈ [t]`.
    pub evals: Vec<Vec<RingK>>,
}

/// Π_RLC proof: the strong-sampling-set challenges `ρ_i ∈ R_F`.
#[derive(Clone, Debug)]
pub struct PiRlcProof {
    /// `K + k` challenges.
    pub rhos: Vec<RingElem>,
}

/// Π_DEC proof: child commitments and their evaluation claims.
#[derive(Clone, Debug)]
pub struct PiDecProof {
    /// `k` child commitments `c_i = L(ẑ_i)`.
    pub child_commitments: Vec<Commitment>,
    /// `child_evals[i][j] = y_{i,j} ∈ R_K`, for `i ∈ [k]`, `j ∈ [t]`.
    pub child_evals: Vec<Vec<RingK>>,
}

/// A full fold proof: `Π_DEC ∘ Π_RLC ∘ Π_CCS` (Theorem 3).
#[derive(Clone, Debug)]
pub struct FoldProof {
    /// Π_CCS sub-proof.
    pub ccs: PiCcsProof,
    /// Π_RLC sub-proof.
    pub rlc: PiRlcProof,
    /// Π_DEC sub-proof.
    pub dec: PiDecProof,
}

/// Commit to a field witness via its SuperNeo embedding: `L(z) = Commit(embed(z))`.
pub fn commit_witness(pp: &PublicParams, z: &[Fp]) -> Commitment {
    pp.commit(&embed_witness(z))
        .expect("embedded witness has n_R ring elements")
}

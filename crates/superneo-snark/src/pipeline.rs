//! End-to-end integration (M8): fold a multi-step IVC run and compress the result.
//!
//! This module ties the whole stack together — the fold (`superneo-fold`), the native
//! IVC loop (`superneo-ivc`), and the final compression (this crate) — behind one
//! [`run_pipeline`] entry point that returns a [`PipelineReport`] (timings, the IVC
//! digest, an estimated proof size, and the verification verdict). It also constructs
//! the plan's demonstration circuits (a counter step, a Fibonacci transition, and a
//! multiplication gate) through the recursive-verifier CCS builder, showing that
//! ordinary computations compile to the SuperNeo `t = 4` CCS shape.
//!
//! Foldable vs. demonstration: the folding scheme requires *low-norm* witnesses
//! (`‖z‖∞ < b`), which the IVC pipeline uses. The counter/Fibonacci circuits are shown
//! as satisfied CCS instances (their witnesses are arbitrary field elements); folding
//! such a circuit through the norm-bounded scheme needs the augmented-witness
//! decomposition noted as the M7 residual.

use std::time::{Duration, Instant};

use superneo_commit::PublicParams;
use superneo_field::fp::Fp;
use superneo_fold::{CcsStructure, CcsWitness, GlobalParams, SparsePoly};
use superneo_ivc::circuit::{ConstraintSystem, Lc};
use superneo_ivc::{prove_ivc_with_witnesses, verify_ivc};

use crate::basefold::OpenProof;
use crate::error::SnarkError;
use crate::spartan::{compress, verify, SnarkProof};

/// A summary of one end-to-end pipeline run.
#[derive(Clone, Debug)]
pub struct PipelineReport {
    /// Number of IVC steps folded.
    pub steps: usize,
    /// Accumulator width (carried CE instances).
    pub k: usize,
    /// Witness field length per instance (`n_f = d·n_R`).
    pub n_f: usize,
    /// Time to prove the IVC run (all folds).
    pub fold_time: Duration,
    /// Time to compress the final accumulator.
    pub compress_time: Duration,
    /// Time to verify the compressed proof.
    pub verify_time: Duration,
    /// The final Construction-2 IVC digest.
    pub ivc_digest: [u8; 32],
    /// Estimated serialized size of the compressed proof, in bytes.
    pub proof_bytes: usize,
    /// Whether the compressed proof verified.
    pub verified: bool,
}

/// Run the full pipeline: prove the IVC run, verify it, compress the final accumulator,
/// and verify the compressed proof. Returns a [`PipelineReport`].
pub fn run_pipeline(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    steps: &[Vec<CcsWitness>],
) -> Result<PipelineReport, SnarkError> {
    let t0 = Instant::now();
    let (ivc_proof, final_wit) = prove_ivc_with_witnesses(gp, pp, s, steps)
        .map_err(|e| SnarkError::Malformed(format!("IVC prove: {e}")))?;
    let fold_time = t0.elapsed();

    verify_ivc(gp, pp, s, &ivc_proof)
        .map_err(|e| SnarkError::Verify(format!("IVC verify: {e}")))?;

    let t1 = Instant::now();
    let snark = compress(gp, pp, s, &ivc_proof.final_acc, &final_wit)?;
    let compress_time = t1.elapsed();

    let t2 = Instant::now();
    let verified = verify(gp, pp, s, &ivc_proof.final_acc, &snark).is_ok();
    let verify_time = t2.elapsed();

    Ok(PipelineReport {
        steps: steps.len(),
        k: gp.k,
        n_f: gp.n_f,
        fold_time,
        compress_time,
        verify_time,
        ivc_digest: ivc_proof.final_digest,
        proof_bytes: proof_size_bytes(&snark),
        verified,
    })
}

// ---- proof-size accounting ---------------------------------------------------

/// Bytes for a `K` element (two field elements) and a hash digest.
const K_BYTES: usize = 16;
const HASH_BYTES: usize = 32;

fn open_proof_bytes(p: &OpenProof) -> usize {
    let mut n = K_BYTES; // value
    n += p.round_polys.len() * 3 * K_BYTES;
    n += p.layer_roots.len() * HASH_BYTES;
    n += K_BYTES; // final_constant
    for q in &p.queries {
        n += 8; // pos
        for lyr in &q.layers {
            n += 2 * K_BYTES; // lo, hi
            n += (lyr.lo_path.siblings.len() + lyr.hi_path.siblings.len()) * HASH_BYTES;
        }
    }
    n
}

/// An estimate of the compressed proof's serialized size, in bytes.
pub fn proof_size_bytes(p: &SnarkProof) -> usize {
    let mut n = HASH_BYTES + 8; // z_commit (root + num_vars)
    n += p
        .lin_rounds
        .iter()
        .map(|r| r.len() * K_BYTES)
        .sum::<usize>();
    n += p
        .norm_rounds
        .iter()
        .map(|r| r.len() * K_BYTES)
        .sum::<usize>();
    n += open_proof_bytes(&p.lin_open);
    n += open_proof_bytes(&p.norm_open);
    n
}

// ---- foldable workloads ------------------------------------------------------

/// A single multiplication gate `z[0]·z[1] = z[2]` as a SuperNeo CCS structure
/// (`t = 4`, `f = X_2·X_3 − X_4`), padded to `m`.
pub fn multiplication_structure(m: usize) -> CcsStructure {
    let mut id = vec![vec![Fp::ZERO; m]; m];
    let mut a = vec![vec![Fp::ZERO; m]; m];
    let mut b = vec![vec![Fp::ZERO; m]; m];
    let mut c = vec![vec![Fp::ZERO; m]; m];
    for (i, row) in id.iter_mut().enumerate() {
        row[i] = Fp::ONE;
    }
    a[0][0] = Fp::ONE; // A·z = z[0]
    b[0][1] = Fp::ONE; // B·z = z[1]
    c[0][2] = Fp::ONE; // C·z = z[2]
    let f = SparsePoly::new(
        4,
        vec![(Fp::ONE, vec![0, 1, 1, 0]), (-Fp::ONE, vec![0, 0, 0, 1])],
    );
    CcsStructure::new(vec![id, a, b, c], f)
}

/// A low-norm multiplication-gate witness `b·b = b²` for `b ∈ {0,1}` (so `‖z‖∞ < 2`):
/// `(1,1,1)` or `(0,0,0)`. Used to fold a real constraint end-to-end.
pub fn mul_gate_witness(m: usize, bit: u64) -> CcsWitness {
    let mut z = vec![Fp::ZERO; m];
    z[0] = Fp::new(bit);
    z[1] = Fp::new(bit);
    z[2] = Fp::new(bit); // bit·bit = bit for bit ∈ {0,1}
    CcsWitness { z }
}

// ---- demonstration circuits (CCS encoding via the recursive-verifier builder) ----

/// A satisfied CCS instance for one counter step `out = in + 1`.
pub fn counter_circuit(input: u64) -> ConstraintSystem {
    let mut cs = ConstraintSystem::new();
    let in_var = cs.alloc(Fp::new(input));
    // out = in + 1, computed in the field so the witness matches the field-arithmetic
    // constraint for any input (a raw `input + 1` u64 add would overflow / wrap).
    let out_var = cs.alloc(Fp::new(input) + Fp::ONE);
    // (in + 1)·1 = out.
    let in_plus_one = Lc::from_var(in_var).add(&cs.constant(Fp::ONE));
    cs.enforce(in_plus_one, Lc::from_var(cs.one()), Lc::from_var(out_var));
    cs
}

/// A satisfied CCS instance for one Fibonacci transition `(a, b) → (b, a + b)`.
pub fn fibonacci_circuit(a: u64, b: u64) -> ConstraintSystem {
    let mut cs = ConstraintSystem::new();
    let a_var = cs.alloc(Fp::new(a));
    let b_var = cs.alloc(Fp::new(b));
    let next_a = cs.alloc(Fp::new(b));
    // a + b in the field (avoids a u64 overflow in the witness value).
    let next_b = cs.alloc(Fp::new(a) + Fp::new(b));
    // next_a = b.
    cs.assert_eq(&Lc::from_var(next_a), &Lc::from_var(b_var));
    // next_b = a + b.
    cs.assert_eq(
        &Lc::from_var(next_b),
        &Lc::from_var(a_var).add(&Lc::from_var(b_var)),
    );
    cs
}

/// A satisfied CCS instance for a multiplication gate `x·y = z` (exercises the quadratic
/// `f`-path: the builder emits a genuine rank-1 constraint).
pub fn multiplier_circuit(x: u64, y: u64) -> ConstraintSystem {
    let mut cs = ConstraintSystem::new();
    let x_var = cs.alloc(Fp::new(x));
    let y_var = cs.alloc(Fp::new(y));
    let _z = cs.mul(&Lc::from_var(x_var), &Lc::from_var(y_var));
    cs
}

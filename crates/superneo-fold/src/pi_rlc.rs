//! Π_RLC: the weak interactive reduction `CE(b)^{K+k} → CE(B)` (§7.4, Lemma 4).
//!
//! Sample `ρ_1, …, ρ_{K+k} ← C` and take the ring linear combination of the input
//! commitments, evaluation claims, and witnesses. By the commitment and evaluation
//! homomorphisms (Theorem 7) the result is a single CE instance of norm bound
//! `B = b^k` (the norm grows by at most `(K+k)·T·(b−1) < B`, Definition 14).

use superneo_commit::Commitment;
use superneo_ring::{RingElem, RingK};

use crate::error::FoldError;
use crate::transcript::Transcript;
use crate::types::{CeInstance, CeWitness, GlobalParams, PiRlcProof};

/// Absorb the input CE instances (binding) and squeeze the `ρ` challenges.
fn absorb_and_sample(tr: &mut Transcript, gp: &GlobalParams, ce: &[CeInstance]) -> Vec<RingElem> {
    for inst in ce {
        tr.absorb_commitment(b"pirlc/c", &inst.c);
        for &ri in &inst.r {
            tr.absorb_ext(b"pirlc/r", ri);
        }
        for yj in &inst.y {
            for c in yj.coeffs() {
                tr.absorb_ext(b"pirlc/y", *c);
            }
        }
    }
    tr.challenge_ring_set(b"pirlc/rho", ce.len(), gp.chal_bound)
}

/// Combine the public parts (commitments, eval claims) under challenges `ρ`.
fn combine_public(ce: &[CeInstance], rhos: &[RingElem]) -> (Commitment, Vec<RingK>) {
    let kappa = ce[0].c.dim();
    let t = ce[0].y.len();
    let mut c_out = Commitment::zero(kappa);
    let mut y_out = vec![RingK::ZERO; t];
    for (inst, rho) in ce.iter().zip(rhos.iter()) {
        c_out = c_out.add(&inst.c.scale_ring(rho));
        for (j, yj) in inst.y.iter().enumerate() {
            y_out[j] = y_out[j].add(&yj.mul_rf(rho));
        }
    }
    (c_out, y_out)
}

/// Π_RLC prover.
pub fn pi_rlc_prove(
    tr: &mut Transcript,
    gp: &GlobalParams,
    ce: &[CeInstance],
    ce_wit: &[CeWitness],
) -> (CeInstance, CeWitness, PiRlcProof) {
    let rhos = absorb_and_sample(tr, gp, ce);
    let (c_out, y_out) = combine_public(ce, &rhos);

    // z_out = Σ ρ_i · ẑ_i  (ring-vector linear combination).
    let n_r = ce_wit[0].z_ring.len();
    let mut z_out = vec![RingElem::ZERO; n_r];
    for (w, rho) in ce_wit.iter().zip(rhos.iter()) {
        for (col, zc) in w.z_ring.iter().enumerate() {
            z_out[col] = z_out[col] + (*zc * *rho);
        }
    }

    let inst = CeInstance {
        c: c_out,
        r: ce[0].r.clone(),
        y: y_out,
    };
    (inst, CeWitness { z_ring: z_out }, PiRlcProof { rhos })
}

/// Π_RLC verifier: recomputes the combination from the public instances.
pub fn pi_rlc_verify(
    tr: &mut Transcript,
    gp: &GlobalParams,
    ce: &[CeInstance],
    proof: &PiRlcProof,
) -> Result<CeInstance, FoldError> {
    let rhos = absorb_and_sample(tr, gp, ce);
    if rhos != proof.rhos {
        return Err(FoldError::Reduction("Π_RLC challenge mismatch".into()));
    }
    let (c_out, y_out) = combine_public(ce, &rhos);
    Ok(CeInstance {
        c: c_out,
        r: ce[0].r.clone(),
        y: y_out,
    })
}

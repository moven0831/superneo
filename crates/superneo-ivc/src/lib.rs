//! `superneo-ivc` — incrementally verifiable computation on the folding scheme.
//!
//! **Phase 1 (native, this module).** An IVC driver that folds a sequence of steps
//! into a running accumulator of `k` CE instances, binding each iteration to the next
//! with a Construction-2 IO digest `digest_i = H(digest_{i-1}, i, acc_i)`. The genesis
//! accumulator is the trivial all-zero CE tuple (`c = L(0)`, `y = 0`), which lies in
//! `CE(b, L)`. Verification replays the per-step `verify_fold` and the digest chain.
//!
//! Each step contributes one or more fresh CCS instances (the step's computation); the
//! single shared transcript across steps mirrors the Fiat–Shamir IVC pattern. The CCS
//! constraints enforce per-step correctness; the digest binds the public transcript of
//! the whole run.
//!
//! **Phase 2 (recursive verifier circuit, [`circuit`]).** The folding verifier's work
//! expressed as CCS constraints — closing a true recursive IVC loop. This PoC implements
//! the dominant component (the sum-check verifier) plus the gadget framework; see the
//! module for the documented residuals.

pub mod circuit;
pub mod error;

pub use error::IvcError;

use superneo_commit::{Commitment, PublicParams};
use superneo_field::ext2::Ext2;
use superneo_ring::{RingElem, RingK};

use superneo_fold::types::commit_witness;
use superneo_fold::{
    fold, verify_fold, CcsInstance, CcsStructure, CcsWitness, CeInstance, CeWitness, FoldProof,
    GlobalParams, Transcript,
};

/// Domain-separation label for the IVC transcript.
const IVC_DOMAIN: &[u8] = b"superneo/ivc/v1";
/// The genesis IO digest.
const GENESIS_DIGEST: [u8; 32] = [0u8; 32];

/// A full native-IVC proof for a run of `steps` folding steps.
#[derive(Clone, Debug)]
pub struct IvcProof {
    /// Fresh CCS instances contributed at each step.
    pub fresh: Vec<Vec<CcsInstance>>,
    /// The fold proof produced at each step.
    pub fold_proofs: Vec<FoldProof>,
    /// The final accumulator (`k` CE instances).
    pub final_acc: Vec<CeInstance>,
    /// The final Construction-2 IO digest.
    pub final_digest: [u8; 32],
}

/// The genesis accumulator instances: trivial CE tuples `(c = 0, r = 0, y = 0)`.
fn genesis_instances(gp: &GlobalParams, pp: &PublicParams, t: usize) -> Vec<CeInstance> {
    let r0 = vec![Ext2::ZERO; gp.log_m];
    let y0 = vec![RingK::ZERO; t];
    (0..gp.k)
        .map(|_| CeInstance {
            c: Commitment::zero(pp.kappa),
            r: r0.clone(),
            y: y0.clone(),
        })
        .collect()
}

/// The genesis accumulator witnesses: all-zero ring vectors.
fn genesis_witnesses(gp: &GlobalParams) -> Vec<CeWitness> {
    (0..gp.k)
        .map(|_| CeWitness {
            z_ring: vec![RingElem::ZERO; gp.n_r],
        })
        .collect()
}

/// Construction-2 IO digest update `H(prev, step, acc)`.
fn fold_digest(prev: &[u8; 32], step: usize, acc: &[CeInstance]) -> [u8; 32] {
    let mut h = blake3::Hasher::new();
    h.update(b"superneo/ivc/digest");
    h.update(prev);
    h.update(&(step as u64).to_le_bytes());
    for inst in acc {
        for ring in &inst.c.0 {
            for c in ring.coeffs() {
                h.update(&c.to_u64().to_le_bytes());
            }
        }
    }
    *h.finalize().as_bytes()
}

/// Prove a native IVC run. Each entry of `steps` lists the fresh CCS witnesses folded
/// in at that step (commitments are recomputed internally).
pub fn prove_ivc(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    steps: &[Vec<CcsWitness>],
) -> Result<IvcProof, IvcError> {
    prove_ivc_with_witnesses(gp, pp, s, steps).map(|(proof, _)| proof)
}

/// Prove a native IVC run, additionally returning the final accumulator's witnesses.
///
/// The witnesses are secret (the prover keeps them across steps); they are surfaced
/// here so the final-compression layer (`superneo-snark`) can compress the accumulator
/// it produced. They are *not* part of [`IvcProof`].
pub fn prove_ivc_with_witnesses(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    steps: &[Vec<CcsWitness>],
) -> Result<(IvcProof, Vec<CeWitness>), IvcError> {
    let mut acc = genesis_instances(gp, pp, s.t());
    let mut acc_wit = genesis_witnesses(gp);
    let mut digest = GENESIS_DIGEST;
    let mut tr = Transcript::new(IVC_DOMAIN);

    let mut fresh_all = Vec::with_capacity(steps.len());
    let mut proofs = Vec::with_capacity(steps.len());

    for (step_idx, wits) in steps.iter().enumerate() {
        let fresh: Vec<CcsInstance> = wits
            .iter()
            .map(|w| CcsInstance {
                c: commit_witness(pp, &w.z),
            })
            .collect();
        let (next_acc, next_wit, proof) = fold(&mut tr, gp, pp, s, &fresh, wits, &acc, &acc_wit)
            .map_err(|e| IvcError::Step {
                step: step_idx,
                source: e,
            })?;
        acc = next_acc;
        acc_wit = next_wit;
        digest = fold_digest(&digest, step_idx, &acc);
        fresh_all.push(fresh);
        proofs.push(proof);
    }

    Ok((
        IvcProof {
            fresh: fresh_all,
            fold_proofs: proofs,
            final_acc: acc,
            final_digest: digest,
        },
        acc_wit,
    ))
}

/// Verify a native IVC run: replay each `verify_fold` and the digest chain, and check
/// the final accumulator and digest match the proof.
pub fn verify_ivc(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    proof: &IvcProof,
) -> Result<(), IvcError> {
    if proof.fresh.len() != proof.fold_proofs.len() {
        return Err(IvcError::AccumulatorMismatch);
    }
    let mut acc = genesis_instances(gp, pp, s.t());
    let mut digest = GENESIS_DIGEST;
    let mut tr = Transcript::new(IVC_DOMAIN);

    for (step_idx, (fresh, fp)) in proof.fresh.iter().zip(proof.fold_proofs.iter()).enumerate() {
        acc = verify_fold(&mut tr, gp, s, fresh, &acc, fp).map_err(|e| IvcError::Step {
            step: step_idx,
            source: e,
        })?;
        digest = fold_digest(&digest, step_idx, &acc);
    }

    if acc != proof.final_acc {
        return Err(IvcError::AccumulatorMismatch);
    }
    if digest != proof.final_digest {
        return Err(IvcError::DigestMismatch);
    }
    Ok(())
}

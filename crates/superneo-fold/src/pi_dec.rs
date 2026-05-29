//! Π_DEC: the decomposition reduction of knowledge `CE(B) → CE(b)^k` (§7.5, Theorem 9).
//!
//! Split the high-norm witness `ẑ` (with `‖ẑ‖∞ < B = b^k`) into `k` low-norm pieces
//! `ẑ = Σ b^{i-1} ẑ_i`, `‖ẑ_i‖∞ < b`. The verifier checks the two additivity relations
//! `c = Σ b^{i-1} c_i` and `y_j = Σ b^{i-1} y_{i,j}`; everything else is derived.

use superneo_commit::PublicParams;
use superneo_field::fp::Fp;
use superneo_ring::decomp::split_b;

use crate::error::FoldError;
use crate::transcript::Transcript;
use crate::types::{
    compute_evals, embed_witness, flatten_ring, CcsStructure, CeInstance, CeWitness, GlobalParams,
    PiDecProof,
};

fn absorb_children(tr: &mut Transcript, proof: &PiDecProof) {
    for c in &proof.child_commitments {
        tr.absorb_commitment(b"pidec/c", c);
    }
    for row in &proof.child_evals {
        for y in row {
            for c in y.coeffs() {
                tr.absorb_ext(b"pidec/y", *c);
            }
        }
    }
}

/// Π_DEC prover.
pub fn pi_dec_prove(
    tr: &mut Transcript,
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    ce: &CeInstance,
    ce_wit: &CeWitness,
) -> Result<(Vec<CeInstance>, Vec<CeWitness>, PiDecProof), FoldError> {
    // Decompose the flattened ring witness into k base-b pieces.
    let flat = flatten_ring(&ce_wit.z_ring);
    let parts = split_b(&flat, gp.b, gp.k)
        .map_err(|e| FoldError::Reduction(format!("Π_DEC split_b: {e}")))?;

    let mut child_wit = Vec::with_capacity(gp.k);
    let mut child_commitments = Vec::with_capacity(gp.k);
    let mut child_evals = Vec::with_capacity(gp.k);
    for part in &parts {
        let zr = embed_witness(part);
        let c = pp
            .commit(&zr)
            .map_err(|e| FoldError::Reduction(format!("Π_DEC commit: {e}")))?;
        let y = compute_evals(s, &zr, &ce.r);
        child_commitments.push(c);
        child_evals.push(y);
        child_wit.push(CeWitness { z_ring: zr });
    }

    let proof = PiDecProof {
        child_commitments,
        child_evals,
    };
    absorb_children(tr, &proof);

    let out_inst = (0..gp.k)
        .map(|i| CeInstance {
            c: proof.child_commitments[i].clone(),
            r: ce.r.clone(),
            y: proof.child_evals[i].clone(),
        })
        .collect();
    Ok((out_inst, child_wit, proof))
}

/// Π_DEC verifier: checks the two additivity relations and reconstructs the children.
pub fn pi_dec_verify(
    tr: &mut Transcript,
    gp: &GlobalParams,
    ce: &CeInstance,
    proof: &PiDecProof,
) -> Result<Vec<CeInstance>, FoldError> {
    let t = ce.y.len();
    if proof.child_commitments.len() != gp.k || proof.child_evals.len() != gp.k {
        return Err(FoldError::Malformed(format!(
            "Π_DEC expects {} children",
            gp.k
        )));
    }
    if proof.child_evals.iter().any(|row| row.len() != t) {
        return Err(FoldError::Malformed(
            "Π_DEC child eval arity mismatch".into(),
        ));
    }
    absorb_children(tr, proof);

    // c =? Σ b^{i-1} c_i.
    let kappa = ce.c.dim();
    let mut c_check = superneo_commit::Commitment::zero(kappa);
    let mut weight = Fp::ONE;
    let b_fp = Fp::new(gp.b);
    for ci in &proof.child_commitments {
        c_check = c_check.add(&ci.scale_fp(weight));
        weight *= b_fp;
    }
    if c_check != ce.c {
        return Err(FoldError::Reduction(
            "Π_DEC commitment additivity c = Σ bⁱ⁻¹·cᵢ failed".into(),
        ));
    }

    // y_j =? Σ b^{i-1} y_{i,j}.
    for j in 0..t {
        let mut y_check = superneo_ring::RingK::ZERO;
        let mut w = Fp::ONE;
        for row in &proof.child_evals {
            y_check = y_check.add(&row[j].scale_fp(w));
            w *= b_fp;
        }
        if y_check != ce.y[j] {
            return Err(FoldError::Reduction(format!(
                "Π_DEC evaluation additivity failed at j={j}"
            )));
        }
    }

    Ok((0..gp.k)
        .map(|i| CeInstance {
            c: proof.child_commitments[i].clone(),
            r: ce.r.clone(),
            y: proof.child_evals[i].clone(),
        })
        .collect())
}

//! In-circuit Π_DEC verifier (M7, Phase P3): the decomposition additivity checks.
//!
//! The native verifier ([`superneo_fold::pi_dec::pi_dec_verify`]) accepts iff the `k`
//! child instances recompose the parent CE instance under the fixed base-`b` weights:
//!
//! ```text
//! c   =? Σ_{i∈[k]} b^{i-1} · c_i        (commitment additivity, over R_F)
//! y_j =? Σ_{i∈[k]} b^{i-1} · y_{i,j}    (evaluation additivity,  over R_K, each j)
//! ```
//!
//! Both checks are *linear*: the weights `b^{i-1}` are base-field constants, so the
//! reconstruction is a constant-scaled sum of advice (the child commitments / evals) tied
//! by equality to the public parent. No multiplication gadgets and no Fiat–Shamir
//! challenges are involved — this is the cheapest phase, and a circuit built from a real
//! proof is satisfied iff `pi_dec_verify` accepts it.
//!
//! Commitments live over the base ring `R_F` (field coefficients), handled with plain
//! [`Lc`] wires; evaluation claims live over `R_K` (extension coefficients), handled with
//! the [`KVar`] gadget whose `scale` by a field constant is free.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::{CeInstance, GlobalParams, PiDecProof};
use superneo_ring::D;

use super::cs::{ConstraintSystem, Lc};
use super::gadgets::KVar;

/// The base-`b` weights `b^0, b^1, …, b^{k-1}` as field constants.
fn base_weights(b: u64, k: usize) -> Vec<Fp> {
    let b_fp = Fp::new(b);
    let mut weights = Vec::with_capacity(k);
    let mut w = Fp::ONE;
    for _ in 0..k {
        weights.push(w);
        w *= b_fp;
    }
    weights
}

/// Build the in-circuit Π_DEC verifier from a real proof and its parent CE instance.
///
/// The parent `(c, y)` is the public statement (pinned as constant wires); the `k` child
/// commitments and evaluations are the proof (advice). The system
/// [`is_satisfied`](ConstraintSystem::is_satisfied) iff
/// [`pi_dec_verify`](superneo_fold::pi_dec::pi_dec_verify) accepts the same proof.
pub fn build_pi_dec_verifier_circuit(
    gp: &GlobalParams,
    ce: &CeInstance,
    proof: &PiDecProof,
) -> ConstraintSystem {
    let t = ce.y.len();
    let kappa = ce.c.dim();
    let weights = base_weights(gp.b, gp.k);

    let mut cs = ConstraintSystem::new();

    // Commitment additivity: c[s][l] == Σ_i b^{i-1} · c_i[s][l] (field-level, linear).
    for s in 0..kappa {
        for l in 0..D {
            let mut acc = Lc::zero();
            for (i, &wi) in weights.iter().enumerate() {
                let coeff = proof.child_commitments[i].0[s].coeffs()[l];
                let var = cs.alloc(coeff);
                acc = acc.add(&Lc::from_var(var).scale(wi));
            }
            let parent = cs.constant(ce.c.0[s].coeffs()[l]);
            cs.assert_eq(&parent, &acc);
        }
    }

    // Evaluation additivity: y[j][l] == Σ_i b^{i-1} · y_{i,j}[l] (K-level, linear).
    for j in 0..t {
        for l in 0..D {
            let mut acc = KVar::constant(&cs, Ext2::ZERO);
            for (i, &wi) in weights.iter().enumerate() {
                let kv = KVar::alloc(&mut cs, proof.child_evals[i][j].coeff(l));
                acc = acc.add(&kv.scale(wi));
            }
            let parent = KVar::constant(&cs, ce.y[j].coeff(l));
            parent.assert_eq(&mut cs, &acc);
        }
    }

    cs
}

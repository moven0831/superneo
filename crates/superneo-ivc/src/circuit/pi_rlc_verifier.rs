//! In-circuit Π_RLC verifier (M7, Phase P2): the ring linear combination.
//!
//! The native verifier ([`superneo_fold::pi_rlc::pi_rlc_verify`]) recomputes the single
//! output CE instance from the `K+k` inputs under the strong-sampling challenges
//! `ρ_i ∈ R_F`:
//!
//! ```text
//! c_out   = Σ_i ρ_i · c_i        (commitment combine, over R_F per κ-slot)
//! y_out_j = Σ_i ρ_i · y_{i,j}    (evaluation combine, over R_K via mul_rf, per j)
//! ```
//!
//! This module synthesizes that combination over the [`RFVar`]/[`RKVar`] ring gadgets and
//! asserts it equals the claimed output. The `ρ_i` are advice (Fiat–Shamir is native in
//! this PoC); the inputs and the output are the public statement (pinned constants). A
//! circuit built from a real proof is satisfied iff the combination is the honest one, so
//! perturbing a `ρ` (which the native verifier rejects as a challenge mismatch) breaks it.
//!
//! This is the heavy phase: each ring product is `D²` (`R_F`) or `2D²` (`R_K`) base mults,
//! so the test runs at scaled-down `κ`, `t`, and input count.

use superneo_fold::{CeInstance, PiRlcProof};

use super::cs::ConstraintSystem;
use super::ring_gadgets::{RFVar, RKVar};

/// Build the in-circuit Π_RLC verifier from the input instances, a real proof (its `ρ`
/// challenges), and the claimed output instance. The system
/// [`is_satisfied`](ConstraintSystem::is_satisfied) iff the output is the honest ring
/// linear combination of the inputs under `ρ` — i.e. iff
/// [`pi_rlc_verify`](superneo_fold::pi_rlc::pi_rlc_verify) would accept these `ρ`.
pub fn build_pi_rlc_verifier_circuit(
    ce: &[CeInstance],
    proof: &PiRlcProof,
    out: &CeInstance,
) -> ConstraintSystem {
    assert!(!ce.is_empty(), "Π_RLC needs at least one input instance");
    let kappa = ce[0].c.dim();
    let t = ce[0].y.len();
    let mut cs = ConstraintSystem::new();

    // The ρ challenges as advice (one R_F element each).
    let rhos: Vec<RFVar> = proof.rhos.iter().map(|rho| RFVar::alloc(&mut cs, rho)).collect();

    // Commitment additivity: c_out[s] == Σ_i ρ_i · c_i[s].
    for s in 0..kappa {
        let mut acc: Option<RFVar> = None;
        for (inst, rho) in ce.iter().zip(&rhos) {
            let ci = RFVar::constant(&cs, &inst.c.0[s]);
            let term = ci.mul(&mut cs, rho);
            acc = Some(match acc {
                None => term,
                Some(a) => a.add(&term),
            });
        }
        let parent = RFVar::constant(&cs, &out.c.0[s]);
        acc.expect("at least one input").assert_eq(&mut cs, &parent);
    }

    // Evaluation additivity: y_out[j] == Σ_i ρ_i · y_{i,j}.
    for j in 0..t {
        let mut acc: Option<RKVar> = None;
        for (inst, rho) in ce.iter().zip(&rhos) {
            let yij = RKVar::constant(&cs, &inst.y[j]);
            let term = yij.mul_rf(&mut cs, rho);
            acc = Some(match acc {
                None => term,
                Some(a) => a.add(&term),
            });
        }
        let parent = RKVar::constant(&cs, &out.y[j]);
        acc.expect("at least one input").assert_eq(&mut cs, &parent);
    }

    cs
}

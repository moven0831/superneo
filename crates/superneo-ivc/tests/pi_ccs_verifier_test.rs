//! Phase P1 of the recursive verifier (M7): the in-circuit Π_CCS `Q(r')` reconstruction.
//!
//! Each test builds the verifier circuit from a *real* `pi_ccs_prove` proof and its
//! `pi_ccs_verify_traced` trace, then asserts the circuit is satisfied **iff** the native
//! Π_CCS verifier accepts the same proof — closing the documented "final-`Q(r')`
//! reconstruction" residual. Tampering an evaluation or a round polynomial must break both.

use superneo_commit::Commitment;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::{
    pi_ccs_prove, pi_ccs_verify, pi_ccs_verify_traced, CcsInstance, CcsStructure, CcsWitness,
    CeInstance, CeWitness, GlobalParams, Transcript,
};
use superneo_ivc::circuit::{build_pi_ccs_verifier_circuit, ConstraintSystem, Lc};
use superneo_ring::{RingElem, RingK};

const KAPPA: usize = 18;
const DOMAIN: &[u8] = b"test/piccs/circuit";

/// A real R1CS structure `x·y = z` with a satisfying, low-norm witness (`‖z‖∞ ≤ 1 < b`).
/// Built with the in-circuit constraint builder, so it finalizes into the `t = 4` CCS
/// shape `fold`/`pi_ccs` consume, with `m` a multiple of the ring degree.
fn mul_structure_and_witness(gp: &GlobalParams) -> (CcsStructure, Vec<CcsInstance>, Vec<CcsWitness>) {
    let mut b = ConstraintSystem::new();
    let x = b.alloc(Fp::ONE);
    let y = b.alloc(Fp::ONE);
    let _z = b.mul(&Lc::from_var(x), &Lc::from_var(y));
    let (s, zf) = b.finalize();
    assert_eq!(zf.len(), gp.n_f, "PoC witness length must equal n_f");
    assert_eq!(s.t(), 4);
    // Π_CCS does not inspect the commitment (that is Π_RLC/Π_DEC), so a zero placeholder
    // is sufficient for the reconstruction check.
    let fresh = vec![CcsInstance {
        c: Commitment::zero(KAPPA),
    }];
    let fresh_wit = vec![CcsWitness { z: zf }];
    (s, fresh, fresh_wit)
}

/// The genesis (trivial, consistent) carried accumulator: a single CE instance with
/// `(c, r, y) = 0` and an all-zero ring witness.
fn trivial_carried(gp: &GlobalParams, t: usize) -> (Vec<CeInstance>, Vec<CeWitness>) {
    let inst = vec![CeInstance {
        c: Commitment::zero(KAPPA),
        r: vec![Ext2::ZERO; gp.log_m],
        y: vec![RingK::ZERO; t],
    }];
    let wit = vec![CeWitness {
        z_ring: vec![RingElem::ZERO; gp.n_r],
    }];
    (inst, wit)
}

#[test]
fn pi_ccs_verifier_circuit_accepts_honest() {
    let gp = GlobalParams::goldilocks(1);
    let (s, fresh, fresh_wit) = mul_structure_and_witness(&gp);
    let (carried, carried_wit) = trivial_carried(&gp, s.t());

    let mut tp = Transcript::new(DOMAIN);
    let (_out, _outw, proof) =
        pi_ccs_prove(&mut tp, &gp, &s, &fresh, &fresh_wit, &carried, &carried_wit);

    let mut tv = Transcript::new(DOMAIN);
    let (_inst, trace) = pi_ccs_verify_traced(&mut tv, &gp, &s, &fresh, &carried, &proof)
        .expect("honest Π_CCS proof must verify natively");

    let cs = build_pi_ccs_verifier_circuit(&gp, &s, &fresh, &carried, &proof, &trace);
    assert!(
        cs.is_satisfied(),
        "honest Π_CCS proof must satisfy the verifier circuit"
    );
    println!(
        "Π_CCS verifier circuit (K=1, k=1 trivial): {} vars, {} constraints",
        cs.num_vars(),
        cs.num_constraints()
    );
}

#[test]
fn pi_ccs_verifier_circuit_rejects_tampered_eval() {
    let gp = GlobalParams::goldilocks(1);
    let (s, fresh, fresh_wit) = mul_structure_and_witness(&gp);
    let (carried, carried_wit) = trivial_carried(&gp, s.t());

    let mut tp = Transcript::new(DOMAIN);
    let (_o, _ow, proof) =
        pi_ccs_prove(&mut tp, &gp, &s, &fresh, &fresh_wit, &carried, &carried_wit);
    let mut tv = Transcript::new(DOMAIN);
    let (_i, trace) = pi_ccs_verify_traced(&mut tv, &gp, &s, &fresh, &carried, &proof).expect("verify");

    // Corrupt the constant term of a claimed evaluation (matrix j = 1, used by the F term).
    // The output evals are absorbed AFTER the sum-check, so the FS-derived (α, γ, r') trace
    // is unchanged — the same situation the native verifier faces, so reusing the honest
    // trace for the circuit is faithful.
    let mut bad = proof.clone();
    let mut coeffs = *bad.evals[0][1].coeffs();
    coeffs[0] += Ext2::ONE;
    bad.evals[0][1] = RingK::from_coeffs(coeffs);

    let mut tv2 = Transcript::new(DOMAIN);
    assert!(
        pi_ccs_verify(&mut tv2, &gp, &s, &fresh, &carried, &bad).is_err(),
        "native Π_CCS must reject a tampered evaluation"
    );
    let cs = build_pi_ccs_verifier_circuit(&gp, &s, &fresh, &carried, &bad, &trace);
    assert!(
        !cs.is_satisfied(),
        "tampered evaluation must break the verifier circuit's Q(r') tie"
    );
}

#[test]
fn pi_ccs_verifier_circuit_rejects_tampered_round_poly() {
    let gp = GlobalParams::goldilocks(1);
    let (s, fresh, fresh_wit) = mul_structure_and_witness(&gp);
    let (carried, carried_wit) = trivial_carried(&gp, s.t());

    let mut tp = Transcript::new(DOMAIN);
    let (_o, _ow, proof) =
        pi_ccs_prove(&mut tp, &gp, &s, &fresh, &fresh_wit, &carried, &carried_wit);
    let mut tv = Transcript::new(DOMAIN);
    let (_i, trace) = pi_ccs_verify_traced(&mut tv, &gp, &s, &fresh, &carried, &proof).expect("verify");

    // Corrupt the first round polynomial: natively this breaks `g(0)+g(1) = T` before the
    // FS-derived point diverges; in-circuit the same round invariant fails.
    let mut bad = proof.clone();
    bad.sumcheck.rounds[0].0[0] += Ext2::ONE;

    let mut tv2 = Transcript::new(DOMAIN);
    assert!(
        pi_ccs_verify(&mut tv2, &gp, &s, &fresh, &carried, &bad).is_err(),
        "native Π_CCS must reject a tampered round polynomial"
    );
    let cs = build_pi_ccs_verifier_circuit(&gp, &s, &fresh, &carried, &bad, &trace);
    assert!(
        !cs.is_satisfied(),
        "tampered round polynomial must break the verifier circuit's sum-check chain"
    );
}

#[test]
fn pi_ccs_verifier_circuit_matches_native_nontrivial_carried() {
    // Exercise the Eval term with a genuine (non-trivial) carried accumulator: a first
    // Π_CCS round produces consistent CE instances + witnesses, which a second round folds.
    let gp = GlobalParams::goldilocks(1);
    let (s, fresh, fresh_wit) = mul_structure_and_witness(&gp);
    let (triv, triv_wit) = trivial_carried(&gp, s.t());

    let mut t1 = Transcript::new(DOMAIN);
    let (carried, carried_wit, _p1) =
        pi_ccs_prove(&mut t1, &gp, &s, &fresh, &fresh_wit, &triv, &triv_wit);
    assert_eq!(carried.len(), 2, "round 1 outputs K + k = 2 CE instances");

    let mut t2p = Transcript::new(DOMAIN);
    let (_o, _ow, proof) =
        pi_ccs_prove(&mut t2p, &gp, &s, &fresh, &fresh_wit, &carried, &carried_wit);
    let mut t2v = Transcript::new(DOMAIN);
    let (_i, trace) = pi_ccs_verify_traced(&mut t2v, &gp, &s, &fresh, &carried, &proof)
        .expect("honest round-2 proof must verify natively");

    let cs = build_pi_ccs_verifier_circuit(&gp, &s, &fresh, &carried, &proof, &trace);
    assert!(
        cs.is_satisfied(),
        "non-trivial carried: circuit must accept the honest proof"
    );
    println!(
        "Π_CCS verifier circuit (K=1, k=2 non-trivial): {} vars, {} constraints",
        cs.num_vars(),
        cs.num_constraints()
    );
}

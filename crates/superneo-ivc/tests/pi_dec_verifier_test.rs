//! Phase P3 of the recursive verifier (M7): the in-circuit Π_DEC verifier.
//!
//! Builds the Π_DEC verifier circuit from a *real* `pi_dec_prove` proof on a non-trivial
//! CE instance and asserts the circuit is satisfied **iff** the native `pi_dec_verify`
//! accepts. The two additivity relations are linear (constant base-`b` weights), so this
//! phase has no challenges and no multiplication gadgets. Tampering a child commitment or
//! a child evaluation must break both the native check and the circuit.

use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::types::{compute_evals, embed_witness};
use superneo_fold::{pi_dec_prove, pi_dec_verify, CeInstance, CeWitness, GlobalParams, Transcript};
use superneo_ivc::circuit::{build_pi_dec_verifier_circuit, ConstraintSystem, Lc};
use superneo_ring::{RingElem, RingK};

const KAPPA: usize = 18;
const N_R: usize = 1;
const DOMAIN: &[u8] = b"test/pidec/circuit";

/// A `t = 4` structure built with the constraint builder (its witness is unused here).
fn structure() -> superneo_fold::CcsStructure {
    let mut b = ConstraintSystem::new();
    let x = b.alloc(Fp::ONE);
    let y = b.alloc(Fp::ONE);
    let _z = b.mul(&Lc::from_var(x), &Lc::from_var(y));
    b.finalize().0
}

/// A non-trivial, high-norm-but-decomposable CE instance and its witness: values are well
/// under `B = b^k = 16384`, so `split_b` produces several non-zero base-`b` pieces.
fn nontrivial_ce(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &superneo_fold::CcsStructure,
) -> (CeInstance, CeWitness) {
    let zf: Vec<Fp> = (0..gp.n_f)
        .map(|i| Fp::new(((i * 37 + 11) % 4000) as u64))
        .collect();
    let zr = embed_witness(&zf);
    let r: Vec<Ext2> = (0..gp.log_m)
        .map(|i| Ext2::from_base(Fp::new(i as u64 + 2)))
        .collect();
    let y = compute_evals(s, &zr, &r);
    let c = pp.commit(&zr).expect("commit");
    (CeInstance { c, r, y }, CeWitness { z_ring: zr })
}

fn setup() -> (GlobalParams, PublicParams, superneo_fold::CcsStructure) {
    (
        GlobalParams::goldilocks(N_R),
        PublicParams::setup_seeded([3u8; 32], KAPPA, N_R),
        structure(),
    )
}

#[test]
fn pi_dec_verifier_circuit_accepts_honest() {
    let (gp, pp, s) = setup();
    let (ce, ce_wit) = nontrivial_ce(&gp, &pp, &s);

    let mut tp = Transcript::new(DOMAIN);
    let (_out, _outw, proof) =
        pi_dec_prove(&mut tp, &gp, &pp, &s, &ce, &ce_wit).expect("honest Π_DEC prove");

    let mut tv = Transcript::new(DOMAIN);
    pi_dec_verify(&mut tv, &gp, &ce, &proof).expect("honest Π_DEC proof must verify natively");

    let cs = build_pi_dec_verifier_circuit(&gp, &ce, &proof);
    assert!(
        cs.is_satisfied(),
        "honest Π_DEC proof must satisfy the verifier circuit"
    );
    println!(
        "Π_DEC verifier circuit (k={}): {} vars, {} constraints",
        gp.k,
        cs.num_vars(),
        cs.num_constraints()
    );
}

#[test]
fn pi_dec_verifier_circuit_rejects_tampered_child_eval() {
    let (gp, pp, s) = setup();
    let (ce, ce_wit) = nontrivial_ce(&gp, &pp, &s);
    let mut tp = Transcript::new(DOMAIN);
    let (_o, _ow, proof) = pi_dec_prove(&mut tp, &gp, &pp, &s, &ce, &ce_wit).expect("prove");

    // Break evaluation additivity: perturb one coefficient of a child evaluation.
    let mut bad = proof.clone();
    let mut coeffs = *bad.child_evals[0][0].coeffs();
    coeffs[0] += Ext2::ONE;
    bad.child_evals[0][0] = RingK::from_coeffs(coeffs);

    let mut tv = Transcript::new(DOMAIN);
    assert!(
        pi_dec_verify(&mut tv, &gp, &ce, &bad).is_err(),
        "native Π_DEC must reject a tampered child evaluation"
    );
    let cs = build_pi_dec_verifier_circuit(&gp, &ce, &bad);
    assert!(
        !cs.is_satisfied(),
        "tampered child evaluation must break evaluation additivity"
    );
}

#[test]
fn pi_dec_verifier_circuit_rejects_tampered_child_commitment() {
    let (gp, pp, s) = setup();
    let (ce, ce_wit) = nontrivial_ce(&gp, &pp, &s);
    let mut tp = Transcript::new(DOMAIN);
    let (_o, _ow, proof) = pi_dec_prove(&mut tp, &gp, &pp, &s, &ce, &ce_wit).expect("prove");

    // Break commitment additivity: perturb one coefficient of a child commitment slot.
    let mut bad = proof.clone();
    let mut coeffs = *bad.child_commitments[0].0[0].coeffs();
    coeffs[0] += Fp::ONE;
    bad.child_commitments[0].0[0] = RingElem::from_coeffs(coeffs);

    let mut tv = Transcript::new(DOMAIN);
    assert!(
        pi_dec_verify(&mut tv, &gp, &ce, &bad).is_err(),
        "native Π_DEC must reject a tampered child commitment"
    );
    let cs = build_pi_dec_verifier_circuit(&gp, &ce, &bad);
    assert!(
        !cs.is_satisfied(),
        "tampered child commitment must break commitment additivity"
    );
}

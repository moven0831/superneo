//! Phase P2 of the recursive verifier (M7): the in-circuit Π_RLC verifier.
//!
//! Builds the Π_RLC verifier circuit from a *real* `pi_rlc_prove` proof and asserts the
//! ring linear combination `c_out = Σ ρ_i·c_i`, `y_out = Σ ρ_i·y_{i,j}` is reproduced
//! in-circuit. Run at scaled-down parameters (κ = 1, t = 1, two inputs) because each ring
//! product is `D²`/`2D²` base multiplications — the dominant recursion cost. Perturbing a
//! ρ challenge (which the native verifier rejects as a challenge mismatch) breaks the
//! circuit's combination.

use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::types::{compute_evals, embed_witness};
use superneo_fold::{
    pi_rlc_prove, pi_rlc_verify, CcsStructure, CeInstance, CeWitness, GlobalParams, Transcript,
};
use superneo_ivc::circuit::build_pi_rlc_verifier_circuit;
use superneo_ring::RingElem;

const KAPPA: usize = 1; // scaled down: each ring product is D² base mults
const N_R: usize = 1;
const N_IN: usize = 2; // K + k inputs
const DOMAIN: &[u8] = b"test/pirlc/circuit";

/// `N_IN` consistent CE instances (real commitments + evals) sharing one evaluation point,
/// and their witnesses. The identity structure gives `t = 1` to keep the test small.
fn inputs(
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
) -> (Vec<CeInstance>, Vec<CeWitness>) {
    let r: Vec<Ext2> = (0..gp.log_m)
        .map(|i| Ext2::from_base(Fp::new(i as u64 + 1)))
        .collect();
    let mut insts = Vec::with_capacity(N_IN);
    let mut wits = Vec::with_capacity(N_IN);
    for seed in 0..N_IN {
        let zf: Vec<Fp> = (0..gp.n_f)
            .map(|i| Fp::new(((i * 13 + seed * 7 + 3) % 100) as u64))
            .collect();
        let zr = embed_witness(&zf);
        let c = pp.commit(&zr).expect("commit");
        let y = compute_evals(s, &zr, &r);
        insts.push(CeInstance { c, r: r.clone(), y });
        wits.push(CeWitness { z_ring: zr });
    }
    (insts, wits)
}

fn setup() -> (GlobalParams, PublicParams, CcsStructure) {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([4u8; 32], KAPPA, N_R);
    let s = CcsStructure::identity(gp.m);
    (gp, pp, s)
}

#[test]
fn pi_rlc_verifier_circuit_accepts_honest() {
    let (gp, pp, s) = setup();
    let (ce, wit) = inputs(&gp, &pp, &s);

    let mut tp = Transcript::new(DOMAIN);
    let (out, _out_wit, proof) = pi_rlc_prove(&mut tp, &gp, &ce, &wit);

    let mut tv = Transcript::new(DOMAIN);
    let out_v = pi_rlc_verify(&mut tv, &gp, &ce, &proof).expect("honest Π_RLC must verify");
    assert_eq!(out_v, out, "prover and verifier agree on the combined instance");

    let cs = build_pi_rlc_verifier_circuit(&ce, &proof, &out);
    assert!(
        cs.is_satisfied(),
        "honest Π_RLC combination must satisfy the verifier circuit"
    );
    println!(
        "Π_RLC verifier circuit (κ={KAPPA}, t={}, {N_IN} inputs): {} vars, {} constraints",
        out.y.len(),
        cs.num_vars(),
        cs.num_constraints()
    );
}

#[test]
fn pi_rlc_verifier_circuit_rejects_tampered_rho() {
    let (gp, pp, s) = setup();
    let (ce, wit) = inputs(&gp, &pp, &s);
    let mut tp = Transcript::new(DOMAIN);
    let (out, _ow, proof) = pi_rlc_prove(&mut tp, &gp, &ce, &wit);

    // Perturb one coefficient of a ρ challenge.
    let mut bad = proof.clone();
    let mut coeffs = *bad.rhos[0].coeffs();
    coeffs[0] += Fp::ONE;
    bad.rhos[0] = RingElem::from_coeffs(coeffs);

    let mut tv = Transcript::new(DOMAIN);
    assert!(
        pi_rlc_verify(&mut tv, &gp, &ce, &bad).is_err(),
        "native Π_RLC must reject mismatched ρ challenges"
    );
    // The honest output no longer equals Σ ρ'·c_i, so the circuit is unsatisfiable.
    let cs = build_pi_rlc_verifier_circuit(&ce, &bad, &out);
    assert!(
        !cs.is_satisfied(),
        "a perturbed ρ must break the in-circuit ring combination"
    );
}

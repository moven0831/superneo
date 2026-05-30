//! M8 integration: the full pipeline verifies for both the identity (recursion-overhead)
//! workload and a real multiplication-gate workload, and the demonstration circuits
//! compile to satisfied CCS instances.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::fp::Fp;
use superneo_fold::{CcsStructure, CcsWitness, GlobalParams};
use superneo_snark::pipeline::{
    counter_circuit, fibonacci_circuit, mul_gate_witness, multiplication_structure,
    multiplier_circuit, run_pipeline,
};

const KAPPA: usize = 18;
const N_R: usize = 1;
const STEPS: usize = 4;

fn setup() -> (GlobalParams, PublicParams) {
    (
        GlobalParams::goldilocks(N_R),
        PublicParams::setup_seeded([13u8; 32], KAPPA, N_R),
    )
}

fn low_norm_field(r: &mut ChaCha8Rng, n: usize) -> Vec<Fp> {
    (0..n)
        .map(|_| Fp::from_i64((r.next_u64() % 3) as i64 - 1))
        .collect()
}

#[test]
fn identity_pipeline_verifies() {
    let (gp, pp) = setup();
    let s = CcsStructure::identity(gp.m);
    let mut r = ChaCha8Rng::seed_from_u64(0xABCD_0001);
    let steps: Vec<Vec<CcsWitness>> = (0..STEPS)
        .map(|_| {
            vec![CcsWitness {
                z: low_norm_field(&mut r, gp.n_f),
            }]
        })
        .collect();
    let rep = run_pipeline(&gp, &pp, &s, &steps).expect("pipeline");
    assert!(rep.verified, "identity pipeline must verify");
    assert_eq!(rep.steps, STEPS);
    assert_eq!(rep.k, gp.k);
    assert!(rep.proof_bytes > 0);
}

#[test]
fn multiplication_pipeline_verifies() {
    let (gp, pp) = setup();
    let s = multiplication_structure(gp.m);
    let steps: Vec<Vec<CcsWitness>> = (0..STEPS)
        .map(|i| vec![mul_gate_witness(gp.m, (i % 2) as u64)])
        .collect();
    let rep = run_pipeline(&gp, &pp, &s, &steps).expect("pipeline");
    assert!(rep.verified, "x·y=z pipeline must verify");
}

#[test]
fn demonstration_circuits_are_satisfied_ccs() {
    for cs in [
        counter_circuit(41),
        fibonacci_circuit(8, 13),
        multiplier_circuit(6, 7),
    ] {
        assert!(cs.is_satisfied(), "demonstration circuit must be satisfied");
        let (structure, z) = cs.finalize();
        // Finalizes into the fold's t=4 CCS shape with a d-aligned dimension.
        assert_eq!(structure.t(), 4);
        assert_eq!(structure.m % superneo_ring::D, 0);
        assert_eq!(z.len(), structure.m);
        assert!(cs.finalized_relation_holds());
    }
}

#[test]
fn counter_witness_is_correct() {
    // out = in + 1 for a concrete input.
    let cs = counter_circuit(100);
    assert!(cs.is_satisfied());
}

#[test]
fn invalid_witness_rejected_end_to_end() {
    // The headline M8 claim: a low-norm but UNSATISFYING witness must not yield a verified
    // run. z0·z1 = 1·1 = 1 ≠ 0 = z2 violates f = X₂·X₃ − X₄, while ‖z‖∞ < b = 2 passes the
    // norm gate (so it is the CCS check, not the norm check, that must catch it).
    use superneo_ivc::{prove_ivc, verify_ivc};
    let (gp, pp) = setup();
    let s = multiplication_structure(gp.m);
    let mut z = vec![Fp::ZERO; gp.n_f];
    z[0] = Fp::ONE;
    z[1] = Fp::ONE;
    z[2] = Fp::ZERO;
    let steps = vec![vec![CcsWitness { z }]];

    // The prover does not pre-check satisfaction; the fold verifier must reject.
    match prove_ivc(&gp, &pp, &s, &steps) {
        Err(_) => {} // rejected already at prove time — acceptable
        Ok(proof) => assert!(
            verify_ivc(&gp, &pp, &s, &proof).is_err(),
            "fold verifier accepted an unsatisfying witness"
        ),
    }
    // End to end, the pipeline must not report a verified run.
    let verified = run_pipeline(&gp, &pp, &s, &steps)
        .map(|r| r.verified)
        .unwrap_or(false);
    assert!(!verified, "unsatisfying witness produced a verified pipeline run");
}

#[test]
fn demonstration_circuits_handle_large_inputs() {
    // Witness values are computed in the field, so inputs near u64::MAX neither overflow
    // nor desync from the field-arithmetic constraints.
    assert!(counter_circuit(u64::MAX).is_satisfied());
    assert!(fibonacci_circuit(u64::MAX, u64::MAX - 1).is_satisfied());
}

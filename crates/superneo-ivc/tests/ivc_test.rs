//! Native IVC: a multi-step run folds and verifies, and tampering with either the
//! digest chain or a step's fold proof is rejected.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::{CcsStructure, CcsWitness, GlobalParams};
use superneo_ivc::{prove_ivc, verify_ivc};

const KAPPA: usize = 18;
const N_R: usize = 1;
const STEPS: usize = 4;

fn low_norm_field(r: &mut ChaCha8Rng, n: usize) -> Vec<Fp> {
    (0..n)
        .map(|_| Fp::from_i64((r.next_u64() % 3) as i64 - 1))
        .collect()
}

fn setup() -> (
    GlobalParams,
    PublicParams,
    CcsStructure,
    Vec<Vec<CcsWitness>>,
) {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([5u8; 32], KAPPA, N_R);
    let s = CcsStructure::identity(gp.m);
    let mut r = ChaCha8Rng::seed_from_u64(0xACE0_1234_5678_9ABC);
    let steps = (0..STEPS)
        .map(|_| {
            vec![CcsWitness {
                z: low_norm_field(&mut r, gp.n_f),
            }]
        })
        .collect();
    (gp, pp, s, steps)
}

#[test]
fn ivc_chain_verifies() {
    let (gp, pp, s, steps) = setup();
    let proof = prove_ivc(&gp, &pp, &s, &steps).expect("prove");
    assert_eq!(proof.fold_proofs.len(), STEPS);
    assert_eq!(proof.final_acc.len(), gp.k);
    verify_ivc(&gp, &pp, &s, &proof).expect("verify");
}

#[test]
fn ivc_rejects_tampered_digest() {
    let (gp, pp, s, steps) = setup();
    let mut proof = prove_ivc(&gp, &pp, &s, &steps).unwrap();
    proof.final_digest[0] ^= 1;
    assert!(verify_ivc(&gp, &pp, &s, &proof).is_err());
}

#[test]
fn ivc_rejects_tampered_fold() {
    let (gp, pp, s, steps) = setup();
    let mut proof = prove_ivc(&gp, &pp, &s, &steps).unwrap();
    proof.fold_proofs[1].ccs.sumcheck.rounds[0].0[0] += Ext2::ONE;
    assert!(verify_ivc(&gp, &pp, &s, &proof).is_err());
}

#[test]
fn ivc_rejects_dropped_step() {
    let (gp, pp, s, steps) = setup();
    let mut proof = prove_ivc(&gp, &pp, &s, &steps).unwrap();
    // Drop the last step's proof but keep its accumulator/digest claim.
    proof.fold_proofs.pop();
    proof.fresh.pop();
    assert!(verify_ivc(&gp, &pp, &s, &proof).is_err());
}

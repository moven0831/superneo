//! End-to-end folding: `fold` / `verify_fold` round-trip, two-step composition
//! (Theorem 3), and tamper rejection.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::types::{commit_witness, compute_evals};
use superneo_fold::{
    fold, verify_fold, CcsInstance, CcsStructure, CcsWitness, CeInstance, CeWitness, GlobalParams,
};
use superneo_ring::{RingElem, D};

const KAPPA: usize = 18;
const N_R: usize = 1;
const CAP_K: usize = 2; // fresh CCS instances per fold

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0xF01D_5CE0_AACC_0DEF)
}

fn low_norm_field(r: &mut ChaCha8Rng, n: usize) -> Vec<Fp> {
    (0..n)
        .map(|_| Fp::from_i64((r.next_u64() % 3) as i64 - 1))
        .collect()
}

fn low_norm_ring(r: &mut ChaCha8Rng, n_r: usize) -> Vec<RingElem> {
    (0..n_r)
        .map(|_| {
            let mut c = [Fp::ZERO; D];
            for x in c.iter_mut() {
                *x = Fp::from_i64((r.next_u64() % 3) as i64 - 1);
            }
            RingElem::from_coeffs(c)
        })
        .collect()
}

fn rand_point(r: &mut ChaCha8Rng, ell: usize) -> Vec<Ext2> {
    (0..ell)
        .map(|_| Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64())))
        .collect()
}

struct Setup {
    gp: GlobalParams,
    pp: PublicParams,
    s: CcsStructure,
}

fn setup() -> Setup {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([11u8; 32], KAPPA, N_R);
    let s = CcsStructure::identity(gp.m);
    Setup { gp, pp, s }
}

/// Build `k` genesis CE instances (all at a shared point) and `K` fresh CCS instances.
fn make_inputs(
    st: &Setup,
    r: &mut ChaCha8Rng,
) -> (
    Vec<CcsInstance>,
    Vec<CcsWitness>,
    Vec<CeInstance>,
    Vec<CeWitness>,
) {
    let point = rand_point(r, st.gp.log_m);
    let mut carried = Vec::new();
    let mut carried_wit = Vec::new();
    for _ in 0..st.gp.k {
        let zr = low_norm_ring(r, st.gp.n_r);
        let c = st.pp.commit(&zr).unwrap();
        let y = compute_evals(&st.s, &zr, &point);
        carried.push(CeInstance {
            c,
            r: point.clone(),
            y,
        });
        carried_wit.push(CeWitness { z_ring: zr });
    }
    let mut fresh = Vec::new();
    let mut fresh_wit = Vec::new();
    for _ in 0..CAP_K {
        let z = low_norm_field(r, st.gp.n_f);
        fresh.push(CcsInstance {
            c: commit_witness(&st.pp, &z),
        });
        fresh_wit.push(CcsWitness { z });
    }
    (fresh, fresh_wit, carried, carried_wit)
}

#[test]
fn fold_round_trip() {
    let st = setup();
    let mut r = rng();
    let (fresh, fresh_wit, carried, carried_wit) = make_inputs(&st, &mut r);

    let mut tp = superneo_fold::Transcript::new(b"superneo/fold/v1");
    let (acc, _acc_wit, proof) = fold(
        &mut tp,
        &st.gp,
        &st.pp,
        &st.s,
        &fresh,
        &fresh_wit,
        &carried,
        &carried_wit,
    )
    .expect("prove");
    assert_eq!(acc.len(), st.gp.k, "accumulator must have k instances");

    let mut tv = superneo_fold::Transcript::new(b"superneo/fold/v1");
    let acc_v = verify_fold(&mut tv, &st.gp, &st.s, &fresh, &carried, &proof).expect("verify");
    assert_eq!(acc, acc_v, "verifier accumulator must match prover's");
}

#[test]
fn fold_two_steps() {
    let st = setup();
    let mut r = rng();
    let (fresh1, fresh1_wit, carried0, carried0_wit) = make_inputs(&st, &mut r);

    let mut tp = superneo_fold::Transcript::new(b"superneo/fold/v1");
    // Step 1.
    let (acc1, acc1_wit, proof1) = fold(
        &mut tp,
        &st.gp,
        &st.pp,
        &st.s,
        &fresh1,
        &fresh1_wit,
        &carried0,
        &carried0_wit,
    )
    .expect("prove step1");
    // Step 2: fold fresh instances against the step-1 accumulator.
    let mut fresh2 = Vec::new();
    let mut fresh2_wit = Vec::new();
    for _ in 0..CAP_K {
        let z = low_norm_field(&mut r, st.gp.n_f);
        fresh2.push(CcsInstance {
            c: commit_witness(&st.pp, &z),
        });
        fresh2_wit.push(CcsWitness { z });
    }
    let (acc2, _acc2_wit, proof2) = fold(
        &mut tp,
        &st.gp,
        &st.pp,
        &st.s,
        &fresh2,
        &fresh2_wit,
        &acc1,
        &acc1_wit,
    )
    .expect("prove step2");

    // Verify both steps on a parallel transcript.
    let mut tv = superneo_fold::Transcript::new(b"superneo/fold/v1");
    let acc1_v =
        verify_fold(&mut tv, &st.gp, &st.s, &fresh1, &carried0, &proof1).expect("verify step1");
    assert_eq!(acc1, acc1_v);
    let acc2_v =
        verify_fold(&mut tv, &st.gp, &st.s, &fresh2, &acc1_v, &proof2).expect("verify step2");
    assert_eq!(acc2, acc2_v);
}

#[test]
fn rejects_tampered_sumcheck() {
    let st = setup();
    let mut r = rng();
    let (fresh, fresh_wit, carried, carried_wit) = make_inputs(&st, &mut r);
    let mut tp = superneo_fold::Transcript::new(b"superneo/fold/v1");
    let (_, _, mut proof) = fold(
        &mut tp,
        &st.gp,
        &st.pp,
        &st.s,
        &fresh,
        &fresh_wit,
        &carried,
        &carried_wit,
    )
    .unwrap();
    proof.ccs.sumcheck.rounds[1].0[0] += Ext2::ONE;
    let mut tv = superneo_fold::Transcript::new(b"superneo/fold/v1");
    assert!(verify_fold(&mut tv, &st.gp, &st.s, &fresh, &carried, &proof).is_err());
}

#[test]
fn rejects_tampered_dec_child() {
    let st = setup();
    let mut r = rng();
    let (fresh, fresh_wit, carried, carried_wit) = make_inputs(&st, &mut r);
    let mut tp = superneo_fold::Transcript::new(b"superneo/fold/v1");
    let (_, _, mut proof) = fold(
        &mut tp,
        &st.gp,
        &st.pp,
        &st.s,
        &fresh,
        &fresh_wit,
        &carried,
        &carried_wit,
    )
    .unwrap();
    // Corrupt the first child commitment so c = Σ bⁱ⁻¹ cᵢ no longer holds.
    proof.dec.child_commitments[0].0[0] = proof.dec.child_commitments[0].0[0] + RingElem::one();
    let mut tv = superneo_fold::Transcript::new(b"superneo/fold/v1");
    assert!(verify_fold(&mut tv, &st.gp, &st.s, &fresh, &carried, &proof).is_err());
}

#[test]
fn rejects_tampered_fresh_commitment() {
    let st = setup();
    let mut r = rng();
    let (mut fresh, fresh_wit, carried, carried_wit) = make_inputs(&st, &mut r);
    let mut tp = superneo_fold::Transcript::new(b"superneo/fold/v1");
    let (_, _, proof) = fold(
        &mut tp,
        &st.gp,
        &st.pp,
        &st.s,
        &fresh,
        &fresh_wit,
        &carried,
        &carried_wit,
    )
    .unwrap();
    // Verifier sees a different fresh commitment than the prover bound.
    fresh[0].c.0[0] = fresh[0].c.0[0] + RingElem::one();
    let mut tv = superneo_fold::Transcript::new(b"superneo/fold/v1");
    assert!(verify_fold(&mut tv, &st.gp, &st.s, &fresh, &carried, &proof).is_err());
}

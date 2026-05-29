//! R1CS-as-CCS folding: exercises the CCS-satisfaction (`F`) path of Π_CCS with a
//! genuine multiplication gate `x·y = p`, encoded with `t=4` matrices
//! `(M_1=I, M_2=A, M_3=B, M_4=C)` and `f = X_2·X_3 − X_4`.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::types::{commit_witness, compute_evals};
use superneo_fold::{
    fold, verify_fold, CcsInstance, CcsStructure, CcsWitness, CeInstance, CeWitness, GlobalParams,
    SparsePoly,
};
use superneo_ring::{RingElem, D};

const KAPPA: usize = 18;
const N_R: usize = 1;

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x1234_ABCD_5678_EF01)
}

/// R1CS for a single multiplication gate `z[0]·z[1] = z[2]`, the rest padding.
fn r1cs_structure(m: usize) -> CcsStructure {
    let mut id = vec![vec![Fp::ZERO; m]; m];
    for (i, row) in id.iter_mut().enumerate() {
        row[i] = Fp::ONE;
    }
    let mut a = vec![vec![Fp::ZERO; m]; m];
    let mut b = vec![vec![Fp::ZERO; m]; m];
    let mut c = vec![vec![Fp::ZERO; m]; m];
    a[0][0] = Fp::ONE; // A·z = z[0]
    b[0][1] = Fp::ONE; // B·z = z[1]
    c[0][2] = Fp::ONE; // C·z = z[2]
    let f = SparsePoly::new(
        4,
        vec![
            (Fp::ONE, vec![0, 1, 1, 0]),  // +X_2·X_3
            (-Fp::ONE, vec![0, 0, 0, 1]), // −X_4
        ],
    );
    CcsStructure::new(vec![id, a, b, c], f)
}

fn witness(m: usize, x: u64, y: u64, p: u64) -> Vec<Fp> {
    let mut z = vec![Fp::ZERO; m];
    z[0] = Fp::new(x);
    z[1] = Fp::new(y);
    z[2] = Fp::new(p);
    z
}

fn carried(
    s: &CcsStructure,
    pp: &PublicParams,
    gp: &GlobalParams,
    r: &mut ChaCha8Rng,
) -> (Vec<CeInstance>, Vec<CeWitness>) {
    let point: Vec<Ext2> = (0..gp.log_m)
        .map(|_| Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64())))
        .collect();
    let mut inst = Vec::new();
    let mut wit = Vec::new();
    for _ in 0..gp.k {
        let mut c = [Fp::ZERO; D];
        for x in c.iter_mut() {
            *x = Fp::from_i64((r.next_u64() % 3) as i64 - 1);
        }
        let zr = vec![RingElem::from_coeffs(c)];
        let comm = pp.commit(&zr).unwrap();
        let y = compute_evals(s, &zr, &point);
        inst.push(CeInstance {
            c: comm,
            r: point.clone(),
            y,
        });
        wit.push(CeWitness { z_ring: zr });
    }
    (inst, wit)
}

#[test]
fn r1cs_valid_witness_folds_and_verifies() {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([21u8; 32], KAPPA, N_R);
    let s = r1cs_structure(gp.m);
    let mut r = rng();

    // Valid: 1 · 1 = 1.
    let z = witness(gp.m, 1, 1, 1);
    let fresh = vec![CcsInstance {
        c: commit_witness(&pp, &z),
    }];
    let fresh_wit = vec![CcsWitness { z }];
    let (carried_i, carried_w) = carried(&s, &pp, &gp, &mut r);

    let mut tp = superneo_fold::Transcript::new(b"r1cs");
    let (acc, _, proof) = fold(
        &mut tp, &gp, &pp, &s, &fresh, &fresh_wit, &carried_i, &carried_w,
    )
    .unwrap();
    let mut tv = superneo_fold::Transcript::new(b"r1cs");
    let acc_v = verify_fold(&mut tv, &gp, &s, &fresh, &carried_i, &proof).unwrap();
    assert_eq!(acc, acc_v);
}

#[test]
fn r1cs_invalid_witness_is_rejected() {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([21u8; 32], KAPPA, N_R);
    let s = r1cs_structure(gp.m);
    let mut r = rng();

    // Invalid: claims 1 · 1 = 0. Low-norm, so NC passes, but F does not vanish.
    let z = witness(gp.m, 1, 1, 0);
    let fresh = vec![CcsInstance {
        c: commit_witness(&pp, &z),
    }];
    let fresh_wit = vec![CcsWitness { z }];
    let (carried_i, carried_w) = carried(&s, &pp, &gp, &mut r);

    let mut tp = superneo_fold::Transcript::new(b"r1cs");
    let (_, _, proof) = fold(
        &mut tp, &gp, &pp, &s, &fresh, &fresh_wit, &carried_i, &carried_w,
    )
    .unwrap();
    // The honest prover's Σ_x Q(x) ≠ T, so the sum-check round-0 check fails.
    let mut tv = superneo_fold::Transcript::new(b"r1cs");
    assert!(verify_fold(&mut tv, &gp, &s, &fresh, &carried_i, &proof).is_err());
}

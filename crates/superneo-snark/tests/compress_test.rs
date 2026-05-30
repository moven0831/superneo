//! Final compression of a CE accumulator: an honest accumulator compresses and
//! verifies; tampering with the public claims, the witness opening, or a sum-check is
//! rejected; and a high-norm witness fails the norm sum-check.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::types::compute_evals;
use superneo_fold::{CcsStructure, CeInstance, CeWitness, GlobalParams};
use superneo_ring::{RingElem, D};
use superneo_snark::{compress, verify};

const KAPPA: usize = 18;
const N_R: usize = 1;

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0xC0FF_EE00_1234_5678)
}

/// A low-norm (`‖·‖∞ < b = 2`) ring witness: coefficients in `{−1, 0, 1}`.
fn low_norm_ring(r: &mut ChaCha8Rng) -> Vec<RingElem> {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::from_i64((r.next_u64() % 3) as i64 - 1);
    }
    vec![RingElem::from_coeffs(c)]
}

struct Setup {
    gp: GlobalParams,
    pp: PublicParams,
    s: CcsStructure,
}

fn setup() -> Setup {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([42u8; 32], KAPPA, N_R);
    let s = CcsStructure::identity(gp.m);
    Setup { gp, pp, s }
}

/// Build a `k`-instance accumulator at a shared point with honest low-norm witnesses.
fn accumulator(
    st: &Setup,
    r: &mut ChaCha8Rng,
    low_norm: bool,
) -> (Vec<CeInstance>, Vec<CeWitness>) {
    let point: Vec<Ext2> = (0..st.gp.log_m)
        .map(|_| Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64())))
        .collect();
    let mut inst = Vec::new();
    let mut wit = Vec::new();
    for idx in 0..st.gp.k {
        let mut zr = low_norm_ring(r);
        if !low_norm && idx == 0 {
            // Inject an out-of-range coefficient (‖·‖∞ = 5 ≥ b).
            let mut c = *zr[0].coeffs();
            c[3] = Fp::from_i64(5);
            zr[0] = RingElem::from_coeffs(c);
        }
        let c = st.pp.commit(&zr).unwrap();
        let y = compute_evals(&st.s, &zr, &point);
        inst.push(CeInstance {
            c,
            r: point.clone(),
            y,
        });
        wit.push(CeWitness { z_ring: zr });
    }
    (inst, wit)
}

#[test]
fn compress_and_verify_round_trip() {
    let st = setup();
    let mut r = rng();
    let (acc, acc_wit) = accumulator(&st, &mut r, true);
    let proof = compress(&st.gp, &st.pp, &st.s, &acc, &acc_wit).expect("compress");
    verify(&st.gp, &st.pp, &st.s, &acc, &proof).expect("verify");
}

#[test]
fn rejects_tampered_eval_claim() {
    let st = setup();
    let mut r = rng();
    let (mut acc, acc_wit) = accumulator(&st, &mut r, true);
    let proof = compress(&st.gp, &st.pp, &st.s, &acc, &acc_wit).unwrap();
    // The verifier is shown a different evaluation claim than was proved.
    acc[0].y[0] = acc[0].y[0].add(&superneo_ring::RingK::from_const(Ext2::ONE));
    assert!(verify(&st.gp, &st.pp, &st.s, &acc, &proof).is_err());
}

#[test]
fn rejects_tampered_commitment_claim() {
    let st = setup();
    let mut r = rng();
    let (mut acc, acc_wit) = accumulator(&st, &mut r, true);
    let proof = compress(&st.gp, &st.pp, &st.s, &acc, &acc_wit).unwrap();
    // Corrupt a public Ajtai commitment coefficient.
    acc[1].c.0[0] = acc[1].c.0[0] + RingElem::one();
    assert!(verify(&st.gp, &st.pp, &st.s, &acc, &proof).is_err());
}

#[test]
fn rejects_tampered_opening() {
    let st = setup();
    let mut r = rng();
    let (acc, acc_wit) = accumulator(&st, &mut r, true);
    let mut proof = compress(&st.gp, &st.pp, &st.s, &acc, &acc_wit).unwrap();
    proof.lin_open.value += Ext2::ONE;
    assert!(verify(&st.gp, &st.pp, &st.s, &acc, &proof).is_err());
}

#[test]
fn rejects_tampered_norm_sumcheck() {
    let st = setup();
    let mut r = rng();
    let (acc, acc_wit) = accumulator(&st, &mut r, true);
    let mut proof = compress(&st.gp, &st.pp, &st.s, &acc, &acc_wit).unwrap();
    proof.norm_rounds[0][0] += Ext2::ONE;
    assert!(verify(&st.gp, &st.pp, &st.s, &acc, &proof).is_err());
}

#[test]
fn rejects_high_norm_witness() {
    let st = setup();
    let mut r = rng();
    let (acc, acc_wit) = accumulator(&st, &mut r, false);
    // The prover can still build a proof, but the norm sum-check no longer sums to 0.
    let proof = compress(&st.gp, &st.pp, &st.s, &acc, &acc_wit).unwrap();
    assert!(verify(&st.gp, &st.pp, &st.s, &acc, &proof).is_err());
}

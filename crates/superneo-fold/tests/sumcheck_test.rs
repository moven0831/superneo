//! Sum-check protocol (Definition 6): completeness, the `Q(r')` reduction, and
//! soundness (wrong claim / tampered round polynomial are rejected).

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::multilinear::mle_eval_ext;
use superneo_fold::sumcheck::{sumcheck_prove, sumcheck_verify};
use superneo_fold::transcript::Transcript;

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x5C_5C_5C_5C_DE_AD_BE_EF)
}

fn rand_ext(r: &mut ChaCha8Rng) -> Ext2 {
    Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64()))
}

/// combine = t0·t1 + t2 (per-variable degree 2).
fn combine(v: &[Ext2]) -> Ext2 {
    v[0] * v[1] + v[2]
}

fn setup(ell: usize, r: &mut ChaCha8Rng) -> (Vec<Vec<Ext2>>, Ext2) {
    let n = 1usize << ell;
    let tables: Vec<Vec<Ext2>> = (0..3)
        .map(|_| (0..n).map(|_| rand_ext(r)).collect())
        .collect();
    let mut claimed = Ext2::ZERO;
    for x in 0..n {
        claimed += combine(&[tables[0][x], tables[1][x], tables[2][x]]);
    }
    (tables, claimed)
}

#[test]
fn completeness_and_reduction() {
    let mut r = rng();
    for ell in 1..=8 {
        let (tables, claimed) = setup(ell, &mut r);
        let degree = 2;

        let mut tp = Transcript::new(b"sc-test");
        let (proof, rp) = sumcheck_prove(tables.clone(), combine, degree, ell, &mut tp);
        assert_eq!(proof.rounds[0].sum_over_bit(), claimed, "round-0 sum != T");

        let mut tv = Transcript::new(b"sc-test");
        let (rp2, claim) = sumcheck_verify(&proof, claimed, degree, ell, &mut tv).unwrap();
        assert_eq!(rp, rp2, "prover/verifier challenge points diverge");

        // The reduced claim must equal Q(r') = combine(MLE of each table at r').
        let q = combine(&[
            mle_eval_ext(&tables[0], &rp),
            mle_eval_ext(&tables[1], &rp),
            mle_eval_ext(&tables[2], &rp),
        ]);
        assert_eq!(claim, q, "reduced claim != Q(r')");
    }
}

#[test]
fn rejects_wrong_claimed_sum() {
    let mut r = rng();
    let ell = 5;
    let (tables, claimed) = setup(ell, &mut r);
    let degree = 2;
    let mut tp = Transcript::new(b"sc-test");
    let (proof, _) = sumcheck_prove(tables, combine, degree, ell, &mut tp);

    let mut tv = Transcript::new(b"sc-test");
    let bad = claimed + Ext2::ONE;
    assert!(sumcheck_verify(&proof, bad, degree, ell, &mut tv).is_err());
}

#[test]
fn rejects_tampered_round_poly() {
    let mut r = rng();
    let ell = 6;
    let (tables, claimed) = setup(ell, &mut r);
    let degree = 2;
    let mut tp = Transcript::new(b"sc-test");
    let (mut proof, _) = sumcheck_prove(tables, combine, degree, ell, &mut tp);

    // Corrupt one evaluation in round 2.
    proof.rounds[2].0[0] += Ext2::ONE;
    let mut tv = Transcript::new(b"sc-test");
    assert!(sumcheck_verify(&proof, claimed, degree, ell, &mut tv).is_err());
}

#[test]
fn transcript_is_deterministic() {
    let mut a = Transcript::new(b"dom");
    let mut b = Transcript::new(b"dom");
    a.absorb_fp(b"x", Fp::new(42));
    b.absorb_fp(b"x", Fp::new(42));
    assert_eq!(a.challenge_ext(b"c"), b.challenge_ext(b"c"));
    // Divergent absorb ⇒ divergent challenge.
    let mut c = Transcript::new(b"dom");
    c.absorb_fp(b"x", Fp::new(43));
    assert_ne!(a.challenge_ext(b"c"), c.challenge_ext(b"c"));
}

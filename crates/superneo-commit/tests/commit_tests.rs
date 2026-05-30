//! Ajtai commitment: module-homomorphism, parameter guard (Appendix D.8), challenge
//! set bounds, and seeded determinism.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::challenge::{sample_challenge, sample_challenges};
use superneo_commit::{CommitError, Params, PublicParams};
use superneo_field::fp::Fp;
use superneo_ring::{RingElem, D};

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0xCAFE_F00D_1234_5678)
}

fn rand_ring(r: &mut ChaCha8Rng) -> RingElem {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::new(r.next_u64());
    }
    RingElem::from_coeffs(c)
}

fn rand_witness(r: &mut ChaCha8Rng, n: usize) -> Vec<RingElem> {
    (0..n).map(|_| rand_ring(r)).collect()
}

#[test]
fn additive_homomorphism() {
    let mut r = rng();
    let n = 4;
    let pp = PublicParams::setup_seeded([7u8; 32], 18, n);
    for _ in 0..500 {
        let z1 = rand_witness(&mut r, n);
        let z2 = rand_witness(&mut r, n);
        let zsum: Vec<RingElem> = z1.iter().zip(&z2).map(|(a, b)| *a + *b).collect();
        let lhs = pp.commit(&z1).unwrap().add(&pp.commit(&z2).unwrap());
        let rhs = pp.commit(&zsum).unwrap();
        assert_eq!(lhs, rhs);
    }
}

#[test]
fn ring_scalar_homomorphism() {
    // Commit(ρ·z) = ρ·Commit(z) for ρ ∈ R_F.
    let mut r = rng();
    let n = 4;
    let pp = PublicParams::setup_seeded([9u8; 32], 18, n);
    for _ in 0..500 {
        let z = rand_witness(&mut r, n);
        let rho = rand_ring(&mut r);
        let scaled: Vec<RingElem> = z.iter().map(|zj| *zj * rho).collect();
        assert_eq!(
            pp.commit(&scaled).unwrap(),
            pp.commit(&z).unwrap().scale_ring(&rho)
        );
    }
}

#[test]
fn field_scalar_homomorphism() {
    // Commit(s·z) = s·Commit(z) for s ∈ F.
    let mut r = rng();
    let n = 5;
    let pp = PublicParams::setup_seeded([3u8; 32], 18, n);
    for _ in 0..500 {
        let z = rand_witness(&mut r, n);
        let s = Fp::new(r.next_u64());
        let scaled: Vec<RingElem> = z.iter().map(|zj| zj.scale(s)).collect();
        assert_eq!(
            pp.commit(&scaled).unwrap(),
            pp.commit(&z).unwrap().scale_fp(s)
        );
    }
}

#[test]
fn commit_rejects_wrong_length() {
    let pp = PublicParams::setup_seeded([1u8; 32], 4, 3);
    assert!(matches!(
        pp.commit(&[RingElem::ZERO; 2]),
        Err(CommitError::DimMismatch(_))
    ));
}

#[test]
fn seeded_setup_is_deterministic() {
    let a = PublicParams::setup_seeded([42u8; 32], 18, 6);
    let b = PublicParams::setup_seeded([42u8; 32], 18, 6);
    assert_eq!(a, b);
    let c = PublicParams::setup_seeded([43u8; 32], 18, 6);
    assert_ne!(a, c);
}

#[test]
fn goldilocks_params_validate_paper_d8() {
    let p = Params::GOLDILOCKS_B2;
    // (K+k)·T·(b−1) = (61+14)·216·1 = 16200 < 16384 = B   [Appendix D.8 prints True].
    let lhs = (p.max_fresh + p.k as u64) * p.t_exp * (p.b - 1);
    assert_eq!(lhs, 16200);
    assert_eq!(p.big_b, 16384);
    assert!(lhs < p.big_b);
    assert!(p.validate().is_ok());
    // Expansion factor T = 2·φ(81)·max_coeff = 2·54·2 = 216 (Theorem 11).
    assert_eq!(p.t_exp, 2 * (p.d as u64) * 2);
    assert_eq!(p.big_b, p.b.pow(p.k as u32));
}

#[test]
fn tampered_params_fail_guard() {
    // Increasing K past the guard threshold must fail validation.
    let mut p = Params::GOLDILOCKS_B2;
    p.max_fresh = 100; // (100+14)·216 = 24624 ≥ 16384
    assert!(matches!(p.validate(), Err(CommitError::GuardFailed { .. })));
    // Breaking B = b^k must also fail.
    let mut q = Params::GOLDILOCKS_B2;
    q.big_b = 12345;
    assert!(matches!(q.validate(), Err(CommitError::DimMismatch(_))));
}

#[test]
fn from_params_binds_kappa_and_validates() {
    // κ is derived from the validated param set, not a free argument: a caller cannot
    // silently downgrade the Module-SIS level.
    let pp = PublicParams::from_params([1u8; 32], &Params::GOLDILOCKS_B2, 6).unwrap();
    assert_eq!(pp.kappa, Params::GOLDILOCKS_B2.kappa);
    assert_eq!(pp.kappa, 18);

    // An invalid parameter set is rejected by the security-bearing constructor.
    let mut bad = Params::GOLDILOCKS_B2;
    bad.max_fresh = 10_000; // breaks (K+k)·T·(b−1) < B
    assert!(matches!(
        PublicParams::from_params([1u8; 32], &bad, 6),
        Err(CommitError::GuardFailed { .. })
    ));
}

#[test]
fn challenge_set_membership() {
    let mut r = rng();
    let p = Params::GOLDILOCKS_B2;
    for _ in 0..10_000 {
        let c = sample_challenge(&mut r, &p);
        for coeff in c.coeffs() {
            assert!(
                coeff.abs_balanced() <= 2,
                "challenge coeff outside {{-2..2}}"
            );
        }
    }
    assert_eq!(sample_challenges(&mut r, &p, 7).len(), 7);
}

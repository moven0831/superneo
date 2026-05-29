//! Goldilocks `Fp` correctness: field axioms, Solinas known vectors, Fermat, and a
//! reduce128-vs-reference cross-check (Definition 1, Appendix B.2).

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::{fp::Fp, Q};

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x5142_3370_DEAD_BEEF)
}

fn rand_fp(r: &mut ChaCha8Rng) -> Fp {
    Fp::new(r.next_u64())
}

#[test]
fn solinas_known_vectors() {
    // (q − 1) ≡ −1, so (q − 1)² ≡ 1.
    assert_eq!(Fp::new(Q - 1) * Fp::new(Q - 1), Fp::ONE);
    // 2^64 ≡ 2^32 − 1 (mod q): (2^32)·(2^32) = 2^64 ≡ 2^32 − 1.
    let two32 = Fp::new(1 << 32);
    assert_eq!(two32 * two32, Fp::new((1 << 32) - 1));
    // 2^96 ≡ −1: (2^32)·(2^64) ≡ (2^32)·(2^32 − 1) ≡ −1.
    let two64 = two32 * two32; // = 2^32 − 1 (mod q)
    assert_eq!(two32 * two64, Fp::new(Q - 1));
    // 1⁻¹ = 1.
    assert_eq!(Fp::ONE.inv().unwrap(), Fp::ONE);
}

#[test]
fn reduce128_matches_reference() {
    let mut r = rng();
    for _ in 0..200_000 {
        let x = ((r.next_u64() as u128) << 64) | (r.next_u64() as u128);
        let got = Fp::reduce128(x).to_u64();
        let want = (x % (Q as u128)) as u64;
        assert_eq!(got, want, "reduce128 mismatch for x = {x:#034x}");
    }
}

#[test]
fn field_axioms() {
    let mut r = rng();
    for _ in 0..50_000 {
        let (a, b, c) = (rand_fp(&mut r), rand_fp(&mut r), rand_fp(&mut r));
        // commutativity
        assert_eq!(a + b, b + a);
        assert_eq!(a * b, b * a);
        // associativity
        assert_eq!((a + b) + c, a + (b + c));
        assert_eq!((a * b) * c, a * (b * c));
        // distributivity
        assert_eq!(a * (b + c), a * b + a * c);
        // identities & inverses
        assert_eq!(a + Fp::ZERO, a);
        assert_eq!(a * Fp::ONE, a);
        assert_eq!(a + (-a), Fp::ZERO);
        assert_eq!(a - b, a + (-b));
    }
}

#[test]
fn inverse_roundtrip_and_fermat() {
    let mut r = rng();
    for _ in 0..20_000 {
        let a = rand_fp(&mut r);
        if a.is_zero() {
            assert!(a.inv().is_err());
            continue;
        }
        // a · a⁻¹ = 1
        assert_eq!(a * a.inv().unwrap(), Fp::ONE);
        // Fermat: a^{q−1} = 1
        assert_eq!(a.pow(Q - 1), Fp::ONE);
    }
}

#[test]
fn balanced_representative() {
    assert_eq!(Fp::ZERO.balanced(), 0);
    assert_eq!(Fp::ONE.balanced(), 1);
    assert_eq!(Fp::new(Q - 1).balanced(), -1);
    assert_eq!(Fp::new(Q - 1).abs_balanced(), 1);
    // Just above the midpoint is negative; just below is positive.
    assert!(Fp::new(Q / 2 + 1).balanced() < 0);
    assert!(Fp::new(Q / 2).balanced() > 0);
    // |balanced| < q/2 always.
    let mut r = rng();
    for _ in 0..10_000 {
        assert!(rand_fp(&mut r).abs_balanced() <= Q / 2);
    }
}

#[test]
fn pow_small_consistency() {
    let mut r = rng();
    for _ in 0..5_000 {
        let a = rand_fp(&mut r);
        assert_eq!(a.pow(0), Fp::ONE);
        assert_eq!(a.pow(1), a);
        assert_eq!(a.pow(2), a.square());
        assert_eq!(a.pow(3), a * a * a);
    }
}

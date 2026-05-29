//! Ring `R_F = F[X]/(X^54+X^27+1)`: reduction correctness, arithmetic, and the
//! rotation-matrix identity `rot(a)·cf(b) = cf(a·b)`.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::fp::Fp;
use superneo_ring::ring::{RingElem, D};
use superneo_ring::s_action::{matvec, rot};

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x9E37_79B9_7F4A_7C15)
}

fn rand_ring(r: &mut ChaCha8Rng) -> RingElem {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::new(r.next_u64());
    }
    RingElem::from_coeffs(c)
}

#[test]
fn x_pow_81_is_one() {
    // Φ_81 | X^81 − 1, so X^81 ≡ 1 (mod Φ).
    assert_eq!(RingElem::x_pow_mod_phi(81), RingElem::one());
}

#[test]
fn phi_annihilates() {
    // X^54 + X^27 + 1 ≡ 0 (mod Φ).
    let sum = RingElem::x_pow_mod_phi(54) + RingElem::x_pow_mod_phi(27) + RingElem::one();
    assert!(sum.is_zero(), "Φ does not annihilate: {sum:?}");
}

#[test]
fn x54_reduces_correctly() {
    // X^54 = −X^27 − 1.
    let x54 = RingElem::x_pow_mod_phi(54);
    let mut expect = [Fp::ZERO; D];
    expect[0] = -Fp::ONE;
    expect[27] = -Fp::ONE;
    assert_eq!(x54, RingElem::from_coeffs(expect));
}

#[test]
fn ring_axioms() {
    let mut r = rng();
    for _ in 0..20_000 {
        let (a, b, c) = (rand_ring(&mut r), rand_ring(&mut r), rand_ring(&mut r));
        assert_eq!(a + b, b + a);
        assert_eq!(a * b, b * a);
        assert_eq!((a + b) + c, a + (b + c));
        assert_eq!((a * b) * c, a * (b * c));
        assert_eq!(a * (b + c), a * b + a * c);
        assert_eq!(a * RingElem::one(), a);
        assert_eq!(a + RingElem::ZERO, a);
    }
}

#[test]
fn multiplication_by_x_matches_shift() {
    let mut r = rng();
    for _ in 0..5_000 {
        let a = rand_ring(&mut r);
        assert_eq!(a * RingElem::x(), a.mul_x_pow(1));
        assert_eq!(a * RingElem::x_pow_mod_phi(5), a.mul_x_pow(5));
    }
}

#[test]
fn rotation_matrix_realizes_multiplication() {
    // cf(a·b) = rot(a)·cf(b) for all a, b.
    let mut r = rng();
    for _ in 0..5_000 {
        let (a, b) = (rand_ring(&mut r), rand_ring(&mut r));
        let lhs = (a * b).into_coeffs();
        let rhs = matvec(&rot(&a), b.coeffs());
        assert_eq!(lhs, rhs);
    }
}

#[test]
fn scale_matches_constant_multiplication() {
    let mut r = rng();
    for _ in 0..5_000 {
        let a = rand_ring(&mut r);
        let s = Fp::new(r.next_u64());
        let s_ring = RingElem::from_coeffs({
            let mut c = [Fp::ZERO; D];
            c[0] = s;
            c
        });
        assert_eq!(a.scale(s), a * s_ring);
    }
}

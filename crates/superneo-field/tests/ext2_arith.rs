//! `Ext2 = F_q[u]/(u² − 7)` correctness: non-residue check, field axioms, inverse,
//! and the Frobenius identity `z^q = conj(z)` (Definition 1).

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::{ext2::Ext2, fp::Fp, Q};

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x00C0_FFEE_1234_5678)
}

fn rand_fp(r: &mut ChaCha8Rng) -> Fp {
    Fp::new(r.next_u64())
}

fn rand_ext2(r: &mut ChaCha8Rng) -> Ext2 {
    Ext2::new(rand_fp(r), rand_fp(r))
}

#[test]
fn seven_is_a_non_residue() {
    // 7 is a QNR iff 7^{(q−1)/2} = −1.
    let legendre = Fp::new(7).pow((Q - 1) / 2);
    assert_eq!(legendre, Fp::new(Q - 1), "7 must be a quadratic non-residue");
}

#[test]
fn u_squared_is_seven() {
    let u = Ext2::new(Fp::ZERO, Fp::ONE);
    assert_eq!(u.square(), Ext2::from_base(Fp::new(7)));
}

#[test]
fn field_axioms() {
    let mut r = rng();
    for _ in 0..50_000 {
        let (a, b, c) = (rand_ext2(&mut r), rand_ext2(&mut r), rand_ext2(&mut r));
        assert_eq!(a + b, b + a);
        assert_eq!(a * b, b * a);
        assert_eq!((a + b) + c, a + (b + c));
        assert_eq!((a * b) * c, a * (b * c));
        assert_eq!(a * (b + c), a * b + a * c);
        assert_eq!(a + Ext2::ZERO, a);
        assert_eq!(a * Ext2::ONE, a);
        assert_eq!(a + (-a), Ext2::ZERO);
    }
}

#[test]
fn base_field_embedding_is_a_homomorphism() {
    let mut r = rng();
    for _ in 0..20_000 {
        let (x, y) = (rand_fp(&mut r), rand_fp(&mut r));
        assert_eq!(Ext2::from_base(x) + Ext2::from_base(y), Ext2::from_base(x + y));
        assert_eq!(Ext2::from_base(x) * Ext2::from_base(y), Ext2::from_base(x * y));
        // mul_base agrees with full multiplication by an embedded scalar.
        let z = rand_ext2(&mut r);
        assert_eq!(z.mul_base(x), z * Ext2::from_base(x));
    }
}

#[test]
fn inverse_roundtrip() {
    let mut r = rng();
    for _ in 0..20_000 {
        let z = rand_ext2(&mut r);
        if z.is_zero() {
            assert!(z.inv().is_err());
            continue;
        }
        assert_eq!(z * z.inv().unwrap(), Ext2::ONE);
    }
}

#[test]
fn frobenius_is_conjugation() {
    // The nontrivial automorphism of F_{q^2} is x ↦ x^q, and it equals conjugation.
    let mut r = rng();
    for _ in 0..2_000 {
        let z = rand_ext2(&mut r);
        assert_eq!(z.pow(Q), z.conjugate());
    }
}

#[test]
fn norm_lands_in_base_field_and_is_multiplicative() {
    let mut r = rng();
    for _ in 0..20_000 {
        let (x, y) = (rand_ext2(&mut r), rand_ext2(&mut r));
        // N(xy) = N(x)·N(y)
        assert_eq!((x * y).norm(), x.norm() * y.norm());
        // N(z) = z·conj(z) (constant term)
        assert_eq!(Ext2::from_base(x.norm()), x * x.conjugate());
    }
}

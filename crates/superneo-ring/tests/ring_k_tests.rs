//! `R_K = K[X]/(Φ)`: embedding `R_F ↪ R_K` is a ring homomorphism and the mixed
//! product `R_K × R_F → R_F` agrees with `R_F` multiplication.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_ring::ring::{RingElem, D};
use superneo_ring::ring_k::RingK;

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x0BAD_F00D_5EED_1234)
}

fn rand_rf(r: &mut ChaCha8Rng) -> RingElem {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::new(r.next_u64());
    }
    RingElem::from_coeffs(c)
}

#[test]
fn embedding_is_ring_homomorphism() {
    // from_rf(a)·from_rf(b) = from_rf(a·b), realized through mul_rf.
    let mut r = rng();
    for _ in 0..5_000 {
        let a = rand_rf(&mut r);
        let b = rand_rf(&mut r);
        // from_rf(a) as RingK, times b ∈ R_F:
        let lhs = RingK::from_rf(&a).mul_rf(&b);
        let rhs = RingK::from_rf(&(a * b));
        assert_eq!(lhs, rhs);
    }
}

#[test]
fn add_and_scale_are_consistent() {
    let mut r = rng();
    for _ in 0..5_000 {
        let a = rand_rf(&mut r);
        let b = rand_rf(&mut r);
        let ka = RingK::from_rf(&a);
        let kb = RingK::from_rf(&b);
        // addition matches R_F addition under the embedding
        assert_eq!(ka.add(&kb), RingK::from_rf(&(a + b)));
        // field scaling matches
        let s = Fp::new(r.next_u64());
        assert_eq!(ka.scale_fp(s), RingK::from_rf(&a.scale(s)));
    }
}

#[test]
fn ct_picks_constant_coefficient() {
    let mut r = rng();
    for _ in 0..2_000 {
        let a = rand_rf(&mut r);
        assert_eq!(RingK::from_rf(&a).ct(), Ext2::from_base(a.coeffs()[0]));
    }
}

#[test]
fn ext_scaling_distributes_over_constant() {
    // (k · 1_RK) scaled is just multiplication of constants.
    let mut r = rng();
    for _ in 0..2_000 {
        let k = Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64()));
        let m = Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64()));
        assert_eq!(RingK::from_const(k).scale_ext(m).ct(), k * m);
    }
}

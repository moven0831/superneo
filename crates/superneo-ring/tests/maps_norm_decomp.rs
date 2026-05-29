//! Embedding (Definition 7), the balanced ℓ∞ norm (Definition 3), and the base-`b`
//! decomposition `split_b` (round-trip + norm bound).

use rand::{Rng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::fp::Fp;
use superneo_ring::decomp::{recompose_b, split_b};
use superneo_ring::maps::embed_vector;
use superneo_ring::norm::{norm_inf_fp, norm_inf_ring_vec};
use superneo_ring::ring::D;
use superneo_ring::RingError;

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x1357_9BDF_2468_ACE0)
}

#[test]
fn embed_chunks_into_d_blocks() {
    let mut r = rng();
    let n_r = 5usize;
    let z: Vec<Fp> = (0..D * n_r).map(|_| Fp::new(r.next_u64())).collect();
    let ring = embed_vector(&z).unwrap();
    assert_eq!(ring.len(), n_r);
    // block i holds coordinates [i·D, (i+1)·D)
    for (i, elem) in ring.iter().enumerate() {
        for j in 0..D {
            assert_eq!(elem.coeffs()[j], z[i * D + j]);
        }
    }
}

#[test]
fn embed_rejects_misaligned() {
    let z = vec![Fp::ONE; D + 1];
    assert!(matches!(embed_vector(&z), Err(RingError::Misaligned(_))));
}

#[test]
fn norm_inf_is_max_balanced_magnitude() {
    let v = vec![Fp::ONE, -Fp::new(5), Fp::new(3), Fp::ZERO];
    assert_eq!(norm_inf_fp(&v), 5);
    let ring = embed_vector(&{
        let mut z = vec![Fp::ZERO; D * 2];
        z[0] = Fp::new(7);
        z[D] = -Fp::new(9);
        z
    })
    .unwrap();
    assert_eq!(norm_inf_ring_vec(&ring), 9);
}

#[test]
fn split_b_roundtrip_and_norm_bound() {
    // b = 2, k = 14 ⇒ bound = 2^14 = 16384, matching the Goldilocks parameter set.
    let b = 2u64;
    let k = 14usize;
    let bound = b.pow(k as u32) as i64;
    let mut r = rng();
    for _ in 0..2_000 {
        // coordinates with balanced magnitude < bound
        let z: Vec<Fp> = (0..32)
            .map(|_| {
                let v = r.random_range(-(bound - 1)..bound);
                Fp::from_i64(v)
            })
            .collect();
        let parts = split_b(&z, b, k).unwrap();
        assert_eq!(parts.len(), k);
        // each part is low-norm: ‖z_i‖∞ < b
        for part in &parts {
            assert!(
                norm_inf_fp(part) < b,
                "decomposition part exceeds norm bound"
            );
        }
        // exact reconstruction
        assert_eq!(recompose_b(&parts, b), z);
    }
}

#[test]
fn split_b_rejects_out_of_range() {
    // A coordinate at exactly b^k cannot be represented in k digits.
    let b = 2u64;
    let k = 4usize; // bound = 16
    let z = vec![Fp::from_i64(16)];
    assert!(matches!(
        split_b(&z, b, k),
        Err(RingError::NormBound { .. })
    ));
    // ... but b^k − 1 can.
    let z_ok = vec![Fp::from_i64(15)];
    assert!(split_b(&z_ok, b, k).is_ok());
}

#[test]
fn split_b_general_base() {
    // Sanity for a non-power-of-two base.
    let b = 5u64;
    let k = 6usize;
    let bound = b.pow(k as u32) as i64;
    let mut r = rng();
    for _ in 0..2_000 {
        let z: Vec<Fp> = (0..16)
            .map(|_| Fp::from_i64(r.random_range(-(bound - 1)..bound)))
            .collect();
        let parts = split_b(&z, b, k).unwrap();
        for part in &parts {
            assert!(norm_inf_fp(part) < b);
        }
        assert_eq!(recompose_b(&parts, b), z);
    }
}

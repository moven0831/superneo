//! BaseFold PCS: an honest evaluation proof round-trips, the proved value matches the
//! true multilinear evaluation, and tampering with the codeword, a layer root, the
//! claimed value, or a sum-check message is rejected.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::Transcript;
use superneo_snark::basefold::{commit, open, verify};

const DOMAIN: &[u8] = b"basefold/test";

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0xBA5E_F01D_0000_2026)
}

fn rand_point(r: &mut ChaCha8Rng, nu: usize) -> Vec<Ext2> {
    (0..nu)
        .map(|_| Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64())))
        .collect()
}

/// Direct little-endian multilinear evaluation, for cross-checking the proved value.
fn mle_eval_le(evals: &[Fp], point: &[Ext2]) -> Ext2 {
    let mut cur: Vec<Ext2> = evals.iter().map(|&x| Ext2::from_base(x)).collect();
    for &rk in point {
        cur = cur
            .chunks_exact(2)
            .map(|p| p[0] + rk * (p[1] - p[0]))
            .collect();
    }
    cur[0]
}

fn sample(nu: usize, r: &mut ChaCha8Rng) -> Vec<Fp> {
    (0..(1 << nu)).map(|_| Fp::new(r.next_u64())).collect()
}

#[test]
fn open_verify_round_trip_all_sizes() {
    let mut r = rng();
    for nu in 1..=8 {
        let evals = sample(nu, &mut r);
        let point = rand_point(&mut r, nu);
        let (comm, data) = commit(&evals);
        assert_eq!(comm.num_vars, nu);

        let mut tp = Transcript::new(DOMAIN);
        let proof = open(&mut tp, &comm, &data, &point);
        assert_eq!(
            proof.value,
            mle_eval_le(&evals, &point),
            "proved value ≠ true f(r) at nu={nu}"
        );

        let mut tv = Transcript::new(DOMAIN);
        verify(&mut tv, &comm, &point, &proof).unwrap_or_else(|e| panic!("nu={nu}: {e}"));
    }
}

#[test]
fn rejects_wrong_claimed_value() {
    let mut r = rng();
    let nu = 6;
    let evals = sample(nu, &mut r);
    let point = rand_point(&mut r, nu);
    let (comm, data) = commit(&evals);
    let mut tp = Transcript::new(DOMAIN);
    let mut proof = open(&mut tp, &comm, &data, &point);

    proof.value += Ext2::ONE;
    let mut tv = Transcript::new(DOMAIN);
    assert!(verify(&mut tv, &comm, &point, &proof).is_err());
}

#[test]
fn rejects_tampered_sumcheck() {
    let mut r = rng();
    let nu = 6;
    let evals = sample(nu, &mut r);
    let point = rand_point(&mut r, nu);
    let (comm, data) = commit(&evals);
    let mut tp = Transcript::new(DOMAIN);
    let mut proof = open(&mut tp, &comm, &data, &point);

    proof.round_polys[2][0] += Ext2::ONE;
    let mut tv = Transcript::new(DOMAIN);
    assert!(verify(&mut tv, &comm, &point, &proof).is_err());
}

#[test]
fn rejects_tampered_final_constant() {
    let mut r = rng();
    let nu = 5;
    let evals = sample(nu, &mut r);
    let point = rand_point(&mut r, nu);
    let (comm, data) = commit(&evals);
    let mut tp = Transcript::new(DOMAIN);
    let mut proof = open(&mut tp, &comm, &data, &point);

    proof.final_constant += Ext2::ONE;
    let mut tv = Transcript::new(DOMAIN);
    assert!(verify(&mut tv, &comm, &point, &proof).is_err());
}

#[test]
fn rejects_tampered_query_leaf() {
    let mut r = rng();
    let nu = 6;
    let evals = sample(nu, &mut r);
    let point = rand_point(&mut r, nu);
    let (comm, data) = commit(&evals);
    let mut tp = Transcript::new(DOMAIN);
    let mut proof = open(&mut tp, &comm, &data, &point);

    // Corrupt an opened base-layer leaf: its Merkle path no longer matches the root.
    proof.queries[0].layers[0].lo += Ext2::ONE;
    let mut tv = Transcript::new(DOMAIN);
    assert!(verify(&mut tv, &comm, &point, &proof).is_err());
}

#[test]
fn rejects_tampered_commitment_root() {
    let mut r = rng();
    let nu = 5;
    let evals = sample(nu, &mut r);
    let point = rand_point(&mut r, nu);
    let (mut comm, data) = commit(&evals);
    let mut tp = Transcript::new(DOMAIN);
    // Prover commits honestly...
    let proof = open(&mut tp, &comm, &data, &point);
    // ...but the verifier is handed a different root.
    comm.root[0] ^= 0xFF;
    let mut tv = Transcript::new(DOMAIN);
    assert!(verify(&mut tv, &comm, &point, &proof).is_err());
}

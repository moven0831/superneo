//! Bench for the Ajtai commitment `Commit(A, z) = A·z` (M8).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::fp::Fp;
use superneo_ring::{RingElem, D};

const KAPPA: usize = 18;
const N_COLS: usize = 64;

fn rand_ring(r: &mut ChaCha8Rng) -> RingElem {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::new(r.next_u64());
    }
    RingElem::from_coeffs(c)
}

fn benches(c: &mut Criterion) {
    let pp = PublicParams::setup_seeded([3u8; 32], KAPPA, N_COLS);
    let mut r = ChaCha8Rng::seed_from_u64(0xA171);
    let z: Vec<RingElem> = (0..N_COLS).map(|_| rand_ring(&mut r)).collect();

    c.bench_function("ajtai_commit", |bn| {
        bn.iter(|| pp.commit(black_box(&z)).unwrap())
    });
}

criterion_group!(g, benches);
criterion_main!(g);

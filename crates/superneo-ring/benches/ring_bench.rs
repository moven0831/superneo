//! Benches for the ring layer: cyclotomic multiplication and the bar transform (M8).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::fp::Fp;
use superneo_ring::bar::{bar_block, bar_ring};
use superneo_ring::{RingElem, D};

fn rand_block(r: &mut ChaCha8Rng) -> [Fp; D] {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::new(r.next_u64());
    }
    c
}

fn benches(c: &mut Criterion) {
    let mut r = ChaCha8Rng::seed_from_u64(0x1234_5678);
    let a = RingElem::from_coeffs(rand_block(&mut r));
    let b = RingElem::from_coeffs(rand_block(&mut r));
    let blk = rand_block(&mut r);

    c.bench_function("ring_mul", |bn| bn.iter(|| black_box(a) * black_box(b)));
    c.bench_function("bar_block", |bn| bn.iter(|| bar_block(black_box(&blk))));
    c.bench_function("bar_ring", |bn| bn.iter(|| bar_ring(black_box(&blk))));
}

criterion_group!(g, benches);
criterion_main!(g);

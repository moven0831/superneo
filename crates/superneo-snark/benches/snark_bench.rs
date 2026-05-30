//! Benches for the final compression: `compress` and `verify` of an accumulator (M8).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::types::compute_evals;
use superneo_fold::{CcsStructure, CeInstance, CeWitness, GlobalParams};
use superneo_ring::{RingElem, D};
use superneo_snark::{compress, verify};

const KAPPA: usize = 18;
const N_R: usize = 1;

fn low_norm_ring(r: &mut ChaCha8Rng) -> Vec<RingElem> {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::from_i64((r.next_u64() % 3) as i64 - 1);
    }
    vec![RingElem::from_coeffs(c)]
}

fn benches(c: &mut Criterion) {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([5u8; 32], KAPPA, N_R);
    let s = CcsStructure::identity(gp.m);
    let mut r = ChaCha8Rng::seed_from_u64(0x5_4A6E);

    let point: Vec<Ext2> = (0..gp.log_m)
        .map(|_| Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64())))
        .collect();
    let (mut acc, mut acc_wit) = (Vec::new(), Vec::new());
    for _ in 0..gp.k {
        let zr = low_norm_ring(&mut r);
        acc.push(CeInstance {
            c: pp.commit(&zr).unwrap(),
            r: point.clone(),
            y: compute_evals(&s, &zr, &point),
        });
        acc_wit.push(CeWitness { z_ring: zr });
    }

    c.bench_function("compress", |bn| {
        bn.iter(|| compress(&gp, &pp, &s, black_box(&acc), black_box(&acc_wit)).unwrap())
    });

    let proof = compress(&gp, &pp, &s, &acc, &acc_wit).unwrap();
    c.bench_function("verify", |bn| {
        bn.iter(|| verify(&gp, &pp, &s, black_box(&acc), black_box(&proof)).unwrap())
    });
}

criterion_group!(g, benches);
criterion_main!(g);

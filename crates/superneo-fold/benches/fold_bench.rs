//! Bench for one full folding step `Π_DEC ∘ Π_RLC ∘ Π_CCS` (M8).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::types::{commit_witness, compute_evals};
use superneo_fold::{
    fold, CcsInstance, CcsStructure, CcsWitness, CeInstance, CeWitness, GlobalParams, Transcript,
};
use superneo_ring::{RingElem, D};

const KAPPA: usize = 18;
const N_R: usize = 1;
const CAP_K: usize = 2;

fn low_norm_field(r: &mut ChaCha8Rng, n: usize) -> Vec<Fp> {
    (0..n)
        .map(|_| Fp::from_i64((r.next_u64() % 3) as i64 - 1))
        .collect()
}

fn low_norm_ring(r: &mut ChaCha8Rng) -> Vec<RingElem> {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::from_i64((r.next_u64() % 3) as i64 - 1);
    }
    vec![RingElem::from_coeffs(c)]
}

fn benches(c: &mut Criterion) {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([7u8; 32], KAPPA, N_R);
    let s = CcsStructure::identity(gp.m);
    let mut r = ChaCha8Rng::seed_from_u64(0xF01D_5713);

    let point: Vec<Ext2> = (0..gp.log_m)
        .map(|_| Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64())))
        .collect();
    let (mut carried, mut carried_wit) = (Vec::new(), Vec::new());
    for _ in 0..gp.k {
        let zr = low_norm_ring(&mut r);
        carried.push(CeInstance {
            c: pp.commit(&zr).unwrap(),
            r: point.clone(),
            y: compute_evals(&s, &zr, &point),
        });
        carried_wit.push(CeWitness { z_ring: zr });
    }
    let (mut fresh, mut fresh_wit) = (Vec::new(), Vec::new());
    for _ in 0..CAP_K {
        let z = low_norm_field(&mut r, gp.n_f);
        fresh.push(CcsInstance {
            c: commit_witness(&pp, &z),
        });
        fresh_wit.push(CcsWitness { z });
    }

    c.bench_function("fold_step", |bn| {
        bn.iter(|| {
            let mut tr = Transcript::new(b"superneo/bench/fold");
            fold(
                &mut tr,
                &gp,
                &pp,
                &s,
                black_box(&fresh),
                black_box(&fresh_wit),
                black_box(&carried),
                black_box(&carried_wit),
            )
            .unwrap()
        })
    });
}

criterion_group!(g, benches);
criterion_main!(g);

//! The full pipeline (M8 preview): fold a multi-step IVC run, then compress the final
//! accumulator it produces into a succinct proof that verifies.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::fp::Fp;
use superneo_fold::{CcsStructure, CcsWitness, GlobalParams};
use superneo_ivc::prove_ivc_with_witnesses;
use superneo_snark::{compress, verify};

const KAPPA: usize = 18;
const N_R: usize = 1;
const STEPS: usize = 4;

fn low_norm_field(r: &mut ChaCha8Rng, n: usize) -> Vec<Fp> {
    (0..n)
        .map(|_| Fp::from_i64((r.next_u64() % 3) as i64 - 1))
        .collect()
}

#[test]
fn ivc_run_then_compressed_proof_verifies() {
    let gp = GlobalParams::goldilocks(N_R);
    let pp = PublicParams::setup_seeded([7u8; 32], KAPPA, N_R);
    let s = CcsStructure::identity(gp.m);
    let mut r = ChaCha8Rng::seed_from_u64(0x0E2E_0E2E_0E2E_0E2E);

    let steps: Vec<Vec<CcsWitness>> = (0..STEPS)
        .map(|_| {
            vec![CcsWitness {
                z: low_norm_field(&mut r, gp.n_f),
            }]
        })
        .collect();

    // Prove the IVC run and keep the final accumulator's (secret) witnesses.
    let (ivc_proof, final_wit) = prove_ivc_with_witnesses(&gp, &pp, &s, &steps).expect("prove ivc");
    assert_eq!(ivc_proof.final_acc.len(), gp.k);
    assert_eq!(final_wit.len(), gp.k);

    // Compress the accumulator the IVC produced and verify the succinct proof.
    let snark = compress(&gp, &pp, &s, &ivc_proof.final_acc, &final_wit).expect("compress");
    verify(&gp, &pp, &s, &ivc_proof.final_acc, &snark).expect("verify compressed proof");
}

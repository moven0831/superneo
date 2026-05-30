//! The recursive verifier circuit (M7): a sum-check verifier synthesized as CCS
//! constraints accepts an honest transcript, rejects a tampered one, and finalizes into
//! a fold-compatible CCS instance whose witness satisfies the relation (loop closure).

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_fold::sumcheck::{sumcheck_prove, sumcheck_verify};
use superneo_fold::Transcript;
use superneo_ivc::circuit::build_sumcheck_verifier_circuit;

const ELL: usize = 6; // log_m at the PoC scale
const DEGREE: usize = 4; // a four-table product — the fold's round degree
const N_TABLES: usize = 4;
const SC_DOMAIN: &[u8] = b"recursive/sumcheck";

fn rand_ext(r: &mut ChaCha8Rng) -> Ext2 {
    Ext2::new(Fp::new(r.next_u64()), Fp::new(r.next_u64()))
}

/// Produce a degree-`DEGREE` sum-check transcript: round polynomials, challenges, the
/// claimed sum `T`, and the verifier's final reduced claim.
fn sumcheck_trace() -> (Vec<Vec<Ext2>>, Vec<Ext2>, Ext2, Ext2) {
    let mut r = ChaCha8Rng::seed_from_u64(0x5EC0_0DED_1234_5678);
    let n = 1usize << ELL;
    let tables: Vec<Vec<Ext2>> = (0..N_TABLES)
        .map(|_| (0..n).map(|_| rand_ext(&mut r)).collect())
        .collect();
    let combine = |v: &[Ext2]| v.iter().copied().fold(Ext2::ONE, |a, b| a * b);

    let mut tp = Transcript::new(SC_DOMAIN);
    let (proof, _r_prime) = sumcheck_prove(tables, combine, DEGREE, ELL, &mut tp);

    let init_claim = proof.rounds[0].sum_over_bit();
    let mut tv = Transcript::new(SC_DOMAIN);
    let (challenges, final_claim) =
        sumcheck_verify(&proof, init_claim, DEGREE, ELL, &mut tv).expect("honest verify");

    let round_evals: Vec<Vec<Ext2>> = proof.rounds.iter().map(|p| p.0.clone()).collect();
    (round_evals, challenges, init_claim, final_claim)
}

#[test]
fn recursive_sumcheck_circuit_is_satisfied() {
    let (rounds, challenges, init, final_claim) = sumcheck_trace();
    let cs = build_sumcheck_verifier_circuit(init, &rounds, &challenges, final_claim);
    assert!(
        cs.is_satisfied(),
        "honest transcript must satisfy the circuit"
    );

    // Report the recursion cost.
    println!(
        "recursive sum-check verifier: {} vars, {} constraints",
        cs.num_vars(),
        cs.num_constraints()
    );
    assert!(
        (50..20_000).contains(&cs.num_constraints()),
        "constraint count {} outside the expected PoC band",
        cs.num_constraints()
    );
}

#[test]
fn finalizes_to_satisfying_ccs_instance() {
    let (rounds, challenges, init, final_claim) = sumcheck_trace();
    let cs = build_sumcheck_verifier_circuit(init, &rounds, &challenges, final_claim);
    // Loop closure: the verifier circuit is a CCS instance of the shape `fold` consumes,
    // and its witness satisfies the CCS relation.
    let (structure, z) = cs.finalize();
    assert_eq!(structure.t(), 4);
    assert_eq!(structure.m % superneo_ring::D, 0);
    assert_eq!(z.len(), structure.m);
    assert!(cs.finalized_relation_holds(), "CCS relation must hold");
}

#[test]
fn rejects_tampered_round_polynomial() {
    let (mut rounds, challenges, init, final_claim) = sumcheck_trace();
    // Corrupt one evaluation of a middle round polynomial.
    rounds[2][0] += Ext2::ONE;
    let cs = build_sumcheck_verifier_circuit(init, &rounds, &challenges, final_claim);
    assert!(
        !cs.is_satisfied(),
        "a tampered round polynomial must break the circuit"
    );
}

#[test]
fn rejects_tampered_final_claim() {
    let (rounds, challenges, init, final_claim) = sumcheck_trace();
    // The circuit binds the reduced claim to the wrong final value.
    let cs = build_sumcheck_verifier_circuit(init, &rounds, &challenges, final_claim + Ext2::ONE);
    assert!(!cs.is_satisfied(), "a wrong final claim must be rejected");
}

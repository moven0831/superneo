//! End-to-end SuperNeo demo (M8).
//!
//! Runs the full pipeline — genesis → fold a multi-step IVC run → compress the final
//! accumulator → verify — and demonstrates that ordinary computations compile to the
//! SuperNeo CCS shape. Run with `cargo run -p superneo-snark --example demo --release`.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_commit::PublicParams;
use superneo_field::fp::Fp;
use superneo_fold::{CcsStructure, CcsWitness, GlobalParams};
use superneo_ivc::circuit::ConstraintSystem;
use superneo_snark::pipeline::{
    counter_circuit, fibonacci_circuit, mul_gate_witness, multiplication_structure,
    multiplier_circuit, run_pipeline, PipelineReport,
};

const KAPPA: usize = 18;
const N_R: usize = 1;
const STEPS: usize = 4;

fn low_norm_field(r: &mut ChaCha8Rng, n: usize) -> Vec<Fp> {
    (0..n)
        .map(|_| Fp::from_i64((r.next_u64() % 3) as i64 - 1))
        .collect()
}

fn hex8(d: &[u8; 32]) -> String {
    d[..8].iter().map(|b| format!("{b:02x}")).collect()
}

fn report_line(label: &str, rep: &PipelineReport) {
    println!(
        "  {label:<26}  fold {:>7.2?} | compress {:>7.2?} | verify {:>7.2?} | proof ~{:>5} B | digest {}… | {}",
        rep.fold_time,
        rep.compress_time,
        rep.verify_time,
        rep.proof_bytes,
        hex8(&rep.ivc_digest),
        if rep.verified { "VERIFIED ✓" } else { "FAILED ✗" },
    );
}

fn circuit_line(name: &str, cs: &ConstraintSystem) {
    let (structure, _z) = cs.finalize();
    // `is_satisfied` is an in-memory witness-satisfaction oracle (it re-evaluates the
    // recorded constraints against the prover's witness), not a verified proof.
    println!(
        "  {name:<26}  {:>3} constraints, {:>3} vars  →  CCS t={}, m={:<4}  | witness-satisfies: {}",
        cs.num_constraints(),
        cs.num_vars(),
        structure.t(),
        structure.m,
        if cs.is_satisfied() { "yes" } else { "no" },
    );
}

fn main() {
    let gp = GlobalParams::goldilocks(N_R);
    // Security-bearing setup: κ is derived from the validated parameter set (not a free arg).
    let pp =
        PublicParams::from_params([9u8; 32], &gp.params(), N_R).expect("valid Goldilocks params");
    let mut rng = ChaCha8Rng::seed_from_u64(0x5EED_04E0_9E37_79B9);

    println!("\n╔══════════════════════════════════════════════════════════════════════════╗");
    println!("║  SuperNeo — post-quantum lattice folding for CCS  (ePrint 2026/242, PoC)   ║");
    println!("╚══════════════════════════════════════════════════════════════════════════╝");
    println!(
        "\nParameters (Appendix B.2, Goldilocks): q = 2⁶⁴−2³²+1, Φ = X⁵⁴+X²⁷+1 (d={}),\n\
         κ = {}, b = {}, k = {} (B = b^k = {}), n_R = {}  ⇒  m = n_f = {}, log_m = {}.",
        gp.d, KAPPA, gp.b, gp.k, gp.cap_b, N_R, gp.m, gp.log_m
    );

    // 1. Recursion-overhead pipeline: identity structure (application cost ≈ 0), so the
    //    cost is purely the fold machinery (k carried CE claims + one sum-check per step).
    println!("\n── Full pipeline: fold → IVC → compress → verify ──────────────────────────\n");
    let s_id = CcsStructure::identity(gp.m);
    let id_steps: Vec<Vec<CcsWitness>> = (0..STEPS)
        .map(|_| {
            vec![CcsWitness {
                z: low_norm_field(&mut rng, gp.n_f),
            }]
        })
        .collect();
    let id_rep = run_pipeline(&gp, &pp, &s_id, &id_steps).expect("identity pipeline");
    report_line("identity step ×4", &id_rep);

    // 2. A real constraint folded end-to-end: a multiplication gate b·b = b, b ∈ {0,1}
    //    (low-norm, so it folds through the norm-bounded scheme).
    let s_mul = multiplication_structure(gp.m);
    let mul_steps: Vec<Vec<CcsWitness>> = (0..STEPS)
        .map(|i| vec![mul_gate_witness(gp.m, (i % 2) as u64)])
        .collect();
    let mul_rep = run_pipeline(&gp, &pp, &s_mul, &mul_steps).expect("multiplication pipeline");
    report_line("x·y=z gate ×4", &mul_rep);

    // 3. Ordinary computations compile to the SuperNeo CCS shape (recursive-verifier
    //    builder). Witnesses here are arbitrary field elements, so these are shown as
    //    satisfied CCS instances; folding them needs the augmented-witness decomposition.
    println!("\n── Computations as CCS instances (built via the circuit gadgets) ──────────\n");
    circuit_line("counter step (out=in+1)", &counter_circuit(41));
    circuit_line("fibonacci (a,b)→(b,a+b)", &fibonacci_circuit(8, 13));
    circuit_line("multiplier (x·y=z)", &multiplier_circuit(6, 7));

    println!("\n── Summary ────────────────────────────────────────────────────────────────\n");
    println!(
        "  Folded {STEPS} steps into a {}-instance accumulator, then produced a proof of\n  \
         ~{} B verified in {:.2?}. Its size is *constant in the number of folds* (the IVC\n  \
         length) — that is the succinctness here; it is not a compression of the witness\n  \
         (at these PoC params the FRI query phase dominates, so the proof is larger than\n  \
         the witness). The proof is post-quantum (hash-based BaseFold/FRI over Goldilocks),\n  \
         and the witness was never sent — only its Ajtai commitment, the constant-term\n  \
         evaluations (Theorem 6), and the ℓ∞ norm bound were proved.\n",
        id_rep.k, id_rep.proof_bytes, id_rep.verify_time
    );
    assert!(
        id_rep.verified && mul_rep.verified,
        "demo pipelines must verify"
    );
}

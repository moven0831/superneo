# SuperNeo (Rust PoC)

A from-scratch Rust implementation of **SuperNeo**, a post-quantum folding scheme for CCS
(Nguyen & Setty, *"Neo and SuperNeo"*, [IACR ePrint 2026/242](https://eprint.iacr.org/2026/242)).

Folding lets you take many constraint-system instances and combine them into one, so a
verifier only ever checks a single accumulated instance instead of re-checking each step.
SuperNeo does this with lattice (Ajtai) commitments, so it stays secure against quantum
attackers, and it runs its main check as one sum-check over a small field (Goldilocks).

This repo is a **proof of concept**: everything (field, ring, commitment, the folding
reductions, IVC, and the final proof) is built from the paper, not wrapped around an
existing library. It is meant to be read, run, and built on — not deployed. See
[Status](#status) for what is and isn't done, and [Scope & limits](#scope--limits) for the
honest caveats.

## Quick start

```sh
git clone <this-repo> superneo && cd superneo
cargo test --workspace --release          # run the test suite (99 tests)
cargo run -p superneo-snark --example demo --release   # see the full pipeline run
```

The demo folds a 4-step run, compresses it, verifies it, and prints timings — a good first
thing to read.

> The `context/` directory references the [Nightstream](https://github.com/LFDT-Nightstream/Nightstream)
> implementation as a git submodule. It is **reference material only** — you do **not** need
> it to build or test. Fetch it with `git submodule update --init` if you want to compare
> against it.

## Use it

The simplest entry point is `run_pipeline`: give it parameters, a constraint system, and a
list of per-step witnesses; it folds them, compresses the result, and verifies.

```rust
use superneo_commit::PublicParams;
use superneo_field::fp::Fp;
use superneo_fold::{CcsStructure, CcsWitness, GlobalParams};
use superneo_snark::run_pipeline;

// Parameters: ring-vector length n_R = 1, Ajtai matrix rank κ = 18.
let gp = GlobalParams::goldilocks(1);
let pp = PublicParams::setup_seeded([0u8; 32], 18, 1);

// The constraint system to fold. `identity` accepts any witness; build your own for
// real constraints (see "Extend it" below).
let s = CcsStructure::identity(gp.m);

// One low-norm witness per step (folding requires ‖z‖∞ < b = 2).
let steps: Vec<Vec<CcsWitness>> = (0..4)
    .map(|_| vec![CcsWitness { z: vec![Fp::ZERO; gp.n_f] }])
    .collect();

let report = run_pipeline(&gp, &pp, &s, &steps).unwrap();
assert!(report.verified);
println!("proof ~{} bytes, verified in {:?}", report.proof_bytes, report.verify_time);
```

If you want the steps separately, the same flow is:

```rust
use superneo_ivc::{prove_ivc_with_witnesses, verify_ivc};
use superneo_snark::{compress, verify};

let (ivc, final_wit) = prove_ivc_with_witnesses(&gp, &pp, &s, &steps)?;
verify_ivc(&gp, &pp, &s, &ivc)?;                         // checks every fold + the digest
let proof = compress(&gp, &pp, &s, &ivc.final_acc, &final_wit)?;
verify(&gp, &pp, &s, &ivc.final_acc, &proof)?;           // one succinct check at the end
```

A single fold step (no IVC) is `superneo_fold::fold` / `verify_fold`.

## Crates

A bottom-up stack — each crate depends only on the ones above it.

| Crate | What it does |
|-------|--------------|
| `superneo-field`  | Goldilocks field `Fp` and its degree-2 extension `K`. |
| `superneo-ring`   | The cyclotomic ring `R_F = F[X]/(X⁵⁴+X²⁷+1)`, coefficient maps, and the "bar" transform that turns matrix–vector products into ring products. |
| `superneo-commit` | The Ajtai commitment `Commit(A, z) = A·z` and the parameter set. |
| `superneo-fold`   | The constraint relations, the sum-check, the three folding reductions, and `fold` / `verify_fold`. |
| `superneo-ivc`    | The IVC loop (`prove_ivc` / `verify_ivc`) and a circuit builder that expresses the fold verifier as constraints. |
| `superneo-snark`  | The final proof: a hash-based (post-quantum) polynomial commitment plus the reduction that compresses an accumulator. The `pipeline` module and demo live here. |

## Extend it

- **Add your own constraints.** A constraint system here is a `CcsStructure` with four
  matrices `(I, A, B, C)` and the fixed rule `f = X₂·X₃ − X₄` (i.e. R1CS: `Az ∘ Bz = Cz`).
  Build one directly (see `superneo_snark::pipeline::multiplication_structure` for the
  `x·y = z` gate), or use the constraint builder in `superneo_ivc::circuit`:
  ```rust
  use superneo_ivc::circuit::{ConstraintSystem, Lc};
  let mut cs = ConstraintSystem::new();
  let x = cs.alloc(/* value */);
  let y = cs.alloc(/* value */);
  let z = cs.mul(&Lc::from_var(x), &Lc::from_var(y));   // enforces x·y = z
  let (structure, witness) = cs.finalize();             // → a CcsStructure you can fold
  ```
  `counter_circuit` / `fibonacci_circuit` / `multiplier_circuit` in `pipeline.rs` are small
  worked examples.
- **Folding needs low-norm witnesses** (`‖z‖∞ < b = 2`). Arbitrary field values won't fold
  directly; that's the lattice structure. The demo circuits above are *satisfied* but not
  folded for this reason (see [Scope & limits](#scope--limits)).
- **Swap pieces.** The Fiat–Shamir transcript (`superneo_fold::Transcript`, Blake3) and the
  polynomial commitment (`superneo_snark::basefold`) are the natural seams to replace.

## Status

All eight build milestones are implemented and tested end-to-end:

- **field / ring** — Goldilocks arithmetic and the bar transform (paper Thms 5–6), checked
  against direct computation.
- **commitment** — the Ajtai commitment and its parameter guard.
- **folding** — the sum-check and the `Π_DEC ∘ Π_RLC ∘ Π_CCS` fold (Thm 3), with round-trip,
  two-step, R1CS valid/invalid, and tamper tests.
- **IVC** — `prove_ivc` / `verify_ivc` over a step sequence, bound by a hash chain.
- **compression** — a BaseFold/FRI polynomial commitment + a Spartan-style reduction that
  proves the accumulator's commitment opening, evaluation claims, and norm bound; with
  round-trip and tamper tests.
- **recursive verifier circuit** — the full in-circuit `verify_fold` chain over advice
  challenges (Π_CCS sum-check + `Q(r′)` reconstruction → Π_RLC ring linear combination →
  Π_DEC additivity), expressed as constraints with Karatsuba and ring-multiplication
  gadgets; each reduction checked satisfied iff its native verifier accepts.
- **integration** — `run_pipeline`, the demo, and Criterion benches.

99 tests pass; `cargo clippy --all-targets -- -D warnings` is clean.

## Scope & limits

The soundness-relevant binding gaps an earlier audit surfaced have been closed; a few things
remain deliberately scoped for a PoC:

- **What the final proof binds (closed).** Compression binds the commitment opening, the norm
  bound, and **all `d` coefficients** of each evaluation claim `y_{i,j} ∈ R_K` — the Theorem-6
  constant term *and* the higher coefficients — via the bar-lifted matrix and the rotation
  identity, so a higher-coefficient forgery no longer passes.
- **Fiat–Shamir binding (closed).** The transcript binds the commitment key `A` and the CCS
  structure `s` (a verifier-key digest at each transcript origin, plus the relation inside
  Π_CCS), and challenge squeezes are domain-separated by label.
- **Parameters & IVC digest (closed).** `GlobalParams::goldilocks` enforces the
  folding-admissibility guard, `PublicParams::from_params` derives `κ` from the validated
  parameter set (so security can't be silently downgraded), and the IVC digest binds the full
  CE state `(c, r, y)`.
- **FRI soundness.** The query count is derived from a stated target (`SECURITY_BITS = 100`)
  under the proximity-gap assumption; the provable unique-decoding bound would need ~3× more
  queries. The proof is still *constant in the number of folds* but larger than the witness at
  these parameters, and the code isn't independently audited.
- **Recursive verifier (in progress).** The in-circuit Π_CCS verifier is complete over advice
  challenges: the sum-check verifier **and** the final `Q(r′)` reconstruction (the combined
  oracle `eq·(F + γ^K·NC) + γ^{2K+k}·Eval`, tied to the sum-check's reduced claim) are
  synthesized as constraints, and a circuit built from a real proof is satisfied **iff** the
  native verifier accepts (tamper an evaluation or a round polynomial → both reject). The
  Π_RLC ring linear combination (`c_out = Σ ρᵢ·cᵢ`, `y_out = Σ ρᵢ·y_{i,j}`, on in-circuit
  `R_F`/`R_K` multiplication gadgets) and the Π_DEC additivity check (`c = Σ bⁱ⁻¹cᵢ`,
  `y_j = Σ bⁱ⁻¹y_{i,j}`) are in-circuit too — so the **whole `verify_fold` chain** is
  synthesized and checked the same way. Still open: the Fiat–Shamir hashing that produces
  the challenges is native (α, γ, r′, ρ are passed as advice), and closing the loop needs
  the witness-decomposition step — the remaining effort (an arithmetization-friendly
  transcript for in-circuit Fiat–Shamir → loop closure). The ring-combination cost (`D²`
  base mults per product) means the Π_RLC step runs at scaled-down `κ, K, n_R`.

## Build & test

```sh
cargo build  --workspace
cargo test   --workspace --release
cargo clippy --workspace --all-targets -- -D warnings
cargo bench  --workspace                              # ring_mul, ajtai_commit, fold_step, compress, ...
cargo run -p superneo-snark --example demo --release
```

## License

Apache-2.0.

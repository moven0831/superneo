# SuperNeo (Rust PoC)

A clean-room Rust proof-of-concept of **SuperNeo**, the post-quantum lattice folding
scheme for CCS from Wilson Nguyen & Srinath Setty,
*"Neo and SuperNeo: Post-quantum folding with pay-per-bit costs over small fields"*
([IACR ePrint 2026/242](https://eprint.iacr.org/2026/242)).

SuperNeo is the lattice analog of HyperNova: the prover commits to a CCS witness with an
Ajtai/Module-SIS commitment that is *pay-per-bit*, runs a **single sum-check over a small
field extension** (Goldilocks), and folds many instance–witness pairs into one. It is the
first folding scheme to satisfy all six desiderata simultaneously — post-quantum,
pay-per-bit, field-native arithmetic, general (non-SIMD) constraints, small-field, and low
recursion overhead.

This is, to our knowledge, the first public implementation. It is **clean-room**: every
layer (field, ring, commitment, embedding, reductions) is implemented from the paper. The
only external crates are generic primitives (`rand`, `rayon`, `thiserror`, `serde`,
`blake3` for Fiat–Shamir). The parameterization is Goldilocks `q = 2⁶⁴−2³²+1` with the
cyclotomic ring `R_F = F[X]/(X⁵⁴+X²⁷+1)` (Appendix B.2).

## Workspace

A strict bottom-up dependency stack:

| Crate | Responsibility |
|-------|----------------|
| `superneo-field`  | Goldilocks `Fp` (Solinas reduction) and the degree-2 extension `K = F_q[u]/(u²−7)` |
| `superneo-ring`   | `R_F = F[X]/Φ₈₁` and `R_K`, coefficient maps, rotation/S-action, the **bar / inner-product transform** (Thm 5), the balanced ℓ∞ norm, and base-`b` decomposition |
| `superneo-commit` | Ajtai commitment `Commit(A,z)=A·z`, the Goldilocks parameter set with its `(K+k)·T·(b−1) < B` guard, and the strong sampling set `C` |
| `superneo-fold`   | Blake3 Fiat–Shamir transcript, CCS/CE relations, sum-check over `K`, the reductions `Π_CCS`/`Π_RLC`/`Π_DEC`, and the `fold` composition |
| `superneo-ivc`    | native IVC loop + the recursive verifier circuit (sum-check verifier as CCS gadgets) |
| `superneo-snark`  | final compression: BaseFold/FRI multilinear PCS + a Spartan-style reduction of the CE accumulator |

## Status

**Implemented and tested (the core folding scheme):**

- **Field** — Goldilocks arithmetic with Solinas reduction (cross-checked against `x mod q`
  over 200k random inputs), Fermat, and `Ext2` (7 is a verified non-residue; Frobenius
  `z^q = conj(z)`).
- **Ring + bar transform** — trinomial reduction (`X⁸¹ = 1`, `Φ` annihilates), the rotation
  identity `rot(a)·cf(b) = cf(ab)`, and the SuperNeo inner-product transform with
  `det(G) = ±1`, **Theorem 5** (`ct(bar(a)·b) = ⟨a,b⟩`, all 54² basis pairs + 20k randoms),
  and **Theorem 6** (`Mz = ct(M̄z)`). The `G⁻¹·G = I` post-condition is asserted always-on.
- **Commitment** — additive / ring-scalar / field-scalar homomorphism, the Appendix-D.8
  guard (`16200 < 16384`), and seeded determinism.
- **Folding** — sum-check (completeness + soundness), `Π_CCS` (the combined
  `Q = eq(X,α)(F + γ^K·NC) + γ^{2K+k}·Eval` sum-check), `Π_RLC`, `Π_DEC`, and
  `fold`/`verify_fold` = `Π_DEC ∘ Π_RLC ∘ Π_CCS` (**Theorem 3**). Tests cover a `K=2`,
  `k=14` round-trip, two-step composition, an `x·y=z` R1CS (valid folds; invalid is
  rejected), and tamper rejection.
- **Native IVC** — `prove_ivc`/`verify_ivc` chain the fold over a step sequence into a
  running `k`-instance accumulator, bound by a Construction-2 IO digest
  `digest_i = H(digest_{i-1}, i, acc_i)`. Tests cover a 4-step run and rejection of a
  tampered digest, a tampered fold proof, and a dropped step.
- **Final compression** — a post-quantum **BaseFold/FRI multilinear PCS** over Goldilocks
  (foldable Reed–Solomon code with `7^{(q−1)≫k}` roots of unity, radix-2 NTT encoding,
  Blake3 Merkle layers, and an interleaved sum-check ⊗ codeword-fold with a query phase;
  the load-bearing `encode(fold_coeffs)=fold_codeword(encode)` commutation is asserted),
  plus a **Spartan-style reduction** that batches the CE relation — the Ajtai opening
  `c = A·ẑ` (via the rotation identity), the **Theorem 6** evaluation claims
  `ct(y_j)=~(M_j z)(r)`, and the ℓ∞ norm bound — into two sum-checks over a single
  committed witness, opened by the PCS. `compress`/`verify` run over an IVC accumulator;
  tests cover the BaseFold round-trip (+5 tamper vectors), the accumulator compression
  round-trip (+tamper/high-norm rejection), and the full `IVC → compress → verify` path.

- **Recursive verifier circuit** — the folding verifier's dominant component expressed as
  CCS constraints. An R1CS builder (`⟨a,z⟩·⟨b,z⟩=⟨c,z⟩`, one-wire `z[0]=1`) finalizes into
  the same `t=4` CCS shape `fold` consumes; extension-field wires multiply via Karatsuba
  (3 base muls); and the **sum-check verifier** is synthesized in-circuit (per round:
  `g(0)+g(1)=claim` and `claim←g(r)` via in-circuit Lagrange). An honest degree-4, 6-round
  transcript synthesizes to **374 constraints** and is satisfied; it finalizes into a
  satisfying CCS instance (loop closure at the relation level); a tampered round polynomial
  or final claim is rejected.

79 tests pass; the workspace is `clippy -D warnings` clean.

**Scope notes (PoC).** The compression layer enforces the Ajtai opening, the (constant-term)
Theorem-6 evaluations, and the norm bound; the higher `R_K` coefficients of the evaluation
claims are the documented residual (same linear-claim shape, with the bar-lifted matrix).
The PCS uses illustrative parameters (rate `1/4`, 32 queries) and is not security-audited.
The recursive circuit synthesizes the sum-check verifier with Fiat–Shamir challenges as
advice (Blake3-in-circuit is out of PoC range); the Π_CCS final-`Q` reconstruction, the
Π_RLC/Π_DEC gadgets, and the augmented-witness decomposition needed to re-fold the
synthesized instance through the norm-bounded scheme are documented residuals.

**Planned (remaining for full end-to-end, see the implementation plan):**

- **integration** — an end-to-end demo binary and Criterion benches (M8), plus the
  recursive-circuit residuals above.

## Build & test

```sh
cargo build --workspace
cargo test  --workspace --release
cargo clippy --workspace --all-targets -- -D warnings
```

## License

Apache-2.0.

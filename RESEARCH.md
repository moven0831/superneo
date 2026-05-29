# SuperNeo → Rust PoC: Research Dossier

**Target paper:** Wilson Nguyen & Srinath Setty (Microsoft Research),
*"Neo and SuperNeo: Post-quantum folding with pay-per-bit costs over small fields"*,
IACR ePrint **2026/242**. https://eprint.iacr.org/2026/242

**What it is:** A *lattice-based folding scheme* for CCS — the post-quantum analog of
HyperNova. The prover commits to a CCS witness with an Ajtai (Module-SIS) commitment that
is pay-per-bit, runs a *single* sum-check over a small field extension (Goldilocks-class),
and folds. SuperNeo extends Neo to drop the SIMD/data-parallel requirement and is the first
folding scheme to satisfy all six desiderata: post-quantum (PQ), pay-per-bit (PB),
field-native arithmetic (FN), general constraints (CS), small-field (SF), low recursion (LR).

There is **no public reference implementation of Neo or SuperNeo yet** — this PoC would be
the first. The nearest existing Rust code is LatticeFold (Nethermind), Sonobe (folding
architecture), and lattirust (lattice ring arithmetic).

---

## 0. Primary sources (read first)

| Ref | Paper | ePrint | Why |
|----|-------|--------|-----|
| [2026/242] | Neo **and SuperNeo** (target) | https://eprint.iacr.org/2026/242 | The construction to implement |
| [71] | **Neo** (predecessor; SuperNeo subsumes it) | https://eprint.iacr.org/2025/294 | Original embedding + sum-check/eval-homomorphism; SuperNeo builds directly on it |
| [42] | **Cyclo** — Lipmaa et al., reinterpretation of Neo + first benchmark | ePrint 2026 (Cyclo) | "Ultimately equivalent to Neo"; an alternative presentation + initial impl notes |
| [56]/[59] | **HyperNova** — Kothapalli & Setty | CRYPTO'24 / 2023/573 | The group-based blueprint SuperNeo mirrors; defines the folding shape |

---

## 1. Core primitives (must implement)

### Commitment: Ajtai / Module-SIS
- [2] Ajtai, *Generating hard instances of lattice problems*, STOC 1996 — the commitment.
- [76] Peikert & Rosen, *Efficient collision-resistant hashing from worst-case assumptions
  on cyclic lattices*, TCC 2006.
- [64] Lyubashevsky & Micciancio, *Generalized compact knapsacks are collision resistant*,
  ICALP 2006 (Ring-SIS).
- [61] Langlois & Stehlé, *Worst-case to average-case reductions for module lattices*, 2015
  (Module-SIS hardness basis).
- [4] Albrecht, Player, Scott, *On the concrete hardness of LWE* — the **lattice estimator**
  the paper uses for parameter selection (Appendix B/D.8).

### Sum-check + CCS
- [63] Lund, Fortnow, Karloff, Nisan, *Algebraic methods for interactive proof systems*,
  FOCS 1990 — the sum-check protocol.
- [82] Thaler, *The sum-check protocol* (lecture notes) — practical reference.
- [80] Setty, Thaler, Wahby, *Customizable Constraint Systems (CCS)*, 2023/552 — the relation
  (generalizes R1CS / Plonkish / AIR).
- [77] Schwartz, *Fast probabilistic algorithms for verification of polynomial identities*,
  1980 (Schwartz–Zippel, soundness).
- [79] Setty, *Spartan*, CRYPTO 2020 / 2019/550 — sum-check-based zkSNARK; used to **compress**
  the final SuperNeo accumulator with a FRI-based PCS.

### Lattice algebra / ring switching (the "Neo embedding" machinery)
- [19] Bootle, Lyubashevsky, Seiler, *Algebraic techniques for short(er) exact lattice ZK*,
  CRYPTO 2019 — reduces norm constraints to ring constraints; the Neo embedding extends this.
- [36] Esgin, Nguyen, Seiler, *Practical exact proofs from lattices*, ASIACRYPT 2020 —
  **inner-product trick for cyclotomic rings** (used in the SuperNeo embedding).
- [65] Lyubashevsky, Nguyen, Plançon, *Lattice-based ZK: shorter, simpler, more general*,
  CRYPTO 2022 — inner-product / norm-via-ring-product technique.
- [66] Lyubashevsky & Seiler, *Short, invertible elements in partially splitting cyclotomic
  rings*, EUROCRYPT 2018 — **challenge sets / invertibility** (parameter selection).
- [3] Albrecht & Lai, *Subtractive sets over cyclotomic rings*, CRYPTO 2021 — low-norm
  challenge sets (strong sampling set 𝒞).
- [28] Chen et al., *Verifiable computing for approximate computation*, 2019/762 — low-norm
  ring challenge elements.
- [12] Beullens & Seiler, *LaBRADOR: Compact proofs for R1CS from Module-SIS*, CRYPTO 2023 —
  recurring lattice proof techniques; NTT-minimization.

### NTT (needed for ring arithmetic even though sum-check avoids it)
- [62] Longa & Naehrig, *Speeding up the NTT for ideal lattices*, CANS 2016.
- [78] Seiler, *Faster AVX2-optimized NTT multiplication for ring-LWE*, 2018.
- [5] Alkim, Ducas, Pöppelmann, Schwabe, *NewHope* — NTT embedding reference.

### Fields
- [75] Polygon Zero, **Plonky2** — Goldilocks field `q = 2^64 − 2^32 + 1`, high 2-adicity,
  FRI. The paper's recommended SF field; also "Almost Goldilocks" `q = (2^64−2^32+1)−32`.
- [44] Grassi et al., *Poseidon* — recursion-friendly hash (paper aims to *avoid* it, but
  relevant for hash-based comparison / Fiat–Shamir).

---

## 2. Folding / IVC / PCD framework (compiler layer)

- [58] Nova — Kothapalli, Setty, Tzialla, CRYPTO 2022 (recursion from folding).
- [54] SuperNova — Kothapalli & Setty (non-uniform IVC).
- [57] NeutronNova — Kothapalli & Setty (folding from zero-check).
- [55] CycleFold — Kothapalli & Setty, 2023/1192 (cycle-of-curves trick; cross-field reuse).
- [85] MicroNova — Zhao, Setty et al., IEEE S&P 2025 — efficient on-chain verification;
  HyperKZG-style accumulator compression (SuperNeo mirrors this with Spartan+FRI).
- [20] Protostar — Bünz & Chen, 2023/620 (accumulation for special-sound protocols).
- [35] Protogalaxy — Eagen & Gabizon (folding many instances).
- [83] Valiant, *Incrementally verifiable computation (IVC)*, TCC 2008.
- [13] Bitansky, Canetti, Chiesa, Tromer, *Recursive composition / PCD*, STOC 2013.
- [86] Zhou, Zhang, Dong, *Proof-carrying data from multi-folding schemes*, 2023/1282 —
  the PCD compiler SuperNeo plugs into.
- [52]/[53] Kothapalli & Parno, *Algebraic reductions of knowledge*, CRYPTO 2023 — the
  abstraction SuperNeo's **interactive reductions** framework generalizes.
- [51] Kothapalli, *A Theory of Composition for Proofs of Knowledge*, PhD thesis 2024.

---

## 3. Direct competitors / technique siblings (lattice folding & SNARKs)

- [14]/[15] **LatticeFold** — Boneh & Chen, 2024/257 — first lattice folding (NTT embedding;
  the thing SuperNeo improves on). **Has a Rust impl** (see §5).
- [16]/[17] **LatticeFold+** — Boneh & Chen, 2025/247 — faster/shorter; re-interprets Neo's
  technique as "tensor-of-rings". **Has a Rust impl**.
- [37]/[38] **Lova** — Fenzi, Knabenhans, Nguyen, Pham, 2024/1964 — folding from *unstructured*
  SIS (safest assumption, slower; subset-sum only).
- [60] **SALSAA** — Kuriyama, Lai, Osadnik, Tucci, 2025/2124 — sum-check-aided lattice
  arguments (uses vSIS assumption).
- [26] **Symphony** — Chen, 2025/1905 — lattice SNARK from high-arity folding (tensor-of-rings).
- [49] **RoK, paper, SISsors** + [50] **RoK and Roll** — Klooß, Lai, Nguyen, Osadnik —
  lattice succinct-argument toolkit / random projection.
- [69] **Hachi** — Nguyen, O'Rourk, Zhang, 2026/156 — lattice **multilinear PCS** over
  extension fields (natural PCS for SuperNeo accumulator compression).
- [70] **Greyhound** — Nguyen & Seiler, CRYPTO 2024 — fast lattice polynomial commitments.
- [32] Cini, Lai, Malavolta, *Lattice succinct arguments from vSIS* + [46] Jyrkinen & Lai,
  *Vanishing SIS revisited*, PKC 2025 — the vSIS assumption.
- [7]/[8]/[9] Attema et al. — compressed Σ-protocol theory for lattices / adaptive special
  soundness / product proofs (special-soundness machinery for the security proofs).
- [33] Cini, Malavolta, Nguyen, Wee, *Polynomial commitments from lattices*, CRYPTO 2024.

---

## 4. Motivating application (post-quantum signature aggregation)

- [34] Corater & Setty, *Post-Quantum Signature Aggregation: A Folding Approach*,
  ethresear.ch 2025 — the Ethereum use case. https://ethresear.ch/t/23639
- [18] Boneh, Lynn, Shacham, *Short signatures from the Weil pairing* (BLS — what's being
  replaced).
- [45] Hülsing et al., **XMSS**, RFC 8391 — the hash-based PQ signature to aggregate.
- [1] Aardal et al., *Aggregating Falcon signatures with LaBRADOR*, CRYPTO 2024.

---

## 5. Existing Rust implementations to reuse / study

| Repo | What you get | Reuse for |
|------|--------------|-----------|
| **lattirust/lattirust** (`lattirust_arithmetic`) https://github.com/lattirust/lattirust | "arkworks for lattices": cyclotomic rings ℤq[X]/(Xⁿ+1), NTT, norms, challenge sets, linear algebra, arkworks-compatible | **Ring + commitment base layer** (closest foundation) |
| **NethermindEth/latticefold** https://github.com/NethermindEth/latticefold | LatticeFold + LatticeFold+ in Rust on arkworks; Ajtai commitments, NTT-packed rings, range-proof folding, rayon | Reference for lattice folding mechanics, decomposition, norm checks |
| **privacy-ethereum/sonobe** https://github.com/privacy-ethereum/sonobe ( docs: https://sonobe.pse.dev ) | Modular folding library: Nova, HyperNova, CycleFold; IVC loop, decider, Solidity/EVM verifier; arkworks `FCircuit` trait | **Folding/IVC architecture template**; how to structure the API + decider |
| **microsoft/Spartan** & **Spartan2** https://github.com/microsoft/Spartan2 | Setty's sum-check zkSNARK, PCS-generic, R1CS/CCS, Rust | Sum-check + CCS + final-proof compression backend |
| **microsoft/Nova** https://github.com/Microsoft/Nova | Canonical folding impl (Setty) | Folding loop reference / idioms |
| **ajtai** crate https://docs.rs/ajtai | NTT cyclic/negacyclic convolution targeting Ajtai commitments | Fast commitment arithmetic |
| **arkworks-rs** https://github.com/arkworks-rs | `ark-ff`, `ark-poly`, `ark-std`, `ark-crypto-primitives` (fields, multilinear poly, sumcheck, FS transcript) | Field/poly/transcript primitives |
| **plonky2** (Polygon Zero) | Goldilocks field + FRI in Rust | Field impl + FRI-based final compression |
| **lurk-lab/awesome-folding** https://github.com/lurk-lab/awesome-folding | Curated index of folding schemes + code | Discovery / staying current |

**Lattice hardness tooling (not Rust):** the official **lattice-estimator** (Albrecht et al.,
Sage/Python) for choosing Module-SIS parameters — the paper's Appendix B/D.8 ships Sage
scripts for the inversion bounds and challenge-set sizing.

---

## 6. Suggested build order (PoC)

1. **Field layer** — Goldilocks (+ degree-2 extension 𝕂 for sum-check soundness). Use
   plonky2/arkworks field, or port. (Paper: Appendix B params.)
2. **Ring layer** — `R𝔽 = 𝔽[X]/(Φ(X))` for a chosen cyclotomic Φ (e.g. `X^54+X^27+1`, d=54
   for Goldilocks, or `X^64+1` for Almost-Goldilocks). Base on **lattirust**.
3. **Ajtai/Module-SIS commitment** — random `A ∈ R𝔽^{κ×n}`, `Commit(z)=A·z`, ℓ∞-norm binding.
   Params via lattice-estimator. (Refs [2],[76],[4].)
4. **Neo / SuperNeo embedding** — `Coeff(z)` placing field vectors in coefficient slots;
   prove the **evaluation homomorphism** (Thm 1 / Thm 7). This is the novel core.
5. **Sum-check over the field** — combined CCS-check + norm-check + eval-check polynomial
   `Q(X)` (paper §7.3). Use arkworks/Spartan sum-check.
6. **Interactive reductions** — Π_CCS (sum-check), Π_RLC (random linear combination),
   Π_DEC (decomposition); compose via strong⊕weak composition theorem (paper §7, Thm 8).
7. **IVC/PCD wrapper** — mirror Sonobe's loop; recursive verifier circuit in CCS.
8. **Final compression** — Spartan + FRI-based PCS (or Hachi lattice PCS) over the accumulator.

---

*Generated 2026-05-29. Paper PDF cached during research; all ePrint refs verifiable at
https://eprint.iacr.org/<id>.*

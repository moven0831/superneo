# `context/` — SuperNeo reference corpus

Local copy of the papers and codebases needed to build SuperNeo in Rust.
See `../RESEARCH.md` for the architecture, citation→component mapping, and build order.
Filenames use `YYYY-NNN_short-name` = IACR ePrint `YYYY/NNN` (verify at `https://eprint.iacr.org/YYYY/NNN`).

## `papers/` (46 PDFs + XMSS RFC)

### Primary
- [`2026-242_neo-superneo`](https://eprint.iacr.org/2026/242) — **the target paper** (Nguyen & Setty).
- [`2025-294_neo`](https://eprint.iacr.org/2025/294) — predecessor; SuperNeo subsumes it.
- [`2026-359_cyclo`](https://eprint.iacr.org/2026/359) — Eurocrypt'26 reinterpretation of Neo + LatticeFold+; "equivalent to Neo".
- [`2023-573_hypernova`](https://eprint.iacr.org/2023/573) — group-based blueprint SuperNeo mirrors.

### Core primitives
- [`2023-552_ccs`](https://eprint.iacr.org/2023/552) — Customizable Constraint Systems (the relation).
- [`2019-550_spartan`](https://eprint.iacr.org/2019/550) — sum-check zkSNARK; final-proof compression backend.
- [`2017-523_lyubashevsky-seiler-invertible-elements`](https://eprint.iacr.org/2017/523) — **LS18**, challenge sets/invertibility (paper's Sage script cites this).
- [`2021-202_albrecht-lai-subtractive-sets`](https://eprint.iacr.org/2021/202) — **AL21**, low-norm challenge sets (cited in Sage script).
- [`2015-046_albrecht-concrete-hardness-lwe`](https://eprint.iacr.org/2015/046) — lattice-estimator hardness basis (param selection).
- [`2012-090_langlois-stehle-module-lattices`](https://eprint.iacr.org/2012/090) — Module-SIS worst→average-case hardness.
- `thaler_sumcheck_notes.pdf` — sum-check protocol reference notes.

### The embedding / ring-switching machinery (novel core)
- [`2019-532_bootle-algebraic-techniques`](https://eprint.iacr.org/2019/532) — norm→ring constraint reduction the Neo embedding extends.
- [`2020-518_esgin-practical-exact-proofs`](https://eprint.iacr.org/2020/518) — inner-product trick for cyclotomic rings.
- [`2022-284_lyubashevsky-zk-shorter-simpler`](https://eprint.iacr.org/2022/284) — inner-product / norm-via-ring-product.
- [`2019-762_chen-verifiable-approx-computation`](https://eprint.iacr.org/2019/762) — low-norm ring challenges.
- [`2022-1341_labrador`](https://eprint.iacr.org/2022/1341) — Module-SIS R1CS proofs; NTT-minimization techniques.

### NTT / fields
- [`2016-504_longa-naehrig-ntt`](https://eprint.iacr.org/2016/504), [`2018-039_seiler-avx2-ntt`](https://eprint.iacr.org/2018/039), [`2015-1092_newhope`](https://eprint.iacr.org/2015/1092) — NTT.
- [`2019-458_poseidon`](https://eprint.iacr.org/2019/458) — recursion-friendly hash (comparison / Fiat–Shamir).
- `plonky2_paper.pdf` — Goldilocks field + FRI (the SF field SuperNeo recommends).

### Folding / IVC / PCD compiler
- [`2021-370_nova`](https://eprint.iacr.org/2021/370), [`2022-1758_supernova`](https://eprint.iacr.org/2022/1758), [`2024-1606_neutronnova`](https://eprint.iacr.org/2024/1606), [`2023-1192_cyclefold`](https://eprint.iacr.org/2023/1192),
  [`2024-2099_micronova`](https://eprint.iacr.org/2024/2099) — the Nova family (MicroNova = on-chain verification SuperNeo mirrors).
- [`2023-620_protostar`](https://eprint.iacr.org/2023/620), [`2023-1106_protogalaxy`](https://eprint.iacr.org/2023/1106) — accumulation / multi-instance folding.
- [`2023-1282_pcd-multifolding`](https://eprint.iacr.org/2023/1282) — PCD compiler SuperNeo plugs into.
- [`2023-366_reductions-of-knowledge`](https://eprint.iacr.org/2023/366) — abstraction SuperNeo's "interactive reductions" generalizes.
- [`2024-1043_kothapalli-composition-thesis`](https://eprint.iacr.org/2024/1043) — Theory of Composition for Proofs of Knowledge.

### Lattice siblings / competitors / technique sources
- [`2024-257_latticefold`](https://eprint.iacr.org/2024/257), [`2025-247_latticefold-plus`](https://eprint.iacr.org/2025/247), [`2026-721_latticefold-plus-l2-norm`](https://eprint.iacr.org/2026/721) — LatticeFold line.
- [`2024-1964_lova`](https://eprint.iacr.org/2024/1964) — folding from unstructured SIS.
- [`2025-2124_salsaa`](https://eprint.iacr.org/2025/2124) — sum-check-aided lattice arguments (vSIS).
- [`2025-1905_symphony`](https://eprint.iacr.org/2025/1905) — lattice SNARK from high-arity folding.
- [`2024-1972_rok-paper-sissors`](https://eprint.iacr.org/2024/1972), [`2025-1086_rok-and-roll`](https://eprint.iacr.org/2025/1086) — lattice succinct-argument toolkit / random projection.
- [`2026-156_hachi`](https://eprint.iacr.org/2026/156) — lattice multilinear PCS (natural accumulator-compression PCS).
- [`2024-1293_greyhound`](https://eprint.iacr.org/2024/1293) — fast lattice polynomial commitments.
- [`2023-941_cini-lattice-vsis`](https://eprint.iacr.org/2023/941), [`2024-281_cini-malavolta-pq-pcs`](https://eprint.iacr.org/2024/281) — vSIS / lattice PCS.
- [`2020-517_attema-product-proofs`](https://eprint.iacr.org/2020/517), [`2024-2038_attema-adaptive-soundness`](https://eprint.iacr.org/2024/2038) — special-soundness machinery.
- [`2024-1586_whir`](https://eprint.iacr.org/2024/1586), [`2023-1705_basefold`](https://eprint.iacr.org/2023/1705) — field-agnostic PCS backends.

### Application
- [`2024-311_falcon-labrador`](https://eprint.iacr.org/2024/311) — aggregating Falcon sigs with LaBRADOR.
- `rfc8391_xmss.txt` — XMSS (the PQ signature to aggregate).

## `implementations/` (8 repos, shallow clones)
- `lattirust/` — "arkworks for lattices": cyclotomic rings, NTT, norms, challenge sets. **Ring/commitment base.**
- `latticefold/` — LatticeFold(+) in Rust on arkworks. Lattice folding mechanics reference.
- `sonobe/` — Nova/HyperNova/CycleFold + IVC loop + decider + EVM verifier. **Folding architecture template.**
- `Spartan/`, `Spartan2/` — Setty's sum-check/CCS engine + proof compression.
- `Nova/` — canonical folding impl idioms.
- `plonky2/` — Goldilocks field + FRI.
- `awesome-folding/` — curated index of folding schemes + code.

## Not included (no open PDF — cite by DOI)
Ajtai STOC'96 (commitment, doi:10.1145/237814.237838); Lund–Fortnow–Karloff–Nisan FOCS'90
(sum-check, doi:10.1109/FSCS.1990.89518); Boneh–Lynn–Shacham ASIACRYPT'01 (BLS); Valiant
TCC'08 (IVC, doi:10.1007/978-3-540-78524-8_1); Peikert–Rosen TCC'06 (doi:10.1007/11681878_8);
Lyubashevsky–Micciancio ICALP'06 (Ring-SIS, doi:10.1007/11787006_13); Bitansky et al. STOC'13
(PCD); Schwartz J.ACM'80 (Schwartz–Zippel). Also: `arkworks-rs` crates (`ark-ff`/`ark-poly`/
`ark-crypto-primitives`) and the `ajtai` crate — pull via `cargo`, not vendored here.

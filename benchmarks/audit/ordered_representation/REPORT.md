# WP-C Ordered Angular Representation Report (SYNTHETIC)

Evaluation of a candidate representation (residue-ordered (cos τ, sin τ) features with same-sequence distances) on SYNTHETIC fixtures. It covers representation fidelity and fixed-fixture retrieval only. **No mode or state detection is performed or claimed, no production behaviour is changed, and no real-MD data are used.** Protocol: `docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md` (amendments §11, recorded before execution). Every number below is read from `results/SYNTHETIC_ordered_representation_*.json`.

## Provenance

- Git HEAD `295d739a783cc1d178131c4c61293042a0c56f7d` (working tree with 17 uncommitted entries at run time). Python 3.14.4, NumPy 2.5.3, rustc 1.98.1 (48a229cea 2026-09-01).
- Production binding (T1): `/home/artem/TopoFold/benchmarks/audit/angular_characterization/build/venv/lib/python3.14/site-packages/topofold/topofold.cpython-314-x86_64-linux-gnu.so`, SHA-256 `732909cc99ec88835f0b8a1f60152b7f134ae4a353d152be7bd9f048e232ba01`. Identical to the WP-A freshly built extension: **True**. Production sources unchanged since that build: **True** (compared: `Cargo.lock`, `crates/topofold-core/src/bimodality.rs`, `crates/topofold-core/src/circular.rs`, `crates/topofold-core/src/geometry.rs`, `crates/topofold-core/src/lib.rs`, `crates/topofold-core/src/types.rs`, `crates/topofold-python/src/lib.rs`).
- WP-A inputs read from: ['wpa_build_file']. Each was verified against the frozen WP-A manifest (SHA-256 `17d5bd2827688e917eff84f7fed9404929214d6549a7fadec03ad20fab08c928`); nothing was written into the WP-A directory, and the WP-A pipeline was not rerun.
- Historical scalar comparator: the unchanged WP-A Rust example (production functions), writing into the WP-C `build/`. Outputs: JSON `fa65ed5d659b21273bcfbea57b55be07ab637d067ccaff187483bf04af8369e5`, arrays `e446b76658b4f4ff5c200351a6c13240b36498674f416c1bfaf773053cfeeba9`.
- Commands and exit codes:

| Exit | Command |
|---:|---|
| 0 | `git rev-parse HEAD` |
| 0 | `git status --porcelain` |
| 0 | `rustc --version` |
| 0 | `cargo run --release -q -p topofold-core --example angular_characterization -- benchmarks/audit/ordered_representation/build/comparator/compa` |

Paths per case:
- Every case: T0 (prescribed), T2 (NumPy recovery from f64 coordinates) and T1 (production binding, float32 input).
- Retrieval fixtures: additionally the scalar comparator and the Cα-distance baseline.
- Case lists and families are in `results` → `inputs.case_families`.

## 1. Harness integrity (V1)

| Component | Result |
|---|---|
| coverage | PASS |
| input_hashes | PASS |
| binding_identity | PASS |
| recovery | PASS |
| comparator_integrity | PASS |
| tie_self_check | PASS |

Recovery: max T2 error 4.44e-16 rad; max T1 (float32) error 3.77e-07 rad.

## 2. Fidelity and continuity (V2)

**Feature-contract checks** (§5a; distances alone cannot detect these conventions):

| Check | Result |
|---|---|
| FC1 values (asymmetric known torsions, order, cos/sin entries) | PASS |
| FC2 window offset (i=3, W=6 uses τ_3, τ_4, τ_5) | PASS |
| FC3 residue mapping (bonds (4,5), (5,6), (6,7)) | PASS |
| FC4 chirality (sine entry of τ = -0.7 is -sin 0.7) | PASS |

**Distance-geometry checks** (§5b):

| Check | Measured | Result |
|---|---|---|
| AD1 chord identity | max error 3.3e-16 | PASS |
| AD2 (2/π)·D_circ ≤ D_ch ≤ D_circ on every retrieval matrix | failures 0 | PASS |
| AD3 F1a_s0 A vs B | D_ch 1.000000000000000, D_circ 1.047197551196598, scalar 0.000e+00 | PASS |
| AD3 F1b_s0 A vs B | D_ch 0.632455532033676, D_circ 0.662305884386407, scalar 0.000e+00 | PASS |
| AD3 F1c_s0 A vs B | D_ch 1.000000000000000, D_circ 1.047197551196598, scalar 1.047e+00 | PASS |
| AD4 paired-shift invariance (F2) | max diff 1.8e-15 | PASS |
| AD5 Lipschitz continuity (F3) | max ‖Δx‖/δ 1.000000129 | PASS |
| AD6 mirror (N2) | D_ch error 0.0e+00; Cα-distance difference 0.0e+00 Å | PASS |
| AD7 branch | D_ch 2.00e-12 | PASS |
| AD8 F6b equal torsions | max error 2.2e-16 | PASS |
| Recovery-aware distances (T2, T1) | failures 0 | PASS |

**Continuity across ±ε** (F3 noiseless, frame 0). The candidate varies continuously where the historical scalar summary jumps:

| ε | D_ch(+ε, −ε) | expected √2·sin ε | scalar comparator (+ε vs −ε) |
|---|---:|---:|---:|
| 1e-01 | 1.411858e-01 | 1.411858e-01 | 3.041593 |
| 1e-03 | 1.414213e-03 | 1.414213e-03 | 3.140593 |
| 1e-06 | 1.414214e-06 | 1.414214e-06 | 3.141592 |
| 1e-09 | 1.414214e-09 | 1.414214e-09 | 3.141593 |

## 3. Fixed-fixture retrieval (V3, descriptive)

Expected neighbour purity at k = 10 on T0, with label-independent tie handling. "ties" is the fraction of queries whose k-th boundary was tied. These are scores on fixed fixtures, with no confidence intervals and no population inference.

| Fixture | chance | D_ch | D_circ | scalar | Cα-dist | D_κ | ties D_ch / scalar / Cα / κ | D_ch (T1) |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| F1a_s0 | 0.499 | 1.000 | 1.000 | 0.499 | 1.000 | 0.499 | 1.00 / 1.00 / 1.00 / 1.00 | 1.000 |
| F1a_s2deg | 0.499 | 1.000 | 1.000 | 0.511 | 1.000 | 0.499 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| F1b_s0 | 0.499 | 1.000 | 1.000 | 0.499 | 1.000 | 0.499 | 1.00 / 1.00 / 1.00 / 1.00 | 1.000 |
| F1b_s2deg | 0.499 | 1.000 | 1.000 | 0.497 | 1.000 | 0.499 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| F1c_s0 | 0.499 | 1.000 | 1.000 | 1.000 | 1.000 | 0.499 | 1.00 / 1.00 / 1.00 / 1.00 | 1.000 |
| F1c_s2deg | 0.499 | 1.000 | 1.000 | 1.000 | 1.000 | 0.499 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| N1_e+1e-06 | 0.500 | 1.000 | 1.000 | 1.000 | 1.000 | 0.500 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| N1_e-1e-06 | 0.500 | 1.000 | 1.000 | 1.000 | 1.000 | 0.500 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| N1_e+1e-03 | 0.500 | 1.000 | 1.000 | 1.000 | 1.000 | 0.500 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| N1_e-1e-03 | 0.500 | 1.000 | 1.000 | 1.000 | 1.000 | 0.500 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| N1_e+1e-01 | 0.500 | 1.000 | 1.000 | 1.000 | 1.000 | 0.500 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| N1_e-1e-01 | 0.500 | 1.000 | 1.000 | 1.000 | 1.000 | 0.500 | 0.00 / 0.00 / 0.00 / 1.00 | 1.000 |
| N2_pooled_F1a_s0_mirror | 0.250 | 1.000 | 1.000 | 0.499 | 0.499 | 0.250 | 1.00 / 1.00 / 1.00 / 1.00 | 1.000 |

F4 (α = 0; 44 fixtures), range across fixtures:

| Representation | min purity (k=10) | max purity (k=10) | min / max tie frequency |
|---|---:|---:|---|
| D_ch | 1.000 | 1.000 | 0.00 / 1.00 |
| D_circ | 1.000 | 1.000 | 0.00 / 1.00 |
| scalar_comparator | 1.000 | 1.000 | 0.00 / 1.00 |
| ca_distance | 1.000 | 1.000 | 0.00 / 1.00 |
| D_kappa | 0.499 | 0.505 | 1.00 / 1.00 |

k = 1 (D_ch, T0): min purity over all retrieval fixtures 1.000.

**Differences from comparators** (margin 0.05, pre-declared):
- Material deficits vs the scalar comparator: none.
- Material deficits vs the Cα-distance baseline: none.
- D_ch above the scalar comparator by more than the margin: F1a_s0 (+0.501), F1a_s2deg (+0.489), F1b_s0 (+0.501), F1b_s2deg (+0.503), N2_pooled_F1a_s0_mirror (+0.501). These are the designed coordinate-wise collisions (F1a/F1b) and the mirror set.
- D_ch above the Cα-distance baseline by more than the margin: N2_pooled_F1a_s0_mirror (+0.501). This is the pre-stated chirality difference: the baseline is reflection-invariant.
- N2 pooled set: 4 construction labels (state × chirality); overall chance purity 0.250. The Cα-distance and scalar scores of 0.499 and 0.499 are **not** chance-level overall. They are consistent with resolving one label factor and leaving pairs of classes unresolved: the Cα-distance baseline cannot separate mirror images, and the scalar summary cannot separate the two states (the per-frame mean is 60° in both, with sign reversed by the mirror).
- On all other fixtures, D_ch, D_circ and the Cα-distance baseline reach the same purity. **The fixtures saturate**: they confirm the designed differences but cannot rank the candidate against the Cα-distance baseline.

**Historical scan on the paired ±ε control N1** (production BC; shared assignment and noise, so only ε differs):

| ε | scan BC (R2) |
|---|---:|
| -1e-01 | 0.9905 |
| -1e-03 | 0.9906 |
| -1e-06 | 0.5967 |
| +1e-06 | 0.5967 |
| +1e-03 | 0.4863 |
| +1e-01 | 0.9905 |

With draws paired, the scan output is the same for ε = ±1e-6 and for ±1e-1, but differs between ε = +1e-3 and −1e-3. This one paired draw did **not reproduce** the ±1e-6 difference seen in WP-A. Because the WP-A configurations were unpaired, that earlier comparison cannot isolate a sign effect from draw variation, and a single paired draw does not establish that the earlier difference was caused entirely by draw variation. The ±1e-3 difference is likewise observed for this one paired draw only; its mechanism is not established. The ordered representation is not involved in these values.

## 4. Limitations

- **Constant curvature (κ = 1.2 rad).** D_κ carries no information here (purity = chance, ties everywhere). The κ block is untested.
- **Fixed bond lengths (3.8 Å)** and single-window traces. Chain breaks and inter-window context are not represented.
- **Constructed labels and simple geometries.** Two well-separated states, mostly W = 5 (F1b: W = 8), so retrieval saturates for every representation except where an information difference was designed in.
- Descriptive scores on fixed fixtures only. No generalization, no superiority claim, no statement about real proteins.

## 5. Decision record

- V1 harness integrity: **PASS**
- V2 fidelity: **PASS**
- V3 deficits vs scalar comparator: none; vs Cα-distance baseline: none
- **Decision (protocol §9 ordered rules): RETAIN: eligible for a real-MD pilot.**

Pilot eligibility means only the following:
- The candidate removed the demonstrated collisions and branch sensitivity.
- It lost nothing relative to the component it would replace, on these fixtures.

It is **not** evidence that the candidate outperforms local Cα-distance features. The fixtures could not discriminate the two outside the designed chirality case. Conditions for the real-MD pilot:
- the Cα-distance baseline and the historical scan are included as comparators;
- a task is used that can separate representations (e.g. retrieval against independent G1 states);
- WP-B gates G0–G1 pass independently of this result.

## 6. Test corrections during execution (recorded, contracts unchanged)

Two harness tests failed on first execution because of test-design errors, not harness defects. Both were corrected without weakening their contracts:
1. `test_distance_checks_cannot_see_consistent_permutation` demanded bitwise equality of distances computed in different summation orders. It now requires equality to 1e-15, which is still the invariance that motivates amendment A.
2. `test_fully_tied_set_gives_chance_not_index_order` used an example where index-based tie-breaking happened to land near chance. It now checks the contract directly: the candidate scores exactly chance under two frame orderings, while index-based tie-breaking changes with the ordering.
## 7. Post-run editorial changes

Wording corrections made after the run (no rerun, results unchanged) are logged in `results/SYNTHETIC_ordered_representation_editorial_log.json`.


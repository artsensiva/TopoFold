# Benchmarks: data provenance

**None of the scripts in this directory analyses a real molecular-dynamics trajectory.** Most ensembles are generated or interpolated by the script itself, with the conformational states and the labels built in. These are *controlled synthetic tests*: useful for checking that the code recovers a known injected signal, but they are not evidence about real proteins. One benchmark (XCL1) uses real experimental structures from the PDB compared statically. Earlier versions of these docstrings and of the README described several of them as "authentic MD"; that description was wrong (see [`docs/KNOWN_ISSUES.md`](../docs/KNOWN_ISSUES.md), C1).

| Script | What the ensemble actually is |
| --- | --- |
| `generate_bistable_trajectory.py`, `run_benchmark.py`, `benchmark_pca_vs_topofold.ipynb` | Fully synthetic 60-residue chain built from prescribed (κ, τ); two loop states; noisy termini |
| `run_real_bpti_validation.py` | Crystal structure 5PTI; residues 12–16 rigidly rotated by 36° about the 11–17 axis; per-frame i.i.d. reaction coordinate with 3.5% "transition" frames; noise and random global motions. Not MD |
| `run_honest_real_md_suite.py` | BPTI part: the synthetic BPTI ensemble above (the Zenodo 7347434 cluster PDBs are downloaded but not used in the analysis). Mpro part: 16 representative structures from Zenodo 13730633 joined by linear Cartesian interpolation plus 0.12 Å noise; one chain only |
| `benchmark_topofold_vs_tica.py` | Same synthetic BPTI ensemble; frames have no kinetics, so TICA results are not meaningful |
| `test_blind_detector_bpti.py` | Same synthetic BPTI ensemble |
| `benchmark_abl_kinase_dfg.py`, `benchmark_allosteric_network.py` | Linear interpolation of residues 380–385 between 2GQG and 1IEP; noise; random N-lobe rotations; labels by frame index |
| `benchmark_kras_mpro.py` | Reference structures plus prescribed switch motions and noise |
| `benchmark_idp_alphasynuclein.py` | Random torsion chain; a helical pattern written into residues 67–78 in 22% of frames |
| `benchmark_rna_riboswitch.py` | Crystal structure 1Y26 plus noise and a prescribed switch of the P1 region |
| `benchmark_protac_ternary.py` | **Invalid** (C5): wrong PDB entry 5T3E, second complex built from 5T35, coupling and statistics hard-coded |
| `benchmark_rotamer_gating.py` | **Physically invalid** (C4): Cβ rotated around a fixed Cα on a rigid ideal helix |
| `benchmark_metamorphic_xcl1.py` | Two real NMR structures (1J9O, 2JP1) compared statically; its "SE(3)-invariant Fréchet" is actually a Cartesian Fréchet distance without alignment (M5) |

## Audit scripts (`audit/`)

Line-by-line NumPy re-implementations used to verify the issues in `docs/KNOWN_ISSUES.md`. They need only NumPy, SciPy and scikit-learn.

| Script | Checks |
| --- | --- |
| `audit/writhe_check.py` | Writhe vs Klenin–Langowski and a numerical Gauss double integral (M1) |
| `audit/bc_check.py` | Bimodality coefficient on unimodal, rare-event and angle-wrapping cases (M2, M3) |
| `audit/barrier_check.py` | The former PMF "barrier" estimator on unimodal data; label circularity (C2, C3) |
| `audit/fair_baseline.py` | Global vs local baselines on the synthetic bistable benchmark (C2) |

Planned: replace these with tests against the compiled library once the fixes are merged.

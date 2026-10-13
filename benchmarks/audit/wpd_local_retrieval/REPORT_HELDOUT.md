# WP-D held-out evaluation (mdCATH): Cα-only local-conformation retrieval

**THIRD-PARTY MD DERIVATIVE (mdCATH, CC BY 4.0).** One authorized evaluation (User message 2026-10-06: approval of Amendments 2–3 scientific design; authorization of research-harness changes, reserve downloads and one complete held-out evaluation). Scope: one complete held-out evaluation of the 30-domain cohort; no change to methods, weights, endpoints, margins or selection. 30 domains, 4 workers, wall 7437.6 s, finished 2026-10-06T07:28:00.391866+00:00. Authorization SHA-256 `ef48f6b6563f1bf2…`; binding extension `732909cc99ec8883…`.

**What this result is:** generalization of a *geometric* retrieval benchmark to held-out CATH superfamilies within mdCATH (CHARMM22*, TIP3P, 320 K and 348 K; cross-replica retrieval from a shared starting configuration). **What it is not:** biological, functional or commercial validation; evidence about other force fields, temperatures, multi-domain proteins or independent starting conditions.

## 1. Cohort

| Order | Domain | Role | Superfamily | Replaces |
|---|---|---|---|---|
| 1 | 2wp0D00 | primary | 3.30.300.180 |  |
| 2 | 2wscN00 | primary | 4.10.1190.10 |  |
| 3 | 1neiA00 | primary | 3.30.160.220 |  |
| 5 | 1whoA00 | primary | 2.60.40.760 |  |
| 6 | 3s0yA01 | primary | 3.30.1330.60 |  |
| 7 | 1v8cA01 | primary | 3.10.20.30 |  |
| 8 | 2yx1A01 | primary | 3.30.70.2580 |  |
| 9 | 4ea9A02 | primary | 2.160.10.10 |  |
| 10 | 2p08B00 | primary | 3.30.450.260 |  |
| 11 | 1zvuA03 | primary | 1.10.268.10 |  |
| 12 | 1d8bA00 | primary | 1.10.150.80 |  |
| 13 | 1vz0A01 | primary | 3.90.1530.30 |  |
| 14 | 3wx4A00 | primary | 3.30.70.2770 |  |
| 15 | 2h85A02 | primary | 3.40.50.11580 |  |
| 16 | 2kzsA00 | primary | 1.10.8.810 |  |
| 17 | 4u49B02 | primary | 2.160.20.10 |  |
| 19 | 3udcA02 | primary | 2.30.30.60 |  |
| 20 | 1aorA01 | primary | 3.60.9.10 |  |
| 21 | 4dkwA00 | primary | 3.30.420.280 |  |
| 22 | 1u6zA04 | primary | 3.30.70.2260 |  |
| 23 | 2bt9C01 | primary | 2.40.128.190 |  |
| 25 | 3u61D03 | primary | 1.20.272.10 |  |
| 29 | 3tekA00 | primary | 3.30.470.50 |  |
| 30 | 3nwsA02 | primary | 1.10.10.2380 |  |
| 31 | 3av0A01 | reserve | 3.60.21.10 | 2hb0A01 |
| 32 | 2l23A00 | reserve | 2.40.290.30 | 1s14B00 |
| 33 | 2zj8A03 | reserve | 1.10.3380.30 | 2yztA00 |
| 34 | 3lppA02 | reserve | 2.60.40.1760 | 3slpA00 |
| 35 | 4uuwB02 | reserve | 3.30.70.2860 | 4uhvA03 |
| 36 | 2ljuA01 | reserve | 3.30.160.190 | 3zxsA04 |

Excluded primaries (structural eligibility, original §4.3): 2hb0A01: 320/0: numFrames 130 < 200; 1s14B00: backbone residues not one consecutive run; 2yztA00: 348/0: numFrames 130 < 200; 3slpA00: 348/1: numFrames 130 < 200; 4uhvA03: 348/1: numFrames 130 < 200; 3zxsA04: 348/0: numFrames 120 < 200.

**Mandatory checks passed:** binding vs NumPy max difference 3.8e-14 rad; Horn vs Kabsch max 2.9e-13 Å; binding failures (frames excluded) 0; input integrity, extension identity and production sources verified.

## 2a. T_φψ (primary, angular)

NR@10: 0 = target-optimal, 1 = random. G = 1 − median NR (gain over random). L = relative loss vs the best method (Amendment 1). Intervals: 95 % percentile domain bootstrap (2,000), best re-selected per resample. Classes vs ρ = 0.10. C1 uses its coverage-matched population.

**320/W5/phipsi** — 30 defined domains; observed best: **A**; best counts across resamples: {'B2': 650, 'B1': 5, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 1345}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.817 | 0.183 [0.168, 0.210] | 0.463 [0.413, 0.499] | beyond_margin | informative |
| M2 τ circular-RMS | 0.816 | 0.184 [0.168, 0.205] | 0.463 [0.417, 0.499] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.691 | 0.309 [0.272, 0.341] | 0.096 [0.040, 0.158] | inconclusive | informative |
| B2 window Cα RMSD | 0.659 | 0.341 [0.308, 0.374] | 0.002 [0.000, 0.025] | within_margin | informative |
| B3 κ-only | 0.694 | 0.306 [0.258, 0.342] | 0.105 [0.066, 0.178] | inconclusive | informative |
| A residue-aligned κ–τ | 0.658 | 0.342 [0.302, 0.377] | 0.000 [0.000, 0.022] | within_margin | informative |
| C1 scalar circular mean (coverage-matched) | 0.903 | 0.097 [0.084, 0.112] | 0.720 [0.695, 0.745] | beyond_margin | inconclusive |

**320/W8/phipsi** — 30 defined domains; observed best: **A**; best counts across resamples: {'B2': 0, 'B1': 0, 'B3': 11, 'M1': 0, 'M2': 0, 'A': 1989}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.701 | 0.299 [0.263, 0.334] | 0.270 [0.235, 0.343] | beyond_margin | informative |
| M2 τ circular-RMS | 0.701 | 0.299 [0.263, 0.335] | 0.270 [0.235, 0.346] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.672 | 0.328 [0.309, 0.371] | 0.198 [0.163, 0.250] | beyond_margin | informative |
| B2 window Cα RMSD | 0.637 | 0.363 [0.329, 0.399] | 0.114 [0.075, 0.180] | inconclusive | informative |
| B3 κ-only | 0.631 | 0.369 [0.328, 0.425] | 0.099 [0.018, 0.170] | inconclusive | informative |
| A residue-aligned κ–τ | 0.591 | 0.409 [0.379, 0.472] | 0.000 [0.000, 0.000] | within_margin | informative |
| C1 scalar circular mean (coverage-matched) | 0.913 | 0.087 [0.069, 0.096] | 0.788 [0.771, 0.825] | beyond_margin | uninformative |

**348/W5/phipsi** — 30 defined domains; observed best: **B2**; best counts across resamples: {'B2': 1023, 'B1': 0, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 977}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.798 | 0.202 [0.173, 0.254] | 0.467 [0.369, 0.515] | beyond_margin | informative |
| M2 τ circular-RMS | 0.798 | 0.202 [0.173, 0.255] | 0.466 [0.369, 0.516] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.653 | 0.347 [0.301, 0.380] | 0.085 [0.033, 0.130] | inconclusive | informative |
| B2 window Cα RMSD | 0.621 | 0.379 [0.338, 0.417] | 0.000 [0.000, 0.020] | within_margin | informative |
| B3 κ-only | 0.682 | 0.318 [0.284, 0.350] | 0.162 [0.101, 0.233] | beyond_margin | informative |
| A residue-aligned κ–τ | 0.621 | 0.379 [0.337, 0.415] | 0.001 [0.000, 0.030] | within_margin | informative |
| C1 scalar circular mean (coverage-matched) | 0.883 | 0.117 [0.081, 0.140] | 0.691 [0.645, 0.765] | beyond_margin | inconclusive |

**348/W8/phipsi** — 30 defined domains; observed best: **A**; best counts across resamples: {'B2': 0, 'B1': 0, 'B3': 68, 'M1': 0, 'M2': 0, 'A': 1932}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.694 | 0.306 [0.276, 0.374] | 0.328 [0.206, 0.369] | beyond_margin | informative |
| M2 τ circular-RMS | 0.693 | 0.307 [0.276, 0.365] | 0.326 [0.209, 0.368] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.632 | 0.368 [0.332, 0.401] | 0.194 [0.135, 0.231] | beyond_margin | informative |
| B2 window Cα RMSD | 0.603 | 0.397 [0.365, 0.426] | 0.129 [0.079, 0.180] | inconclusive | informative |
| B3 κ-only | 0.567 | 0.433 [0.382, 0.463] | 0.050 [0.000, 0.142] | inconclusive | informative |
| A residue-aligned κ–τ | 0.544 | 0.456 [0.420, 0.477] | 0.000 [0.000, 0.006] | within_margin | informative |
| C1 scalar circular mean (coverage-matched) | 0.912 | 0.088 [0.076, 0.114] | 0.806 [0.758, 0.828] | beyond_margin | inconclusive |

**Overall outcome per method (Amendment 1 §4 precedence):** M1: `fails_criterion_in_all_cells`; M2: `fails_criterion_in_all_cells`; B1: `inconclusive`; B2: `inconclusive`; B3: `inconclusive`; A: `meets_criterion_in_all_cells`; C1: `fails_criterion_in_all_cells`.

## 2b. T_pp (separate, Cartesian peptide-plane)

NR@10: 0 = target-optimal, 1 = random. G = 1 − median NR (gain over random). L = relative loss vs the best method (Amendment 1). Intervals: 95 % percentile domain bootstrap (2,000), best re-selected per resample. Classes vs ρ = 0.10. C1 uses its coverage-matched population.

**320/W5/pp** — 30 defined domains; observed best: **B2**; best counts across resamples: {'B2': 2000, 'B1': 0, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 0}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.673 | 0.327 [0.297, 0.367] | 0.446 [0.396, 0.473] | beyond_margin | informative |
| M2 τ circular-RMS | 0.673 | 0.327 [0.297, 0.374] | 0.446 [0.386, 0.473] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.480 | 0.520 [0.499, 0.564] | 0.118 [0.061, 0.129] | inconclusive | informative |
| B2 window Cα RMSD | 0.411 | 0.589 [0.539, 0.618] | 0.000 [0.000, 0.000] | within_margin | informative |
| B3 κ-only | 0.645 | 0.355 [0.335, 0.377] | 0.398 [0.355, 0.428] | beyond_margin | informative |
| A residue-aligned κ–τ | 0.461 | 0.539 [0.504, 0.571] | 0.086 [0.068, 0.105] | inconclusive | informative |
| C1 scalar circular mean (coverage-matched) | 0.825 | 0.175 [0.153, 0.228] | 0.703 [0.620, 0.732] | beyond_margin | informative |

**320/W8/pp** — 30 defined domains; observed best: **B2**; best counts across resamples: {'B2': 2000, 'B1': 0, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 0}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.475 | 0.525 [0.489, 0.570] | 0.340 [0.293, 0.376] | beyond_margin | informative |
| M2 τ circular-RMS | 0.478 | 0.522 [0.489, 0.568] | 0.343 [0.294, 0.376] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.295 | 0.705 [0.679, 0.726] | 0.113 [0.093, 0.143] | inconclusive | informative |
| B2 window Cα RMSD | 0.206 | 0.794 [0.777, 0.828] | 0.000 [0.000, 0.000] | within_margin | informative |
| B3 κ-only | 0.586 | 0.414 [0.375, 0.442] | 0.479 [0.441, 0.538] | beyond_margin | informative |
| A residue-aligned κ–τ | 0.363 | 0.637 [0.612, 0.658] | 0.198 [0.179, 0.230] | beyond_margin | informative |
| C1 scalar circular mean (coverage-matched) | 0.856 | 0.144 [0.134, 0.165] | 0.819 [0.795, 0.830] | beyond_margin | informative |

**348/W5/pp** — 30 defined domains; observed best: **B2**; best counts across resamples: {'B2': 2000, 'B1': 0, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 0}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.643 | 0.357 [0.301, 0.410] | 0.424 [0.363, 0.480] | beyond_margin | informative |
| M2 τ circular-RMS | 0.643 | 0.357 [0.301, 0.410] | 0.423 [0.362, 0.480] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.457 | 0.543 [0.512, 0.590] | 0.121 [0.064, 0.136] | inconclusive | informative |
| B2 window Cα RMSD | 0.381 | 0.619 [0.569, 0.644] | 0.000 [0.000, 0.000] | within_margin | informative |
| B3 κ-only | 0.627 | 0.373 [0.352, 0.390] | 0.397 [0.348, 0.432] | beyond_margin | informative |
| A residue-aligned κ–τ | 0.425 | 0.575 [0.517, 0.602] | 0.071 [0.055, 0.094] | within_margin | informative |
| C1 scalar circular mean (coverage-matched) | 0.796 | 0.204 [0.161, 0.233] | 0.671 [0.631, 0.725] | beyond_margin | informative |

**348/W8/pp** — 30 defined domains; observed best: **B2**; best counts across resamples: {'B2': 2000, 'B1': 0, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 0}; C1 coverage-matched domains: 30.

| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |
|---|---|---|---|---|---|
| M1 τ chord-RMS | 0.441 | 0.559 [0.475, 0.603] | 0.321 [0.278, 0.398] | beyond_margin | informative |
| M2 τ circular-RMS | 0.440 | 0.560 [0.472, 0.593] | 0.320 [0.286, 0.400] | beyond_margin | informative |
| B1 Cα–Cα distances | 0.270 | 0.730 [0.683, 0.747] | 0.114 [0.089, 0.164] | inconclusive | informative |
| B2 window Cα RMSD | 0.176 | 0.824 [0.793, 0.846] | 0.000 [0.000, 0.000] | within_margin | informative |
| B3 κ-only | 0.583 | 0.417 [0.371, 0.441] | 0.494 [0.452, 0.549] | beyond_margin | informative |
| A residue-aligned κ–τ | 0.334 | 0.666 [0.620, 0.678] | 0.192 [0.168, 0.237] | beyond_margin | informative |
| C1 scalar circular mean (coverage-matched) | 0.841 | 0.159 [0.135, 0.176] | 0.807 [0.785, 0.831] | beyond_margin | informative |

**Overall outcome per method (Amendment 1 §4 precedence):** M1: `fails_criterion_in_all_cells`; M2: `fails_criterion_in_all_cells`; B1: `inconclusive`; B2: `meets_criterion_in_all_cells`; B3: `fails_criterion_in_all_cells`; A: `cell_dependent`; C1: `fails_criterion_in_all_cells`.

## 3. Secondary (reported, not deciding)

| Cell | defined / total domains | windows defined | undefined windows | random overlap@10 | M1 ovl@10 / NR@1 / ρ | M2 ovl@10 / NR@1 / ρ | B1 ovl@10 / NR@1 / ρ | B2 ovl@10 / NR@1 / ρ | B3 ovl@10 / NR@1 / ρ | A ovl@10 / NR@1 / ρ | C1 ovl@10 / NR@1 / ρ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 320/W5/phipsi | 30/30 | 3490 | {} | 0.0222 | 0.0598 / 0.836 / 0.310 | 0.0596 / 0.836 / 0.305 | 0.0785 / 0.706 / 0.369 | 0.0948 / 0.674 / 0.407 | 0.0781 / 0.726 / 0.427 | 0.0994 / 0.666 / 0.393 | 0.0399 / 0.915 / 0.196 |
| 320/W5/pp | 30/30 | 3490 | {} | 0.0222 | 0.0935 / 0.719 / 0.474 | 0.0935 / 0.719 / 0.475 | 0.1667 / 0.501 / 0.594 | 0.2013 / 0.444 / 0.675 | 0.1010 / 0.674 / 0.450 | 0.1658 / 0.479 / 0.581 | 0.0571 / 0.847 / 0.329 |
| 320/W8/phipsi | 30/30 | 3400 | {} | 0.0222 | 0.0907 / 0.727 / 0.399 | 0.0897 / 0.727 / 0.399 | 0.0884 / 0.658 / 0.346 | 0.1012 / 0.639 / 0.369 | 0.1013 / 0.653 / 0.471 | 0.1229 / 0.603 / 0.465 | 0.0369 / 0.925 / 0.179 |
| 320/W8/pp | 30/30 | 3400 | {} | 0.0222 | 0.1664 / 0.501 / 0.568 | 0.1664 / 0.502 / 0.565 | 0.2776 / 0.299 / 0.727 | 0.3778 / 0.225 / 0.851 | 0.1229 / 0.604 / 0.457 | 0.2267 / 0.343 / 0.617 | 0.0481 / 0.867 / 0.279 |
| 348/W5/phipsi | 30/30 | 3490 | {} | 0.0222 | 0.0687 / 0.819 / 0.355 | 0.0686 / 0.819 / 0.353 | 0.0853 / 0.659 / 0.401 | 0.0989 / 0.634 / 0.433 | 0.0812 / 0.717 / 0.463 | 0.1031 / 0.631 / 0.436 | 0.0459 / 0.897 / 0.233 |
| 348/W5/pp | 30/30 | 3490 | {} | 0.0222 | 0.1028 / 0.690 / 0.534 | 0.1027 / 0.690 / 0.528 | 0.1778 / 0.467 / 0.617 | 0.2194 / 0.407 / 0.718 | 0.1033 / 0.657 / 0.471 | 0.1844 / 0.437 / 0.632 | 0.0619 / 0.818 / 0.353 |
| 348/W8/phipsi | 30/30 | 3400 | {} | 0.0222 | 0.0993 / 0.692 / 0.447 | 0.0975 / 0.690 / 0.444 | 0.0974 / 0.624 / 0.382 | 0.1120 / 0.599 / 0.411 | 0.1089 / 0.603 / 0.513 | 0.1326 / 0.556 / 0.503 | 0.0408 / 0.919 / 0.209 |
| 348/W8/pp | 30/30 | 3400 | {} | 0.0222 | 0.1720 / 0.460 / 0.577 | 0.1731 / 0.460 / 0.576 | 0.2907 / 0.278 / 0.750 | 0.4073 / 0.196 / 0.876 | 0.1253 / 0.606 / 0.493 | 0.2423 / 0.318 / 0.624 | 0.0522 / 0.861 / 0.318 |

Target-spread terciles (windows pooled within cell; mean window NR@10 per tercile T1 / T2 / T3):

| Cell | bounds | M1 | M2 | B1 | B2 | B3 | A | C1 |
|---|---|---|---|---|---|---|---|---|
| 320/W5/phipsi | 0.188–0.286; 0.286–0.419; 0.420–2.032 | 0.880 / 0.784 / 0.630 | 0.880 / 0.784 / 0.631 | 0.807 / 0.657 / 0.496 | 0.783 / 0.621 / 0.474 | 0.827 / 0.675 / 0.530 | 0.776 / 0.617 / 0.489 | 0.941 / 0.872 / 0.765 |
| 320/W5/pp | 0.205–0.340; 0.340–0.510; 0.510–2.519 | 0.765 / 0.674 / 0.473 | 0.765 / 0.675 / 0.473 | 0.565 / 0.477 / 0.349 | 0.528 / 0.408 / 0.282 | 0.718 / 0.614 / 0.565 | 0.571 / 0.449 / 0.335 | 0.860 / 0.815 / 0.661 |
| 320/W8/phipsi | 0.198–0.311; 0.311–0.502; 0.503–1.765 | 0.803 / 0.673 / 0.552 | 0.803 / 0.673 / 0.558 | 0.760 / 0.599 / 0.543 | 0.743 / 0.571 / 0.527 | 0.762 / 0.578 / 0.484 | 0.708 / 0.536 / 0.475 | 0.940 / 0.885 / 0.851 |
| 320/W8/pp | 0.274–0.495; 0.495–0.830; 0.830–3.421 | 0.580 / 0.493 / 0.401 | 0.580 / 0.493 / 0.402 | 0.366 / 0.302 / 0.254 | 0.289 / 0.206 / 0.136 | 0.626 / 0.537 / 0.596 | 0.431 / 0.345 / 0.343 | 0.864 / 0.833 / 0.794 |
| 348/W5/phipsi | 0.195–0.312; 0.312–0.522; 0.522–1.985 | 0.873 / 0.765 / 0.601 | 0.873 / 0.765 / 0.601 | 0.798 / 0.626 / 0.462 | 0.774 / 0.590 / 0.447 | 0.820 / 0.649 / 0.515 | 0.766 / 0.588 / 0.467 | 0.937 / 0.861 / 0.737 |
| 348/W5/pp | 0.213–0.368; 0.368–0.633; 0.634–2.200 | 0.758 / 0.655 / 0.432 | 0.758 / 0.655 / 0.432 | 0.560 / 0.456 / 0.311 | 0.520 / 0.388 / 0.251 | 0.713 / 0.604 / 0.545 | 0.563 / 0.429 / 0.309 | 0.856 / 0.801 / 0.612 |
| 348/W8/phipsi | 0.211–0.339; 0.339–0.753; 0.755–1.861 | 0.794 / 0.649 / 0.526 | 0.794 / 0.649 / 0.532 | 0.748 / 0.568 / 0.519 | 0.731 / 0.542 / 0.506 | 0.749 / 0.553 / 0.464 | 0.696 / 0.510 / 0.461 | 0.936 / 0.869 / 0.830 |
| 348/W8/pp | 0.289–0.547; 0.547–1.021; 1.021–3.409 | 0.579 / 0.472 / 0.365 | 0.579 / 0.473 / 0.368 | 0.366 / 0.293 / 0.233 | 0.286 / 0.194 / 0.116 | 0.621 / 0.528 / 0.584 | 0.428 / 0.333 / 0.323 | 0.862 / 0.817 / 0.767 |

**Within-domain dependence diagnostics** (descriptive; AR(1)-based heuristic only, Amendment 3 §4): class counts over all domain × cell × pair series: {'estimated_below_heuristic_cutoff': 2153, 'above_heuristic_cutoff': 1447}. They do not affect the domain-level inference.

## 4. Interpretation limits

- Intervals are marginal and unadjusted (7 methods × 4 cells × 2 statistics per target); no family-wise coverage. Classes are pre-declared summaries, not tests; "within the margin" is not equivalence.
- ρ = 0.10 is a stated value judgment, not derived from a downstream task.
- Replicas share their starting configuration; cross-replica retrieval may be easier than between independent runs. The first 10 % of frames were excluded.
- Method A entered the candidate set after AdK development results; this held-out evaluation is its first test.
- Both targets are geometric reconstructions from the same trajectories; T_φψ may favour angular and T_pp Cartesian methods.
- 30 domains, one per superfamily, from a seeded draw; generalization is to that population within mdCATH.


# Correction Note: descriptive comparison tool v0.2.0 (2026-10-23)

**Scope:** three numerical issues in `python/topofold_inputs/ensemble.py`. The v1 demonstrations (`results/`) are
kept unchanged as history. Corrected demonstrations are in `results_v2/`, verified by `verify_v2.py`
(`results_v2/verification.json`). Jensen–Shannon calculations are unchanged.

| Issue in v1 | Correction in v2 | Effect on the demonstrations |
|---|---|---|
| The circular temporal diagnostic averaged separately normalized cos θ and sin θ autocorrelations, which depends on the angular origin | C(k) = Re Σ conj(u_t) u_{t+k} / Σ \|u_t\|² with u_t = e^{iθ_t} − mean(e^{iθ}); C(0) = 1; invariant to a common rotation (tested at five shifts) | τ lag-1 values changed by up to 0.27 (1ha8A00) and 0.32 (AdK); τ decorrelation summaries changed. Linear channels unchanged |
| Circular summaries had no resultant length; a near-zero resultant gave an arbitrary mean direction, and rounding could give R > 1 | The resultant length R is reported and clipped to [0, 1]. When R < 10⁻⁶ (a numerical conditioning rule, not a biological criterion) the center, the spread and any center difference involving it are reported as unavailable, with the reason | None on these data (all statuses "ok") |
| Single-frame, constant or short series could produce NaN or be reported as "beyond max lag" | Explicit statuses: "fewer than 2 frames", "constant series" (sample variance ≤ 10⁻²⁰), and "not reached within L lags", where L is the number of lags actually available | None on these data; covered by SYNTHETIC tests |

Verified in `results_v2/verification.json`:
- pooled JS distances are identical to v1 at the recorded precision (6 significant digits);
- the corrected τ diagnostic is origin-invariant on 902 real τ series (max |Δ| 1.0 × 10⁻¹⁵);
- unavailable values and reasons reach the tables and the report (SYNTHETIC tests);
- the manifests identify the module (v0.2.0, matching SHA-256) and the native extension loaded (`d90e7ed4…`).

**Pooling (unchanged, now documented):** a group's pooled distribution concatenates its members' retained frames, so
members are weighted by their frame counts (`member_frames` column; `pooling` in the manifest).

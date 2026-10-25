# structpair advisory correction (2026-10-25): failing before, passing after

**The defect (structpair 0.2.0, schema /2).** The default advisory output flagged a position as `separated` when the between-entry value
equalled the within-entry value at zero (identical chains), and when differences were below the numerical tolerance (the comparison
was widened by 1e-9). NaN and infinity were accepted as the within factor and as minimum differences.

**The correction (structpair 0.3.0, schema /3).** The numbers are always reported; the flag is assessed only for a channel with a
configured finite minimum difference above 1e-9; equality at that boundary meets the rule, with no tolerance widening; values within
1e-9 of the threshold are marked; non-finite numbers are refused at spec parsing. Full rule: `docs/STRUCTPAIR_CHAINS_SPEC_2026-10-25.md`
(amendment 2). Records of schema /2 keep their historical meaning; the viewer reads both.

| File | Content |
|---|---|
| `record_regression.py` | Runs `tests/test_structpair_advisory.py` against the committed 0.2.0 code (temporary git worktree) and against the corrected code |
| `results/regression_before_after.json` | The record: 19 cases, **13 fail on 0.2.0 and pass on 0.3.0**, 6 pass on both (non-finite trim fraction, negative infinity, and identical chains once a floor is configured were already handled), 0 fail after |

Numerical compatibility: the numeric columns do not change between 0.2.0 and 0.3.0. The recorded K1 comparison (1,224 values, 0 mismatches at the
declared 1e-4 Å tolerance, flagged set equal to the recorded K1 set, frozen K1 gate result unchanged) is repeated in
`../structpair_chains_demo/results/v0_3_0/verification.json`. In the 0.2.0 K1 record the flag was true at 78 position-channel cells
(8 on the displacement channel, where K1's 1 Å floor was configured, 70 on channels without a floor); in the 0.3.0 record it is true at
the same 8 displacement cells and nowhere else.

```bash
.venv/bin/python benchmarks/audit/structpair_advisory_correction_2026-10-25/record_regression.py
```

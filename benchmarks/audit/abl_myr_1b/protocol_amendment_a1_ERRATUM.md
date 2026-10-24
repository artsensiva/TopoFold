# Erratum to protocol amendment a1 (dated 2026-10-25, written before any trajectory was opened)

`protocol_amendment_a1.json` is frozen and unchanged (its hash is in `results/ABL-MYR-1b_protocol_amendment_a1_freeze.json`).
Its change C1 says the bilateral G5 stop happens "before any extraction or analysis". That wording is inaccurate for the
unchanged frozen runner, which:

1. extracts the 40 trajectory members from `PCA.zip` at G1 (a file operation, no frame is read for analysis),
2. reads only their time stamps and atom counts at G2 and G5.

**Correct reading of C1:** if any of the four groups has fewer than 8 usable runs after the G5 exclusions, the run stops at G5
**before any frame is analysed**: no κ/τ, distance, RMSF or φ/ψ values are computed, no statistic is formed, and no outcome is
produced. The trajectory files will already have been extracted into the git-ignored `data_dryad/extracted/`, and their time
stamps read. The rule, the thresholds (450 frames, 8 runs) and the scientific question are unchanged. The wrapper tests check the
stop at G5 and the absence of any result.

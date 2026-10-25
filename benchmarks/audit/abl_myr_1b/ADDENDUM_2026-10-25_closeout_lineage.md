# Addendum (2026-10-25): lineage of the ABL-MYR-1b closeout and its later supplement

This addendum is a living note. It does not change the closeout, the closeout ledger or any frozen file, and the ledger
does not list it.

**ABL-MYR-1b is closed. The comparison result is final for its protocol.**

| Item | Commit | Status |
|---|---|---|
| Original closeout snapshot: `CLOSEOUT_2026-10-25.md`, `results/ABL-MYR-1b_closeout_inventory.json` (SHA-256 `3506ce36d134da4189358d4e8ddb51e900e81e3dc5da4fc90c9713b10b7c19f0`) and `results/ABL-MYR-1b_closeout_SHA256SUMS.txt` (75 lines, SHA-256 `208e42b578fadd1469d200502411dabb63e2ca85c0b947c421b69efb16aa5d1f`) | `35c699f` | **the original closeout snapshot** |
| Residue-coverage figure `figures/ABL-MYR-1b_residue_coverage_SUPPLEMENT_2026-10-25.png` with caption and provenance | `2b75a0d` | **a later supplement**, not part of the closeout; it only reads recorded results and changes no frozen file |

**What happened in `2b75a0d`, stated plainly.** When the supplement was added, the closeout inventory and checksum list
were regenerated (inventory SHA-256 `94d8431f35ce2016217737cba845d999fff9b8734d35f33ef751001df441b077`; list of 80 lines,
SHA-256 `7aad7e3ff08e95b51f70ab2514a583652b096900ced701c4309e0f0cdd552b08`), because the supplement edit of
`README.md` and the new files had changed the files the original list covers. That rewrote a closeout ledger, which should
not have been done. The original files are unchanged in Git history.

**How to verify each.**
- Original snapshot: `git show 35c699f:benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_closeout_SHA256SUMS.txt` against the
  files of that commit. Against today's working tree, the original 75-line list matches every file except
  `benchmarks/audit/abl_myr_1b/README.md`, a status document that gained the supplement section in `2b75a0d`.
- Ledger regenerated in `2b75a0d`: `sha256sum -c benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_closeout_SHA256SUMS.txt`
  from the repository root verifies all 80 lines at the time of this addendum.

**Rule from here.** The closeout inventory and checksum list are not regenerated again. Changes to living documentation
after the closeout are recorded as dated addenda like this one, not by rewriting the ledger or the files it lists. The
backup gap of the closeout is unchanged: the hashes and the `~/Downloads` copy are not an independent backup, and nothing
was uploaded.

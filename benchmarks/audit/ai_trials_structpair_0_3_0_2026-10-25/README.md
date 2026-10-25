# Two independent AI technical trials of structpair 0.3.0 (commit 19b5268, 2026-10-25): review and close-out

**Independent AI technical trials, not human usability studies, not human validation and not biological validation.** Both testers
cloned commit `19b5268180d533efedcbdef846e7b3f24250a767`, installed it, ran the KRAS smoke test and an Abl comparison of their own
design, used the viewer in a real browser and tampered with a copy of a result. No project developer helped either tester.

| | Trial 1 | Trial 2 |
|---|---|---|
| Tester | GPT-6 family in Codex desktop (exact variant not exposed) | Claude Opus 5.5 (`claude-opus-5-5`), Claude desktop app |
| Original report | `/home/artem/Documents/Codex/2026-10-25/files-pasted-by-the-user-you/outputs/topofold-independent-trial/REPORT.md` | `/home/artem/.config/Claude/scratch-workspaces/e9fd6719-…/scratch-2026-10-25-04801e/trial/REPORT.md` |
| Abl range (predeclared before the run) | author residues 248–534, 287 positions; missing tails kept and excluded (`on_missing=exclude`); no advisory minimum | 248–529, 282 positions; `on_missing=reject`; advisory `min_difference` `{disp: 1.0 Å}` |
| Time | 6.62 h wall (interrupted), about 28 min active | about 20 min |

**The 282- and 287-position outputs are not a numerical disagreement.** The two testers chose different residue ranges and policies
before they saw any result (Trial 1 kept the unobserved tail 530–534 mapped; Trial 2 stopped at 529). Where they overlap they agree: both
report 5.121 Å global Cα RMSD over the same 282 Cα, the same dominant region (residues about 398–419, a displacement of about 33 Å at
409 after the trimmed fit, against about 1 Å between the two 1OPJ copies), and the same confounds.

## Evidence matrix (outcomes of each tester, and what this review did)

| Area | Trial 1 | Trial 2 | Review of this archive |
|---|---|---|---|
| Installation | Completed following the documentation (fresh clone, pinned packages, build succeeded) | Completed following the documentation, no failures | Consistent. Downloader writes into the clone (documented, guide §2) |
| Computation | Completed with a self-devised workaround (relocated downloader); KRAS 158 positions, 0.974 Å; Abl once, no reruns; 21 of 21 hash checks | Completed following the documentation; KRAS reproduced; Abl first attempt; 9 of 9 manifest hashes; self-written checker | Guide stdout and the byte-identical `summary.json` statement were out of date: corrected (guide §2) |
| Viewer use | Completed with a workaround (screenshot helper; server restarted after a usage-limit interruption). Settled states consistent in the checked paths; **no persistent defect found** | Completed with a workaround (port 8501 busy). **One reproducible defect: sidebar focus after a Results-folder change and back.** 3D view not tested (no WebGL in its pane) | **Defect reproduced on 19b5268 and fixed** (below). 3D highlight checked visually in a WebGL-capable headless Chrome |
| Integrity checking | One altered value in a copied `windows.csv`; the viewer refused it, naming the file; originals restored | One altered value in a copy; refused, naming the file; original passes 9 of 9 | Consistent. A documented manifest-check command now exists (`structpair_verify`, guide §3b) |
| Interpretation | Frozen before viewer testing and before any shipped Abl result was opened; attribution to ligand, mutation or resistance explicitly refused | Same, written without opening the shipped demonstration | Consistent with each other and with the shipped demonstration; no new scientific claim made here |
| Limitations stated | Not human; not an audit of every metric; tester mistakes and environment failures listed and not blamed on the tool | Not human; 3D not testable; one platform | Both: one machine, one OS, AI agents, no human usability evidence |

## Findings: reproduced, fixed, documented, open

| ID | Finding | Source | Status |
|---|---|---|---|
| F1 | Sidebar focus shows position 0 while details, chart marker, table row and 3D markers keep the old window, after Results folder → a folder whose result fails its check → back | Trial 2 | **Reproduced** on 19b5268 (3 of 3 browser runs; `../structpair_viewer_folder_switch_2026-10-25/`). **Fixed**; regression test fails on the old code, passes now. A change to another *valid* folder and back did **not** show the defect |
| F2 | Guide: stdout lacked the within-entry clause | Trial 2 | **Fixed** (guide §2) |
| F3 | Guide: `summary.json` "byte-identical" | Trial 2 | **Fixed**: the three CSV files are byte-identical, `summary.json` adds the `within_entry_variation` block (verified against Trial 2's output) |
| F4 | Guide: manifest schema `/2` in §3a | Trial 2 | **Fixed** (`/3` since 0.3.0, `/2` for 0.2.0 results) |
| F5 | No documented manifest check | Trial 2 | **Fixed**: `python -m topofold_inputs.structpair_verify` and guide §3b |
| F6 | 3D view needs WebGL; port 8501 busy | Trial 2 | **Documented** (guide §6) and a one-sentence note in the viewer |
| F7 | Downloader has no output option | Trial 1 | **Documented** (guide §2); no code change |
| F8 | Exclusion `detail` text cites row indices (58 and 59) while `location` gives A:59→A:68 | Trial 2 | **Open**: the text comes from the frozen exclusion records; changing it alters result files |
| F9 | Explanation column truncated in the grid; wide tables clipped | Both | **Open** (interface redesign, kept separate) |
| F10 | Window-position axes are indices, not residue numbers; hover needed | Both | **Open** (redesign) |
| F11 | `chain_variation.png` axis text says one tick label per position; ticks are sparse | Trial 2 | **Open**: changing it alters regenerated figures |
| F12 | `A`/`B` mean entries and also author chain IDs; table `A_residue` is the window centre; 38 exclusion rows are causes of 10 windows; trimmed displacement and local RMSD use different fits | Trial 1 | **Open** (wording; guide already states some of these) |
| F13 | Transient mixed states (1–3 s) after a control change | Both | **Unchanged**: expected Streamlit behaviour, guide §6; every state settled consistently in the fixed code |
| F14 | Trial 1's own folder round trip restored focus 162 in the sidebar (its snapshot 17) | Trial 1 | **Unreconciled**: its record does not isolate the conditions (for example a page reload in between); the defect needs a run that stops before the focus widget is drawn |

## Archive

- `inventory.json`: SHA-256, size and original path of every file of both trials outside their environments, clones, build output and `work/` folders (Trial 1: 115 files, Trial 2: 80 files, including all screenshots, results, corrupted copies and inputs).
- `archive/trial1/`, `archive/trial2/`: compact copies (reports, specs, manifests and summaries of every run, pre-registrations, interpretations, evidence records, captions, small logs, the five screenshots of the viewer finding). The original screenshots remain in place and are hashed in the inventory.
- `build_archive.py --verify-only` re-hashes the archive and the originals. At the time of this review: 115 of 115 and 80 of 80 originals unchanged, 17 and 28 archive copies match the inventory.
- The originals were read only. The trial environments, clones and the Abl demonstration of this repository were not copied.

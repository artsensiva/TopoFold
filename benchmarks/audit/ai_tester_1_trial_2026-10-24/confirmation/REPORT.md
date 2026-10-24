# Bounded TopoFold confirmation — AI tester 1

Commit **049a47f196e39ec46cc8e2c9d0bd63f03600df71**, resolved before execution. Fresh checkout and isolated sibling Python environment, Cargo downloads/build target and caches. Reported model family GPT-6; exact model identifier unknown; application Codex desktop (version unknown). No project code edited. No parallel agent assistance used in this confirmation.

**Outcome 1 — setup instructions sufficient:** the §2 dependency/build/example/viewer sequence worked on this Linux environment with existing Python 3.14.4 and Rust 1.98.1. No missing dependency, extra package, changed package version or adapter-python option was needed. This is not a literal unchanged-command transcript: fixed commit checkout, trial-local output/cache paths and automated-launch environment settings were applied and recorded. The added kras1 run populated the second comparison requested by this test; §2 itself generates only kras2. [Every command and deviation](evidence/COMMANDS.md).

The documented kras2 invocation exited 0 in 7.039 seconds, printing the expected `158 positions; global Cα RMSD 0.974 Å over 158 Cα; 10 exclusion rows`. The additional shipped kras1 spec exited 0 in 6.102 seconds. Viewer dependencies installed at the guide's pinned versions and the launch command reached HTTP 200. No broader scientific validation, extra specification or numerical comparison was repeated.

**Outcome 2 — selection synchronization:** the observed settled states were synchronized; the previous persistent wrong-row checkbox after changing the filter was not observed. The strict requirement that a mismatched row must *never* occur is **not fully confirmed**: immediately after each comparison switch, accessibility snapshots briefly combined the new sidebar focus with old details/selected row (84/58/58 when entering kras1;58/84/84 when returning). The following screenshot/reset and settled snapshots synchronized. This evidence cannot distinguish a visible transient from accessibility-cache lag. Therefore this is a qualified settled-state pass, not an unconditional “never” pass.

| Requested step | Sidebar | Details | Checked table row |
|---|---:|---:|---|
| Load kras2 | 79 | 79 | None visible; focus outside table viewport |
| Sort d_kappa descending | 79 | 79 | None visible; focus outside table viewport |
| Filter without values | 79 | 79 | None; focus filtered out; explanatory message shown |
| Select position58 | 58 | 58 | 58 |
| Switch to all | 58 | 58 | None visible; selected position outside first displayed rows; position33 visibly unchecked |
| Sort by position | 58 | 58 | None visible before interruption; row58 confirmed checked after restoring state and scrolling |
| Switch to kras1 (settled) | 84 | 84 | 84 |
| Switch back to kras2 (settled) | 58 | 58 | 58 |

The table is virtualized: “none visible” means the visible rows were unchecked, not that there was no selection offscreen. Under “without values,” the app explicitly stated focus 79 was excluded by the filter and no row was marked. After selection 58, switching to all left row 33 unchecked. For the resumed position sort and comparison switches, the focused rows were scrolled into view and their actual checked boxes inspected. [Step record](evidence/16-step-outcomes.json), [direct resumed accessibility excerpts](evidence/browser/resume-browser-records.json).

Seven screenshots are preserved in [browser evidence](evidence/browser), especially [selected58](evidence/browser/03-selected-58.jpg), [all after58](evidence/browser/04-all-after-58.jpg), [position sort](evidence/browser/05-position-sort-resumed.jpg), [kras1 row84](evidence/browser/06-kras1-focus-84.jpg), and [return to kras2 row58](evidence/browser/07-return-kras2-focus-58.jpg). Pre-interruption step records were reconstructed from actual tool outputs in the conversation and surviving screenshots; the in-memory browser record array did not survive the reboot.

**Interruption and scope:** initial phase started 2026-10-24 08:28:16 UTC. Usage interruption and computer restart occurred after the position-sort check, before the comparison switches. Resumed 12:22:04 UTC; services stopped 12:25:54 UTC. Elapsed 3h 57m 38s; active work approximately 14 minutes, with exact active time unknown because the suspension instant was not recorded. The resumed measured phase was 3m 50s. This was not an uninterrupted trial. Only the viewer state and services needed restoration; completed install and output files were reused within this same fresh confirmation. [Timing](evidence/17-timing.json).

All recorded finite setup/computation commands exited 0. A canvas accessibility click failed with “Node does not have a layout object”; visible-coordinate checkbox selection succeeded. Sandbox read probes emitted stream-fd warnings but returned exit 0; final approved-access Git checks returned clean. [Errors and limitations](evidence/18-errors-and-limitations.json). Final Git status empty, staged/unstaged diffs exit 0, same full commit. [Final commands, versions and status](evidence/15-final-checks.json). Resumed services stopped and observed absent; detached exit statuses unknown. [Shutdown](evidence/14-shutdown.json).

If the strict never-mismatch condition is required, the remaining issue is transition-time synchronization/accessibility behavior; this confirmation stopped after the requested sequence and did not investigate or fix its implementation.

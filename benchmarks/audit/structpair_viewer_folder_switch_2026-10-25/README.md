# Viewer: focus after a Results-folder change and back (AI trial 2 finding F1), 2026-10-25

Real-browser reproduction (headless Chrome 153 through the DevTools protocol, Streamlit 1.64.0, Ubuntu 26.04), then fix and verification.
Script: `apps/structpair_viewer/browser_check_folder.py` (`.venv/bin/python apps/structpair_viewer/browser_check_folder.py URL OUT FOLDER_A FOLDER_B`).
Test: `tests/test_structpair_viewer_browser.py` (skipped without Chrome, `websocket-client` and Streamlit; about 5 minutes).

Sequence, as reported by Trial 2: select comparison kras2; focus 39; change Results folder to a folder whose only result fails its manifest check
(one value of a copied `windows.csv` altered, manifest unchanged); change it back; wait until the page has settled.

| State after the return | Commit 19b5268 (`old_commit_19b5268/`) | Fixed (`fixed/failing_result_and_back/`) |
|---|---|---|
| Sidebar focus | **0 · A:1 ↔ B A:1** | 39 · A:40 ↔ B A:40 |
| Details heading | position 39 | position 39 |
| Selected table row | 39 | 39 |
| Chart marker (x of the orange rule) | 196 px, as when 39 was set | 196 px |
| 3D markers (Plotly trace text) | A:39–A:41 | A:39–A:41 |

On 19b5268 the failure appeared in 3 of 3 runs (two scripted runs and the new test run against the old code). The values are in the JSON records.
A change to another valid folder and back was also run (`fixed/valid_other_folder_and_back/`): no inconsistency, on 19b5268 or after the fix. In that case the
first comparison of the folder is shown after each change, so the comparison itself is not restored; the remembered focus is.

**Cause.** A run that stops before the focus widget is drawn (a result that fails its manifest check) removes the widget from the browser while its
value stays in the server session. The next run saw a valid stored value and did not push it again, so the browser drew the widget's first option.
**Fix** (`apps/structpair_viewer/app.py`): the focus is kept per comparison after every draw, and a run that follows one without the focus widget assigns it again.

**Verification.** The three browser cases pass on the fixed code (3 passed, 287 s); the failing-result case fails on 19b5268 (sidebar 0 against 39) and the valid-folder
case passes on both. The ordinary-switching check of AI tester 1 (7 steps: sort, filter, row selection, comparison switch) agrees at every step. After the fix, the 3D view was
checked visually in headless Chrome with software WebGL: the orange markers sit on the three residues of the focus window (screenshot `fixed/failing_result_and_back/04-…jpg`).
Limitation: one browser, one machine; the screenshots are cropped JPEG copies of the first 1300 pixel rows.

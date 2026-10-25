"""Figure, caption and short report for the KRAS K1 crystal endpoint, generated from the recorded endpoint only.

    PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_k1_endpoint.py

Reads results/KRAS-K1_endpoint.json (and protocol_params.json for the display of the frozen thresholds). Writes
figures/KRAS-K1_endpoint_residues_2026-10-25.png, its caption and provenance record, and results/KRAS-K1_endpoint_report.md.
Presentation only: nothing is recomputed, and every number comes from the record. No MD coordinate is involved.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RECORD = HERE / "results/KRAS-K1_endpoint.json"
FREEZE = HERE / "results/KRAS-K1_freeze.json"
PARAMS = HERE / "protocol_params.json"
FIG_DIR = HERE / "figures"
STEM = "KRAS-K1_endpoint_residues_2026-10-25"
REPORT = HERE / "results/KRAS-K1_endpoint_report.md"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, BLUE_LIGHT, ORANGE, GREY, HATCH = "#2a78d6", "#a9c8ee", "#d9822b", "#8a8984", "#b9b8b3"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def gaps(scored: list) -> list:
    """(first, last) of each run of residue numbers missing between consecutive scored residues."""
    return [(a + 1, b - 1) for a, b in zip(scored, scored[1:]) if b - a > 1]


def derived(rec: dict, params: dict) -> dict:
    pr = rec["per_residue"]
    md = params["endpoint"]["min_displacement_A"]
    rows = list(zip(pr["resid"], pr["b_mean"], pr["w_max"], pr["threshold"], pr["positive"]))
    segs = rec["gate"]["segments"]
    return {"rows": rows, "high_not_positive": [r for r, b, w, t, p in rows if b >= md and not p],
            "gaps": gaps(rec["scored_residues"]), "segments": segs,
            "qualifying": [s for s in segs if s["length"] >= rec["gate"]["min_segment_length"]]}


def draw(rec: dict, params: dict, d: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    res = [r for r, *_ in d["rows"]]
    b = [x[1] for x in d["rows"]]
    w = [x[2] for x in d["rows"]]
    thr = [x[3] for x in d["rows"]]
    pos = [x[4] for x in d["rows"]]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10.0, 4.9), dpi=200, sharex=True, gridspec_kw={"height_ratios": [4.2, 1.1], "hspace": 0.12})
    for a, z in d["gaps"]:
        for axx in (ax, ax2):
            axx.axvspan(a - 0.5, z + 0.5, facecolor="none", edgecolor=HATCH, hatch="////", linewidth=0, zorder=0)
    ax.bar(res, b, width=0.85, color=[BLUE if p else BLUE_LIGHT for p in pos], zorder=2)
    ax.step(res, w, where="mid", color=GREY, linewidth=1.3, zorder=3)
    ax.step(res, thr, where="mid", color=ORANGE, linewidth=1.3, linestyle="--", zorder=3)
    hi = [(r, bb) for r, bb, p in zip(res, b, pos) if r in d["high_not_positive"]]
    ax.scatter([r for r, _ in hi], [bb for _, bb in hi], s=22, facecolors="white", edgecolors=INK, linewidths=0.9, zorder=4)
    ax.set_ylabel("Cα displacement (Å)", fontsize=8.5)
    top = 1.55 * max(b)
    ax.set_ylim(0, top)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6, zorder=1)
    for a, z in d["gaps"]:
        ax.text((a + z) / 2, 0.9 * top, f"{a}–{z}\nunobserved", ha="center", va="center", fontsize=7, color=MUTED,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    ax.set_title(f"KRAS K1 crystal endpoint, 5US4 (GDP) vs 5USJ (GMPPNP): {rec['n_E']} positive residues in {len(d['segments'])} segments; "
                 f"adequacy gate {'passed' if rec['gate']['passes'] else 'failed'}", fontsize=9.3, loc="left", color=INK)
    handles = [Patch(facecolor=BLUE, label="between-state mean, endpoint-positive"), Patch(facecolor=BLUE_LIGHT, label="between-state mean, not positive"),
               plt.Line2D([0], [0], color=GREY, lw=1.3, label="within-state maximum (copy-to-copy)"),
               plt.Line2D([0], [0], color=ORANGE, lw=1.3, ls="--", label=f"threshold max({params['endpoint']['min_displacement_A']:g} Å, {params['endpoint']['within_state_factor']:g} × within)"),
               plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=INK, label=f"≥ {params['endpoint']['min_displacement_A']:g} Å but not positive"),
               Patch(facecolor="none", edgecolor=HATCH, hatch="////", label="not scored")]
    ax.legend(handles=handles, loc="upper right", ncol=2, fontsize=6.8, frameon=True, framealpha=0.92, edgecolor=GRID)
    min_len = rec["gate"]["min_segment_length"]
    for s in d["segments"]:
        rs = s["residues"]
        ok = s["length"] >= min_len
        for r in rs:
            ax2.broken_barh([(r - 0.5, 1)], (0.2, 0.6), facecolors=BLUE if ok else GREY, edgecolors="none", zorder=2)
        label = f"{rs[0]}–{rs[-1]}: {s['length']} positives" if len(rs) > 1 else f"{rs[0]}: 1"
        ax2.text(rs[-1] + 2.5 if rs[0] < 100 else rs[0] - 2.5, 0.5, label, va="center", ha="left" if rs[0] < 100 else "right", fontsize=7.2, color=INK)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([0.5])
    ax2.set_yticklabels(["segments"], fontsize=8)
    ax2.text(res[-1] + 1, 0.5, f"blue: qualifying (≥ {min_len} positives); grey: not", ha="right", va="center", fontsize=7, color=MUTED)
    ax2.set_xlim(res[0] - 2, res[-1] + 2)
    ax2.set_xlabel("residue number (UniProt P01116-2, KRAS4B numbering; MD topology range 4–164)", fontsize=8.5)
    ax2.tick_params(labelsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
        ax2.spines[sp].set_visible(False)
    fig.text(0.01, 0.005, "External static structural proxy, already exposed during KRAS-STRUCT-2; not held out. Measured from the two crystal "
             "structures only; no MD coordinate was read.", fontsize=6.6, color=MUTED, ha="left", va="bottom")
    fig.subplots_adjust(left=0.1, right=0.985, top=0.93, bottom=0.17)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def caption(rec: dict, params: dict, d: dict) -> str:
    g, e = rec["gate"], params["endpoint"]
    top = 1.55 * max(d["rows"][i][1] for i in range(len(d["rows"])))
    clipped = [r for r, b, w, t, p in d["rows"] if t > top]
    segtxt = "; ".join(f"{s['residues'][0]}–{s['residues'][-1]} ({s['length']})" if len(s["residues"]) > 1 else f"{s['residues'][0]} (1)"
                       for s in d["segments"])
    return f"""# Figure: KRAS K1 crystal endpoint (2026-10-25)

![KRAS K1 endpoint](../figures/{STEM}.png)

**Purpose.** To show, residue by residue, how the frozen crystal-only endpoint of KRAS K1 was formed and why its adequacy gate
{'passed' if g['passes'] else 'failed'}.

**What it shows.** For the {rec['n_scored']} scored residues, the mean Cα displacement between the 5US4 (GDP) and 5USJ (GMPPNP)
chains (bars; {rec['n_between_pairs']} chain pairs), the largest displacement between copies of the same entry (grey line;
{rec['n_within_pairs']} pairs), and the frozen threshold max({e['min_displacement_A']:g} Å, {e['within_state_factor']:g} × within-state) (dashed). Residues whose
bar reaches the threshold are the {rec['n_E']} endpoint-positive residues (dark bars). Merged segments (gap of up to
{e['segments']['merge_gap']} residues): {segtxt}. The gate needs at least {g['min_segments']} segments of at least {g['min_segment_length']} residues; {g['segments_meeting_min_length']} qualify.
Residues {', '.join(f'{a}–{z}' for a, z in d['gaps'])} are unobserved in at least one chain and are not scored. Open circles mark residues with a
mean displacement of at least {e['min_displacement_A']:g} Å that are not positive: for them the within-state value is more than
{1 / e['within_state_factor']:g} of the between-state mean, so the threshold is not reached. The dashed threshold runs above the plotted range at
{('residue ' if len(clipped) == 1 else 'residues ') + ', '.join(str(r) for r in clipped) if clipped else 'no residue'}.

**Main limitation.** The endpoint is an external static structural proxy that the developers already analysed during
KRAS-STRUCT-2; it is not held out and not functional ground truth. It reflects two crystal forms (different space group, pH
and analogue) and omits the unobserved region. The gate tests only whether the proxy contains a sufficiently distributed
comparison target; it does not show that the endpoint can discriminate methods.

**Evidence class.** Measured, from two crystal structures. No MD coordinate was read.

**Source data.** `benchmarks/audit/kras_k1/results/KRAS-K1_endpoint.json` (SHA-256 `{sha256_file(RECORD)}`);
protocol parameters in `protocol_params.json`.

**Reproduce.** `PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_k1_endpoint.py` (about 1.5 CPU-seconds,
matplotlib only). External material: none.
"""


def report(rec: dict, params: dict, d: dict, fz: dict, ver: dict | None = None) -> str:
    g, c, e = rec["gate"], rec["context"], params["endpoint"]
    sens = rec["sensitivity"]
    miss = rec["residues_not_observed_in_all_chains_1_to_169"]

    def rng(v):
        out, cur = [], []
        for r in v:
            if cur and r == cur[-1] + 1:
                cur.append(r)
            else:
                if cur:
                    out.append(cur)
                cur = [r]
        if cur:
            out.append(cur)
        return ", ".join(f"{s[0]}–{s[-1]}" if len(s) > 1 else str(s[0]) for s in out)
    segs = "; ".join(f"{rng(s['residues'])} ({s['length']} positive{'s' if s['length'] != 1 else ''})" for s in g["segments"])
    sens_lines = "\n".join(f"- {k.replace('_', ' ')}: |E| = {v['n']} ({rng(v['E'])}); Jaccard with E {v['jaccard_with_E']:.2f}; "
                           f"segments of at least {g['min_segment_length']} residues: {sum(len(s) >= g['min_segment_length'] for s in v['segments'])}"
                           for k, v in sens.items())
    verdict = rec["decision"]
    return f"""# KRAS K1 crystal endpoint report (step S1, 2026-10-25)

**Decision: {verdict}** — {'application closed' if not g['passes'] else 'owner approval is required before any MD download'}.

Generated from [`KRAS-K1_endpoint.json`](KRAS-K1_endpoint.json) by `plot_k1_endpoint.py`; nothing is typed by hand. Figure:
[`../figures/{STEM}.png`](../figures/{STEM}.png). No MD trajectory was downloaded, opened or scored.

## Result
- **Endpoint E:** {rec['n_E']} residues of {rec['n_scored']} scored: {rng(rec['E'])}.
- **Merged segments** (gap of up to {e['segments']['merge_gap']} residues): {segs}.
- **Gate (frozen, inclusive):** |E| between {g['min_positive']} and {g['max_positive']}: {'met' if g['checks']['n_positive_at_least_min'] and g['checks']['n_positive_at_most_max'] else 'not met'};
  at least {g['min_segments']} segments of at least {g['min_segment_length']} residues: **{'met' if g['checks']['segments_of_min_length_at_least_min'] else 'not met'}** ({g['segments_meeting_min_length']} found).
- **Boundary residues** (within 1e-9 Å of the threshold): {rng(rec['boundary_residues']) if rec['boundary_residues'] else 'none'}.

## Why the endpoint is {'adequate' if g['passes'] else 'not adequate'}
{'The gate passed: the proxy contains a sufficiently distributed comparison target.' if g['passes'] else f"Only {g['segments_meeting_min_length']} of the {len(g['segments'])} segments has at least {g['min_segment_length']} residues, and the frozen gate needs {g['min_segments']}. The proxy therefore does not contain a sufficiently distributed comparison target for a method comparison. K1 is closed with this reason; no alternative chain, threshold, endpoint or system is tried (protocol section 4.2)."}
The gate tests only the spread of the target; it never showed that the endpoint could discriminate methods.

## Displacement context (Å)
- Between-state mean: median {c['median_b']:.2f}, maximum {c['max_b']:.2f}. Within-state (copy-to-copy) maximum: median {c['median_w']:.2f}, maximum {c['max_w']:.2f}.
- {c['n_b_at_least_min_displacement']} residues have a between-state mean of at least {e['min_displacement_A']:g} Å; {c['n_b_at_least_factor_times_w']} have at least {e['within_state_factor']:g} × the within-state value.
- The within-state value exceeds the between-state mean at {c['w_exceeds_b_count']} of {rec['n_scored']} residues; copy-to-copy variation is therefore a limit on what the two entries can say about the nucleotide, although both values are small at most residues (medians above).
- Residues with a mean of at least {e['min_displacement_A']:g} Å that are not positive: {rng(d['high_not_positive'])}.
- Positives whose threshold is set by the within-state value rather than by {e['min_displacement_A']:g} Å: {rng(c['positives_with_threshold_set_by_within_state'])}.

## Exclusions and missing residues
- Scored: {rec['n_scored']} residues (4–59 and 68–164); the observed set within 1–166 equals the frozen KRAS-STRUCT-2 set.
- Residues not observed (of 1–169) per chain: {'; '.join(f'{k} {rng(v)}' for k, v in miss.items())}.
- Not scored because unobserved in at least one chain: {', '.join(f'{a}–{z}' for a, z in d['gaps'])}.

## Sensitivity to the predefined chain handling (reported only; never changes E or the gate)
{sens_lines}

## Not applicable
{'No MD download is proposed: the gate failed and K1 is closed.' if not g['passes'] else 'See the proposal below.'}

## Independent check
{('A separate implementation (' + ver['method'] + ') reproduces the recorded values: maximum absolute difference ' + format(max(ver['max_abs_difference_A'].values()), '.1e') + ' Å over ' + str(ver['n_pairs'][0]) + ' between-state and ' + str(ver['n_pairs'][1]) + ' within-state pairs; endpoint set ' + ('identical' if ver['E_equal'] else 'DIFFERENT') + ' (`KRAS-K1_endpoint_verification.json`).') if ver else 'Not run.'}

## Provenance
- **Frozen commit:** `{rec['frozen_commit']}`; HEAD at run `{rec['head_at_run']}`.
- **SHA-256:** protocol `{fz['protocol']['sha256']}`; parameters `{fz['params']['sha256']}`; endpoint code `{fz['endpoint_code']['sha256']}`;
  tests `{', '.join(f"{k.split('/')[-1]} {v}" for k, v in fz['tests_sha256'].items())}`; freeze record `{rec['freeze_record_sha256']}`.
- **Inputs (equal to the KRAS-STRUCT-2 spec):** {'; '.join(f"{k} `{v['sha256']}`" for k, v in rec['inputs'].items())}.
- **Endpoint record:** `{sha256_file(RECORD)}`.
- **Resources:** {rec['resources']['cpu_s']:.2f} CPU-s, {rec['resources']['wall_s']:.2f} s wall, {rec['resources']['max_rss_kb'] / 1024:.0f} MB peak memory; numpy {rec['software']['numpy']}, Python {rec['software']['python']}.
- **Evidence class:** {rec['evidence_class']}.
"""


def main() -> int:
    rec, params, fz = (json.loads(p.read_text()) for p in (RECORD, PARAMS, FREEZE))
    d = derived(rec, params)
    vp = HERE / "results/KRAS-K1_endpoint_verification.json"
    ver = json.loads(vp.read_text()) if vp.exists() else None
    png = FIG_DIR / f"{STEM}.png"
    draw(rec, params, d, png)
    cap = FIG_DIR / f"{STEM}_caption.md"
    FIG_DIR.mkdir(exist_ok=True)
    cap.write_text(caption(rec, params, d).replace(f"../figures/{STEM}.png", f"{STEM}.png"))
    REPORT.write_text(report(rec, params, d, fz, ver))
    import matplotlib
    prov = {"id": "KRAS-K1-endpoint-figure", "date": "2026-10-25", "kind": "measured (two crystal structures); no MD coordinate",
            "inputs_sha256": {"results/KRAS-K1_endpoint.json": sha256_file(RECORD), "protocol_params.json": sha256_file(PARAMS)},
            "generator": {"path": "plot_k1_endpoint.py", "sha256": sha256_file(Path(__file__).resolve())},
            "output": {"path": f"figures/{STEM}.png", "bytes": png.stat().st_size, "sha256": sha256_file(png)},
            "caption": f"figures/{STEM}_caption.md", "report": "results/KRAS-K1_endpoint_report.md",
            "software": {"python": sys.version.split()[0], "matplotlib": matplotlib.__version__},
            "reproduce": "PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_k1_endpoint.py"}
    (FIG_DIR / f"{STEM}_provenance.json").write_text(json.dumps(prov, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({"png_bytes": prov["output"]["bytes"], "decision": rec["decision"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())

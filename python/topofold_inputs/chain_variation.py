"""Within-entry crystal-copy variation for ``structpair``: several explicitly selected chains per entry.

Pure functions (numpy only); the rules are those of ``docs/STRUCTPAIR_CHAINS_SPEC_2026-10-25.md``. For every mapped position
and channel the table gives the between-entry difference over the chain pairings A x B and the range of the same difference
among the chains inside each entry. The numbers are descriptive.

Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The
within-entry variation is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect
and not an HDX replicate estimate. The advisory flag exists only for a channel whose minimum difference the user configured;
it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.

Output semantics (schema /3, structpair 0.3.0; amendment of 2026-10-25 in docs/STRUCTPAIR_CHAINS_SPEC_2026-10-25.md): the numeric
between-entry and within-entry values are always reported; the flag is left unassessed unless a finite minimum difference above
the numerical tolerance is configured for that channel; equality with the threshold meets the rule, with no tolerance widening;
values within the tolerance of the threshold are marked as at the boundary. Records of schema /2 (0.2.0) keep their historical
meaning (default floor 0, tolerance-widened comparison, column `separated`).
"""

from __future__ import annotations

import itertools
import math
import statistics

import numpy as np

TOL = 1e-9                                  # numerical tolerance: the smallest admissible minimum difference, and the boundary marker
SEMANTICS = "advisory/2"                    # schema /3; "advisory/1" is the meaning of schema /2 records
CHANNELS = ("disp", "d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4")
UNIT = {"disp": "angstrom", "d_kappa": "rad", "d_i_i2": "angstrom", "rmsd3": "angstrom", "d_tau": "rad",
        "d_i_i3": "angstrom", "rmsd4": "angstrom"}
FOUR = ("d_tau", "d_i_i3", "rmsd4")
STATISTICS = ("min", "median", "mean", "max")
DEFAULT_TRIM = 0.2
DEFAULT_ADVISORY = {"statistic": "min", "within_factor": 1.0, "min_difference": {}}
MIN_FIT = 4
STATEMENT = ("Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of "
             "them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the "
             "user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand "
             "effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference "
             "the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not "
             "biological separation.")
DEFINITIONS = {
    "disp": "Cα displacement at the position after a trimmed Kabsch fit of the two chains on the positions present in both "
            "(drop the floor(trim × n) most displaced, ties: lower position first; refit once); the K1 endpoint calculation",
    "d_kappa": "|κ difference|, three-residue window (existing structpair metric, computed for the chain pair)",
    "d_i_i2": "|d(i−1,i+1) difference|, three-residue window (existing metric)",
    "rmsd3": "local Cα RMSD, three-residue window, MDAnalysis QCP (existing metric)",
    "d_tau": "circular |τ difference|, four-residue window (existing metric)",
    "d_i_i3": "|d(i−1,i+2) difference|, four-residue window (existing metric)",
    "rmsd4": "local Cα RMSD, four-residue window, MDAnalysis QCP (existing metric)"}


class VariationError(ValueError):
    pass


# --- geometry (one Kabsch for structpair and for the displacement) ----------------------------------------------------

def kabsch_frame(moving: np.ndarray, reference: np.ndarray) -> tuple:
    """(R, centre of moving, centre of reference): proper rotation (det +1) with (moving - cm) @ R.T + cr closest to reference."""
    cm, cr = moving.mean(0), reference.mean(0)
    u, _, vt = np.linalg.svd((moving - cm).T @ (reference - cr))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    if d == 0:
        d = 1.0
    return vt.T @ np.diag([1.0, 1.0, d]) @ u.T, cm, cr


def trim_indices(disp: np.ndarray, m: int) -> list:
    """Indices of the m most displaced positions; equal displacements drop the lower position first."""
    return sorted(range(len(disp)), key=lambda i: (-float(disp[i]), i))[:m]


def pair_displacement(a: np.ndarray, b: np.ndarray, trim_fraction: float = DEFAULT_TRIM) -> np.ndarray:
    """Per-position Cα displacement (Å) of chain b from chain a after the trimmed fit; NaN where either Cα is missing.

    The fit uses the positions with a finite Cα in both chains (at least MIN_FIT, otherwise every value is NaN)."""
    out = np.full(len(a), np.nan)
    ok = np.isfinite(a).all(1) & np.isfinite(b).all(1)
    idx = np.where(ok)[0]
    if len(idx) < MIN_FIT:
        return out
    ref, mov = a[idx], b[idx]
    r, cm, cr = kabsch_frame(mov, ref)
    first = np.linalg.norm(ref - ((mov - cm) @ r.T + cr), axis=1)
    drop = set(trim_indices(first, int(math.floor(trim_fraction * len(idx) + 1e-9))))
    keep = [i for i in range(len(idx)) if i not in drop]
    if len(keep) < MIN_FIT:
        return out
    r2, cm2, cr2 = kabsch_frame(mov[keep], ref[keep])
    out[idx] = np.linalg.norm(ref - ((mov - cm2) @ r2.T + cr2), axis=1)
    return out


# --- parameters ---------------------------------------------------------------------------------------------------------

def _finite_number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def parse_params(block: dict | None) -> dict:
    """Validated parameters. Every number must be finite (NaN and infinity are refused); a minimum difference must exceed TOL."""
    block = dict(block or {})
    unknown = set(block) - {"displacement_trim_fraction", "advisory"}
    if unknown:
        raise VariationError(f"chain_variation: unknown keys {sorted(unknown)}")
    trim = block.get("displacement_trim_fraction", DEFAULT_TRIM)
    if not _finite_number(trim) or not 0.0 <= trim < 0.5:
        raise VariationError("chain_variation.displacement_trim_fraction must be a finite number in [0, 0.5)")
    adv = dict(DEFAULT_ADVISORY, **dict(block.get("advisory") or {}))
    extra = set(adv) - set(DEFAULT_ADVISORY)
    if extra:
        raise VariationError(f"chain_variation.advisory: unknown keys {sorted(extra)}")
    if adv["statistic"] not in STATISTICS:
        raise VariationError(f"chain_variation.advisory.statistic must be one of {STATISTICS}")
    f = adv["within_factor"]
    if not _finite_number(f) or f < 0:
        raise VariationError("chain_variation.advisory.within_factor must be a finite number >= 0")
    md = adv["min_difference"]
    if not isinstance(md, dict) or any(k not in CHANNELS or not _finite_number(v) or v <= TOL for k, v in md.items()):
        raise VariationError(f"chain_variation.advisory.min_difference maps channels {CHANNELS} to finite numbers above the "
                             f"numerical tolerance {TOL:g} (a flag is assessed only for a channel with a positive minimum "
                             "difference, in that channel's unit)")
    return {"displacement_trim_fraction": float(trim),
            "advisory": {"statistic": adv["statistic"], "within_factor": float(f),
                         "min_difference": {k: float(v) for k, v in md.items()}}}


# --- aggregation ---------------------------------------------------------------------------------------------------------

def aggregate(values) -> dict:
    """n, min, median, mean, max over the defined (not None, finite) values; empty statistics when n = 0."""
    v = [float(x) for x in values if x is not None and math.isfinite(x)]
    if not v:
        return {"n": 0, "min": None, "median": None, "mean": None, "max": None}
    return {"n": len(v), "min": min(v), "median": statistics.median(v), "mean": math.fsum(v) / len(v), "max": max(v)}


def pair_keys(chains_a: list, chains_b: list) -> dict:
    """Cross pairings (A chain, B chain) and the unordered within-entry pairs, in the order of the chain lists."""
    return {"cross": [("cross", a, b) for a in chains_a for b in chains_b],
            "withinA": [("withinA", x, y) for x, y in itertools.combinations(chains_a, 2)],
            "withinB": [("withinB", x, y) for x, y in itertools.combinations(chains_b, 2)]}


def advisory(stat_value, within_max, between_n: int, within_scope: str, channel: str, params: dict) -> dict:
    """The advisory rule of the amendment (semantics advisory/2).

    Assessed only for a channel with a configured minimum difference `m` (finite, above TOL). Then, with `f` = within_factor and
    `S` the chosen between-entry statistic: flag = S >= max(m, f x within_max), an exact float64 comparison (equality meets the
    rule; no tolerance widening). `at_boundary` is true when |S - threshold| <= TOL, so a flag that hinges on rounding is
    visible. Never False when it cannot be assessed: it is left empty with a reason."""
    floor = params["min_difference"].get(channel)
    if floor is None:
        return {"threshold": None, "flag": None, "at_boundary": None,
                "reason": "no minimum difference configured for this channel (the numeric values are reported)"}
    if between_n == 0:
        return {"threshold": None, "flag": None, "at_boundary": None, "reason": "no defined cross pairing"}
    if within_scope == "none" or within_max is None:
        return {"threshold": None, "flag": None, "at_boundary": None, "reason": "no defined within-entry pair"}
    thr = max(floor, params["within_factor"] * within_max)
    return {"threshold": thr, "flag": bool(stat_value >= thr), "at_boundary": bool(abs(stat_value - thr) <= TOL), "reason": ""}


def variation_rows(pair_values: dict, n: int, labels: list, windows: list, chains_a: list, chains_b: list,
                   primary: tuple, params: dict) -> list:
    """Long table: one row per position and channel.

    pair_values[(kind, c1, c2)][channel][i] is the value of one chain pair (None where undefined); ``labels[i]`` is
    (A residue, B residue); ``windows[i]`` is (three-residue window label, four-residue window label)."""
    keys = pair_keys(chains_a, chains_b)
    adv = params["advisory"]
    rows = []
    for i in range(n):
        for ch in CHANNELS:
            def vals(group):
                return [pair_values[k][ch][i] for k in keys[group]]
            x, wa, wb = aggregate(vals("cross")), aggregate(vals("withinA")), aggregate(vals("withinB"))
            maxima = [w["max"] for w in (wa, wb) if w["n"]]
            within_max = max(maxima) if maxima else None
            scope = {(True, True): "both", (True, False): "A_only", (False, True): "B_only", (False, False): "none"}[(bool(wa["n"]), bool(wb["n"]))]
            s_val = x[adv["statistic"]]
            a = advisory(s_val, within_max, x["n"], scope, ch, adv)
            win = labels[i][0] if ch == "disp" else windows[i][1 if ch in FOUR else 0]
            rows.append({"position": i, "A_residue": labels[i][0], "B_residue": labels[i][1], "channel": ch, "unit": UNIT[ch],
                         "window_A": win, "primary": pair_values[("cross",) + primary][ch][i],
                         "between_n": x["n"], "between_min": x["min"], "between_median": x["median"], "between_mean": x["mean"],
                         "between_max": x["max"], "withinA_n": wa["n"], "withinA_min": wa["min"], "withinA_max": wa["max"],
                         "withinB_n": wb["n"], "withinB_min": wb["min"], "withinB_max": wb["max"], "within_max": within_max,
                         "within_scope": scope, "advisory_statistic": s_val, "advisory_threshold": a["threshold"],
                         "advisory_flag": a["flag"], "advisory_at_boundary": a["at_boundary"],
                         "advisory_not_assessed_reason": a["reason"]})
    return rows


def summarize(rows: list, chains_a: list, chains_b: list, primary: tuple, params: dict) -> dict:
    keys = pair_keys(chains_a, chains_b)
    status = ("available" if len(chains_a) >= 2 and len(chains_b) >= 2 else "partial" if len(chains_a) >= 2 or len(chains_b) >= 2
              else "unavailable")
    per = {}
    for ch in CHANNELS:
        r = [x for x in rows if x["channel"] == ch]
        per[ch] = {"positions_with_a_value": sum(x["between_n"] > 0 for x in r),
                   "advisory_configured": ch in params["advisory"]["min_difference"],
                   "advisory_flag_true": sum(x["advisory_flag"] is True for x in r),
                   "advisory_flag_false": sum(x["advisory_flag"] is False for x in r),
                   "advisory_not_assessed": sum(x["advisory_flag"] is None and x["between_n"] > 0 for x in r),
                   "advisory_at_boundary": sum(x["advisory_at_boundary"] is True for x in r)}
    return {"status": status, "chains": {"A": chains_a, "B": chains_b}, "primary": {"A": primary[0], "B": primary[1]},
            "n_pairs": {"cross": len(keys["cross"]), "withinA": len(keys["withinA"]), "withinB": len(keys["withinB"])},
            "parameters": params, "advisory_semantics": SEMANTICS, "channels": per, "definitions": DEFINITIONS,
            "statement": STATEMENT}


def unavailable_block() -> dict:
    return {"status": "unavailable", "reason": "one chain per entry", "statement": STATEMENT}


def report_lines(block: dict, rows: list, top_n: int) -> list:
    """The report section; prominent, and the statement comes first."""
    L = ["## Within-entry crystal-copy variation", "", f"**{block['statement']}**", ""]
    if block["status"] == "unavailable":
        return L + ["**Unavailable: one chain per entry.** Name every chain of an entry in `structures.X.chains` (with "
                    "`primary_chain`) to obtain it.", ""]
    adv = block["parameters"]["advisory"]
    n = block["n_pairs"]
    floors = adv["min_difference"]
    L += [f"- Status: **{block['status']}**. Chains A {block['chains']['A']} (primary {block['primary']['A']}), "
          f"B {block['chains']['B']} (primary {block['primary']['B']}). Pairs: {n['cross']} cross, {n['withinA']} within A, "
          f"{n['withinB']} within B. The legacy tables above are the primary pair.",
          "- The numeric between-entry and within-entry values are always given in `chain_variation.csv`; an empty cell is "
          "unavailable, never zero, and the `*_n` columns count the chain pairs that contributed.",
          ("- Advisory rule (not significance, not biological separation), configured by the user for "
           + ", ".join(f"{c} (minimum difference {v:g} {UNIT[c].replace('angstrom', 'Å')})" for c, v in floors.items())
           + f": the between-entry **{adv['statistic']}** is at least max(minimum difference, {adv['within_factor']:g} × the "
           "within-entry maximum); equality meets the rule and a value within 1e-9 of the threshold is marked as at the boundary."
           if floors else
           "- **No advisory rule is configured, so no flag is assessed.** Set `chain_variation.advisory.min_difference` to a finite "
           "positive number for a channel (in its unit) to obtain a flag for that channel."),
          f"- Displacement trim fraction {block['parameters']['displacement_trim_fraction']:g}.", "",
          "| Channel | Unit | Positions with a value | Advisory rule | Rule met | Rule not met | Not assessed |", "|---|---|---|---|---|---|---|"]
    for ch in CHANNELS:
        c = block["channels"][ch]
        L.append(f"| {ch} | {UNIT[ch].replace('angstrom', 'Å')} | {c['positions_with_a_value']} | "
                 f"{'configured' if c['advisory_configured'] else 'not configured'} | {c['advisory_flag_true']} | "
                 f"{c['advisory_flag_false']} | {c['advisory_not_assessed']} |")
    d = sorted((r for r in rows if r["channel"] == "disp" and r["between_median"] is not None), key=lambda r: (-r["between_median"], r["position"]))[:top_n]
    L += ["", f"Largest between-entry displacement (median over cross pairings) with the within-entry maximum (top {top_n}):", "",
          "| A residue | Between median [min–max] (Å) | Within-entry max (Å) | Scope | Advisory rule |", "|---|---|---|---|---|"]
    for r in d:
        flag = "not assessed" if r["advisory_flag"] is None else ("met" if r["advisory_flag"] else "not met")
        w = "" if r["within_max"] is None else f"{r['within_max']:.3f}"
        L.append(f"| {r['A_residue']} | {r['between_median']:.3f} [{r['between_min']:.3f}–{r['between_max']:.3f}] | {w} | "
                 f"{r['within_scope']} | {flag} |")
    L += ["", "Files: `chain_variation.csv` (long format; every statistic with its count), `chains.csv`, `chain_variation.png`.", ""]
    return L


def figure(rows: list, path) -> str:
    """One compact panel: displacement between entries (min–max band, median) against the within-entry maxima."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return "unavailable: matplotlib is not installed in this interpreter"
    r = [x for x in rows if x["channel"] == "disp"]
    xs = np.arange(len(r))

    def col(k):
        return np.array([np.nan if x[k] is None else x[k] for x in r], dtype=float)
    fig, ax = plt.subplots(figsize=(9.5, 3.9), dpi=150)
    ax.fill_between(xs, col("between_min"), col("between_max"), color="#c7dcf5", lw=0, label="between entries: min–max over cross pairings")
    ax.plot(xs, col("between_median"), color="#2a78d6", lw=1.4, label="between entries: median")
    from matplotlib.patches import Patch
    extra = []
    for k, ent, lab, c in (("withinA_max", "A", "within A: maximum over A chain pairs", "#52514e"),
                           ("withinB_max", "B", "within B: maximum over B chain pairs", "#d9822b")):
        if np.isfinite(col(k)).any():
            ax.plot(xs, col(k), color=c, lw=1.1, ls="--", label=lab)
        else:                                                       # no chain pair inside this entry: missing, never zero
            extra.append(Patch(facecolor="none", edgecolor=c, hatch="////", label=f"within {ent}: unavailable (one chain), missing, not zero"))
    top = np.nanmax(np.concatenate([col(k) for k in ("between_max", "withinA_max", "withinB_max")]))
    top = 1.0 if not np.isfinite(top) else float(top)
    met = [i for i, x in enumerate(r) if x["advisory_flag"] is True]
    if met:
        ax.scatter(met, [1.06 * top] * len(met), marker="v", s=18, color="#0b0b0b", label="advisory rule met (configured)", zorder=3,
                   clip_on=False)
    elif not any(x["advisory_flag"] is not None for x in r):
        ax.text(0.01, 0.96, "no advisory rule configured: no flag assessed", transform=ax.transAxes, fontsize=6.8, ha="left",
                va="top", color="#52514e")
    num = []
    for x in r:
        try:
            num.append(int(str(x["A_residue"]).split(":")[-1]))
        except ValueError:
            num.append(None)
    for i in range(1, len(num)):
        if num[i] is not None and num[i - 1] is not None and num[i] - num[i - 1] != 1:
            ax.axvline(i - 0.5, color="#8a8984", lw=0.8, ls=":")
            ax.text(i - 2.0, 0.8 * top, "jump in\nnumbering", fontsize=6.5, color="#52514e", va="center", ha="right")
    n = len(r)
    ticks = list(range(0, n, max(1, n // 8)))
    ax.set_xticks(ticks)
    ax.set_xticklabels([r[t]["A_residue"].split(":")[-1] for t in ticks], fontsize=8)
    ax.set_xlim(0, max(n - 1, 1))
    ax.set_ylim(0, 1.12 * top)
    ax.set_ylabel("Cα displacement (Å)", fontsize=9)
    ax.set_xlabel("A residue number (one tick label per position; dotted line: jump in residue numbers)", fontsize=9)
    h, lab = ax.get_legend_handles_labels()
    ax.legend(h + extra, lab + [e.get_label() for e in extra], fontsize=7, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.93))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("Within-entry crystal-copy variation vs between-entry displacement · descriptive", fontsize=9, loc="left")
    fig.text(0.01, 0.01, "Chains in one crystal are not independent replicates; not a calibrated noise floor or an uncertainty.",
             fontsize=6.8, color="#52514e")
    fig.subplots_adjust(left=0.08, right=0.985, top=0.9, bottom=0.2)
    fig.savefig(path)
    plt.close(fig)
    return "written"

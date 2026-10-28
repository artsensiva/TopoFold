"""MET exploratory evaluation: position-blocked, paired comparison of feature sets (plan: PLAN.md / plan.json; design constants: DESIGN below).

    python3 evaluate.py freeze   # hash-freezes plan, code, tests and input hashes
    python3 evaluate.py run      # one scoring attempt (resumable from results/progress.jsonl; design unchanged)

Scope: internal, exploratory, descriptive. The target population is variants that pass the measured no-drug eligibility rule; nothing here
shows prediction for variants whose baseline function is unknown. Resource accounting covers this whole process (RUSAGE_SELF plus
RUSAGE_CHILDREN; no child process is started); threads are bounded to one before numpy is imported.
"""

from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_k] = "1"

import csv  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import re  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
EXPORT = HERE.parent / "met_feature_export_2026-10-27"

DESIGN = {
    "drugs": ["Crizo", "Savo", "Mere", "A458", "Tepo", "NVP"],
    "n_positions": 287,
    "outer_blocks": 10,
    "buffer": 8,
    "inner_blocks": 5,
    "min_train_rows": 50,
    "min_inner_train_rows": 20,
    "alphas": [0.1, 1.0, 10.0, 100.0, 1000.0],
    "default_alpha_if_no_inner_fold": 10.0,
    "gbm": {"max_depth": 3, "learning_rate": 0.05, "max_iter": 200, "min_samples_leaf": 20, "l2_regularization": 1.0, "random_state": 0, "early_stopping": False},
    "bootstrap": {"block_length": 20, "resamples": 2000, "seed": 20261028, "level": 0.95},
    "criterion": {"min_blocks_improving": 8, "min_drugs_improving": 5, "no_criterion_max_blocks_improving": 5},
    "continuation_relative_mse_reduction": 0.05,
    "cpu_seconds_cap": 1800,
    "rss_mb_cap": 2048,
    "feature_sets": {"F0": ["M"], "F1": ["M", "L"], "F2": ["M", "L", "C"], "F3": ["M", "L", "T"], "F4": ["M", "L", "ESM"], "F5": ["M", "L", "ESM", "T"]},
    "comparisons": [("C1", "F1", "F0"), ("C2", "F2", "F1"), ("C3", "F3", "F1"), ("C4", "F3", "F2"), ("C5", "F5", "F4"), ("C6", "F4", "F1")],
    "models": ["ridge", "gbm"],
}
M_COLS = ["vol_wt", "vol_mut", "dvol", "kd_wt", "kd_mut", "dkd", "charge_wt", "charge_mut", "dcharge", "blosum62", "to_pro", "from_pro", "to_gly", "from_gly"]
L_COLS = ["L_min_heavy_dist_A", "L_min_sidechain_dist_A", "L_n_contact_atoms_4p5", "L_ca_min_dist_A", "L_no_side_chain"]
T_COLS = [f"kappa_v{i}" for i in range(1, 8)] + [f"tau_sin_h{i}" for i in range(1, 7)] + [f"tau_cos_h{i}" for i in range(1, 7)]
C_COLS = [f"d_{i}_{j}" for i in range(9) for j in range(i + 1, 9)]
BLOCK_COLS = {"M": M_COLS, "L": L_COLS, "T": T_COLS, "C": C_COLS, "ESM": ["esm_llr"]}


def sha256_file(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# --- inputs -------------------------------------------------------------------------------------------------------------------

def read_assay(path, drugs):
    """{(inhibitor, position, mutation): (effect, lof_flag_is_true)} for missense rows with a numeric effect; duplicate keys raise."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["type"] != "missense" or r["inhibitor"] not in drugs + ["DMSO"]:
            continue
        try:
            v = float(r["ROSACE_effects"])
        except ValueError:
            continue
        if not math.isfinite(v):
            continue
        key = (r["inhibitor"], int(r["position"]), r["mutation"])
        if key in out:
            raise ValueError(f"duplicate assay row {key}")
        out[key] = (v, r["ROSACE.LOF"] == "TRUE", r["wildtype"])
    return out


def gamma_table(assay, drug):
    """Eligibility and endpoint: gamma = effect_inh - effect_DMSO for variants with both numeric effects and DMSO ROSACE.LOF == FALSE."""
    out = {}
    for (inh, pos, mut), (eff, _, wt) in assay.items():
        if inh != drug:
            continue
        d = assay.get(("DMSO", pos, mut))
        if d is None or d[1]:
            continue
        out[(pos, mut)] = (eff - d[0], wt)
    return out


def read_esm(path):
    """{(dms_position, mutant): score}: only the pos_mut and score columns are used (score is invariant across inhibitors, checked)."""
    out, seen = {}, {}
    for r in csv.DictReader(open(path)):
        m = re.fullmatch(r"(\d+)([A-Z])", r["pos_mut"])
        key = (int(m.group(1)), m.group(2))
        v = float(r["score"])
        if key in seen and seen[key] != v:
            raise ValueError(f"score differs across rows for {key}")
        seen[key] = v
        out[key] = v
    return out


def read_sd(path, drugs):
    """{(condition, position, mutation): sd} from the filtered score table (approximate uncertainty context only)."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["type"] == "missense" and r["key"] in drugs + ["DMSO"]:
            out[(r["key"], int(r["position"]), r["mutation"])] = float(r["sd"])
    return out


def load_export():
    spec = importlib.util.spec_from_file_location("met_export", EXPORT / "met_export.py")
    me = importlib.util.module_from_spec(spec)
    sys.modules["met_export"] = me
    spec.loader.exec_module(me)
    return me


def num(x):
    return float("nan") if x == "" else float(x)


def structure_lookup(path):
    out = {}
    for r in csv.DictReader(open(path)):
        out[(r["drug"], int(r["dms_position"]))] = r
    return out


def mutation_lookup(path):
    return {(r["wt"], r["mut"]): r for r in csv.DictReader(open(path))}


def build_drug_arrays(drug, gamma, struct, mutation, esm, sd_cond=None, n_positions=287):
    """Rows of one drug: eligible variants with every structural block and an ESM score. Returns dict of arrays plus counts."""
    keys = sorted(gamma)
    rows, counts = [], {"eligible": len(keys), "no_structure": 0, "no_esm": 0, "unknown_mutation": 0}
    for (pos, mut) in keys:
        y, wt = gamma[(pos, mut)]
        s = struct.get((drug, pos))
        m = mutation.get((wt, mut))
        if m is None or s is None:
            counts["unknown_mutation"] += 1
            continue
        if s["T_reason"] != "ok" or s["L_reason"] not in ("ok", "no_side_chain"):
            counts["no_structure"] += 1
            continue
        if (pos, mut) not in esm:
            counts["no_esm"] += 1
            continue
        gly = s["L_reason"] == "no_side_chain"
        lvals = [num(s["L_min_heavy_dist_A"]), num(s["L_min_heavy_dist_A"]) if gly else num(s["L_min_sidechain_dist_A"]), num(s["L_n_contact_atoms_4p5"]), num(s["L_ca_min_dist_A"]), 1.0 if gly else 0.0]
        feats = {"M": [float(m[c]) for c in M_COLS], "L": lvals, "T": [num(s[c]) for c in T_COLS], "C": [num(s[c]) for c in C_COLS], "ESM": [esm[(pos, mut)]]}
        if not all(math.isfinite(v) for b in feats.values() for v in b):
            raise ValueError(f"non-finite feature in an exported-complete row {drug} {pos} {mut}")
        sd = None
        if sd_cond is not None and (drug, pos, mut) in sd_cond and ("DMSO", pos, mut) in sd_cond:
            sd = math.sqrt(sd_cond[(drug, pos, mut)] ** 2 + sd_cond[("DMSO", pos, mut)] ** 2)
        rows.append((pos, mut, y, feats, sd))
    out = {"pos": np.array([r[0] for r in rows], int), "mut": [r[1] for r in rows], "y": np.array([r[2] for r in rows], float),
           "sd": np.array([np.nan if r[4] is None else r[4] for r in rows], float),
           "X": {b: np.array([r[3][b] for r in rows], float).reshape(len(rows), len(BLOCK_COLS[b])) for b in BLOCK_COLS}, "counts": dict(counts, used=len(rows))}
    out["row_id"] = [f"{drug}:{p}{m}" for p, m in zip(out["pos"], out["mut"])]
    return out


# --- folds ----------------------------------------------------------------------------------------------------------------------

def block_assignment(n_positions=287, n_blocks=10):
    """Contiguous blocks of positions 1..n; sizes differ by at most one (the first n % n_blocks blocks are larger). Returns array fold[pos-1]."""
    base, extra = divmod(n_positions, n_blocks)
    sizes = [base + 1 if i < extra else base for i in range(n_blocks)]
    return np.repeat(np.arange(n_blocks), sizes)


def train_test_masks(pos, fold_of_pos, k, buffer):
    """Test = rows in block k. Train = rows whose position is farther than `buffer` residues from every held-out position."""
    held = np.where(fold_of_pos == k)[0] + 1
    test = fold_of_pos[pos - 1] == k
    near = (pos >= held.min() - buffer) & (pos <= held.max() + buffer)
    return ~near, test


def position_mean_sq(pos, err2):
    """{position: mean squared error over its rows}."""
    d = {}
    for p, e in zip(pos, err2):
        d.setdefault(int(p), []).append(float(e))
    return {p: float(np.mean(v)) for p, v in d.items()}


def position_weighted_mse(pos, err2):
    """Mean over positions of the per-position mean squared error (equal weight per position)."""
    pm = position_mean_sq(pos, err2)
    return float(np.mean(list(pm.values()))) if pm else float("nan")


# --- models ---------------------------------------------------------------------------------------------------------------------

def standardize(Xtr, Xte):
    mu = Xtr.mean(0)
    sd = Xtr.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    return (Xtr - mu) / sd, (Xte - mu) / sd


def ridge_predict(Xtr, ytr, Xte, alpha):
    from sklearn.linear_model import Ridge
    A, B = standardize(Xtr, Xte)
    return Ridge(alpha=alpha).fit(A, ytr).predict(B)


def gbm_predict(Xtr, ytr, Xte):
    from sklearn.ensemble import HistGradientBoostingRegressor
    return HistGradientBoostingRegressor(**DESIGN["gbm"]).fit(Xtr, ytr).predict(Xte)


def inner_blocks(train_positions, n_inner):
    """Contiguous groups of the sorted training positions."""
    u = np.array(sorted(set(int(p) for p in train_positions)))
    return [g for g in np.array_split(u, n_inner) if len(g)]


def select_alpha(X, y, pos, train, alphas=None, buffer=None, n_inner=None):
    """Alpha with the smallest position-weighted validation loss over inner position blocks (ties -> larger alpha). Inner folds with an empty
    validation set or fewer than min_inner_train_rows are skipped; no usable inner fold -> default alpha. Returns (alpha, losses|None, n_used)."""
    alphas = alphas or DESIGN["alphas"]
    buffer = DESIGN["buffer"] if buffer is None else buffer
    n_inner = n_inner or DESIGN["inner_blocks"]
    tr_idx = np.where(train)[0]
    losses = {a: [] for a in alphas}
    used = 0
    for g in inner_blocks(pos[tr_idx], n_inner):
        val = np.isin(pos, g) & train
        far = ((pos < g.min() - buffer) | (pos > g.max() + buffer)) & train
        if val.sum() == 0 or far.sum() < DESIGN["min_inner_train_rows"]:
            continue
        used += 1
        for a in alphas:
            p = ridge_predict(X[far], y[far], X[val], a)
            losses[a].append(position_weighted_mse(pos[val], (p - y[val]) ** 2))
    if used == 0:
        return DESIGN["default_alpha_if_no_inner_fold"], None, 0
    mean = {a: float(np.mean(v)) for a, v in losses.items()}
    best = min(mean.values())
    return max(a for a in alphas if mean[a] <= best + 0.0), mean, used


def feature_matrix(arrs, blocks):
    return np.hstack([arrs["X"][b] for b in blocks])


def run_unit(arrs, fold_of_pos, k):
    """One (drug, outer fold): predictions for the held-out rows for every feature set and model. Cells that cannot be fitted use the
    training-mean fallback for every feature set alike and are flagged; zero training rows or zero test rows are flagged and returned empty."""
    pos, y = arrs["pos"], arrs["y"]
    train, test = train_test_masks(pos, fold_of_pos, k, DESIGN["buffer"])
    info = {"train_rows": int(train.sum()), "test_rows": int(test.sum()), "flag": "ok", "alpha": {}, "inner_folds_used": {}}
    if test.sum() == 0:
        info["flag"] = "empty_test"
        return {}, info, test
    if train.sum() == 0:
        info["flag"] = "no_training_rows"
        return {}, info, test
    preds = {}
    if train.sum() < DESIGN["min_train_rows"]:
        info["flag"] = "insufficient_training_fallback_mean"
        for s in DESIGN["feature_sets"]:
            for m in DESIGN["models"]:
                preds[(s, m)] = np.full(int(test.sum()), float(y[train].mean()))
        return preds, info, test
    for s, blocks in DESIGN["feature_sets"].items():
        X = feature_matrix(arrs, blocks)
        a, _, used = select_alpha(X, y, pos, train)
        info["alpha"][s], info["inner_folds_used"][s] = a, used
        preds[(s, "ridge")] = ridge_predict(X[train], y[train], X[test], a)
        preds[(s, "gbm")] = gbm_predict(X[train], y[train], X[test])
    return preds, info, test


# --- aggregation ------------------------------------------------------------------------------------------------------------------

def error_matrix(drug_data, preds, n_positions, mask=None):
    """E[d, p-1] = mean squared error at position p of drug d (NaN where the drug has no row there); P[d, p-1] = presence."""
    E = np.full((len(drug_data), n_positions), np.nan)
    for i, (d, a) in enumerate(drug_data.items()):
        keep = np.ones(len(a["pos"]), bool) if mask is None else mask[d]
        if keep.sum() == 0:
            continue
        pm = position_mean_sq(a["pos"][keep], (preds[d] - a["y"])[keep] ** 2)
        for p, v in pm.items():
            E[i, p - 1] = v
    return E


def overall_mse(E):
    per_drug = np.nanmean(E, axis=1)
    return float(np.nanmean(per_drug)), per_drug


def block_deltas(Ea, Eb, fold_of_pos, n_blocks):
    out = []
    for k in range(n_blocks):
        cols = np.where(fold_of_pos == k)[0]
        a = np.nanmean(np.nanmean(Ea[:, cols], axis=1))
        b = np.nanmean(np.nanmean(Eb[:, cols], axis=1))
        out.append(float(a - b))
    return out


def bootstrap_deltas(Ea, Eb, n_positions, cfg):
    """Moving-block bootstrap over positions (all drugs and models resampled together). Returns the array of Delta over resamples."""
    rng = np.random.default_rng(cfg["seed"])
    L, B = cfg["block_length"], cfg["resamples"]
    nb = math.ceil(n_positions / L)
    W = np.zeros((B, n_positions))
    for b in range(B):
        starts = rng.integers(0, n_positions - L + 1, size=nb)
        idx = np.concatenate([np.arange(s, s + L) for s in starts])[:n_positions]
        np.add.at(W[b], idx, 1.0)

    def overall(E):
        present = (~np.isnan(E)).astype(float)
        num = W @ np.nan_to_num(E).T
        den = W @ present.T
        with np.errstate(invalid="ignore", divide="ignore"):
            per = num / den
        return np.nanmean(per, axis=1)

    return overall(Ea) - overall(Eb)


def classify(delta, hi, blocks_improving, drugs_improving, n_blocks=10, n_drugs=6):
    """'criterion met' / 'criterion not met' / 'inconclusive' (the criterion is statistical consistency; it is not equivalence, and not usefulness)."""
    c = DESIGN["criterion"]
    if delta < 0 and hi < 0 and blocks_improving >= c["min_blocks_improving"] and drugs_improving >= c["min_drugs_improving"]:
        return "criterion met"
    if delta >= 0 or blocks_improving <= c["no_criterion_max_blocks_improving"]:
        return "criterion not met"
    return "inconclusive"


def family_outcome(ridge_class, gbm_class):
    if ridge_class == "criterion met" and gbm_class == "criterion met":
        return "criterion met in both families"
    if "criterion met" in (ridge_class, gbm_class):
        return "criterion met in one family only (family-dependent; not supported)"
    return "criterion not met in either family (not-met or inconclusive)"


def engineering_outcome(res):
    """Provisional continuation rule (fixed before scoring). res[comparison][family] = dict(rel, cls)."""
    thr = DESIGN["continuation_relative_mse_reduction"]
    c3, c4 = res["C3"], res["C4"]
    a = all(c3[f]["cls"] == "criterion met" and c3[f]["rel"] >= thr for f in DESIGN["models"])
    if not a:
        return "do not continue geometric-feature development on this evidence (not equivalence; not proof of no effect)"
    if all(c4[f]["cls"] == "criterion met" for f in DESIGN["models"]):
        return "continue; angular encoding favoured over equal-coordinate distances"
    return "continue with geometry; angular encoding not distinguished from distances"


def evaluate_predictions(drug_data, preds, fold_of_pos, n_positions, mask=None):
    """All paired comparisons from out-of-fold predictions. preds[drug][(set, model)] = array over the drug's rows."""
    drugs = list(drug_data)
    E = {}
    for s in DESIGN["feature_sets"]:
        for m in DESIGN["models"]:
            E[(s, m)] = error_matrix(drug_data, {d: preds[d][(s, m)] for d in drugs}, n_positions, mask)
    per = {f"{s}|{m}": {"overall_mse": overall_mse(E[(s, m)])[0], "per_drug_mse": dict(zip(drugs, overall_mse(E[(s, m)])[1].tolist()))} for (s, m) in E}
    comps = {}
    for cid, a, b in DESIGN["comparisons"]:
        comps[cid] = {"test": f"{a} vs {b}", "families": {}}
        for m in DESIGN["models"]:
            Ea, Eb = E[(a, m)], E[(b, m)]
            ma, pa = overall_mse(Ea)
            mb, pb = overall_mse(Eb)
            delta = ma - mb
            boots = bootstrap_deltas(Ea, Eb, n_positions, DESIGN["bootstrap"])
            lo, hi = np.percentile(boots, [100 * (1 - DESIGN["bootstrap"]["level"]) / 2, 100 * (1 + DESIGN["bootstrap"]["level"]) / 2])
            bd = block_deltas(Ea, Eb, fold_of_pos, DESIGN["outer_blocks"])
            dd = (pa - pb)
            cls = classify(delta, float(hi), sum(x < 0 for x in bd), int((dd < 0).sum()))
            comps[cid]["families"][m] = {"mse_a": ma, "mse_b": mb, "delta": delta, "relative_mse_reduction": (mb - ma) / mb, "interval_95_descriptive": [float(lo), float(hi)],
                                         "blocks_improving": int(sum(x < 0 for x in bd)), "block_deltas": bd, "drugs_improving": int((dd < 0).sum()), "per_drug_delta": dict(zip(drugs, dd.tolist())), "classification": cls}
        comps[cid]["outcome"] = family_outcome(comps[cid]["families"]["ridge"]["classification"], comps[cid]["families"]["gbm"]["classification"])
    eng = engineering_outcome({c: {m: {"rel": comps[c]["families"][m]["relative_mse_reduction"], "cls": comps[c]["families"][m]["classification"]} for m in DESIGN["models"]} for c in ("C3", "C4")})
    return {"per_set_model": per, "comparisons": comps, "engineering_outcome": eng}


def correlations(drug_data, preds):
    out = {}
    from scipy.stats import spearmanr
    for key in next(iter(preds.values())):
        out["|".join(key)] = {d: {"pearson": float(np.corrcoef(preds[d][key], a["y"])[0, 1]), "spearman": float(spearmanr(preds[d][key], a["y"])[0])} for d, a in drug_data.items()}
    return out


# --- resource accounting and I/O --------------------------------------------------------------------------------------------------

def cpu_seconds():
    s, c = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime


def rss_mb():
    return max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024.0


FROZEN = ["plan.json", "evaluate.py", "../../../tests/test_met_eval.py"]


def input_hashes():
    ex = json.loads((EXPORT / "results/manifest.json").read_text())
    return {"WT_rosace_effect_all.tsv": sha256_file(HERE / "data_raw/WT_rosace_effect_all.tsv"), "met_scores_filtered.tsv": sha256_file(HERE / "data_raw/met_scores_filtered.tsv"),
            "all_features_all_data_pre_binning.csv": sha256_file(HERE / "data_raw/all_features_all_data_pre_binning.csv"), "exporter_tables": {k: v["sha256"] for k, v in ex["tables"].items()},
            "exporter_freeze": json.loads((EXPORT / "freeze_exporter.json").read_text())["files"]}


def freeze():
    rec = {"schema": "topofold.met-eval-freeze/1", "files": {f: sha256_file(HERE / f) for f in FROZEN}, "inputs": input_hashes()}
    (HERE / "freeze_eval.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


def check_freeze():
    fz = json.loads((HERE / "freeze_eval.json").read_text())
    bad = [f for f, h in fz["files"].items() if sha256_file(HERE / f) != h]
    if bad or fz["inputs"] != input_hashes():
        raise SystemExit(f"frozen file or input differs: {bad or 'inputs'}")
    return fz


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} cpu={cpu_seconds():.1f}s rss={rss_mb():.0f}MB {msg}"
    print(line, flush=True)
    with open(HERE / "logs/run.log", "a") as f:
        f.write(line + "\n")


def assemble(n_positions=None):
    me = load_export()
    drugs = DESIGN["drugs"]
    assay = read_assay(HERE / "data_raw/WT_rosace_effect_all.tsv", drugs)
    esm = read_esm(HERE / "data_raw/all_features_all_data_pre_binning.csv")
    sd = read_sd(HERE / "data_raw/met_scores_filtered.tsv", drugs)
    struct = structure_lookup(EXPORT / "results/features_structure.csv")
    mutation = mutation_lookup(EXPORT / "results/features_mutation.csv")
    ref = me.read_fasta(EXPORT / "data/P08581.fasta")
    for (inh, pos, mut), (_, _, wt) in assay.items():
        if ref[pos + 1058 - 1] != wt:
            raise ValueError(f"assay wild type differs from the reference at {pos}")
    data = {d: build_drug_arrays(d, gamma_table(assay, d), struct, mutation, esm, sd) for d in drugs}
    return data


def run_units(data, fold_of_pos, prog, cpu_cap, rss_cap, logger=print):
    """Run (drug, fold) units not yet in the progress file; stop when the CPU or memory cap is reached. Completed units are never recomputed."""
    prog = Path(prog)
    done = {}
    if prog.is_file():
        for ln in open(prog):
            r = json.loads(ln)
            done[(r["drug"], r["fold"])] = r
    status = "complete"
    for d, a in data.items():
        for k in range(DESIGN["outer_blocks"]):
            if (d, k) in done:
                continue
            if cpu_seconds() > cpu_cap or rss_mb() > rss_cap:
                return done, "budget exhausted: partial"
            preds, info, test = run_unit(a, fold_of_pos, k)
            rec = {"drug": d, "fold": k, "info": info, "row_ids": [a["row_id"][i] for i in np.where(test)[0]], "preds": {f"{s}|{m}": v.tolist() for (s, m), v in preds.items()}}
            with open(prog, "a") as f:
                f.write(json.dumps(rec) + "\n")
            done[(d, k)] = rec
            logger(f"unit {d} fold {k} {info['flag']} train={info['train_rows']} test={info['test_rows']}")
    return done, status


def run():
    fz = check_freeze()
    resource.setrlimit(resource.RLIMIT_CPU, (DESIGN["cpu_seconds_cap"] + 60, DESIGN["cpu_seconds_cap"] + 120))
    t0 = time.time()
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "logs").mkdir(exist_ok=True)
    data = assemble()
    n = DESIGN["n_positions"]
    fold_of_pos = block_assignment(n, DESIGN["outer_blocks"])
    prog = HERE / "results/progress.jsonl"
    log(f"start; rows per drug " + json.dumps({d: a["counts"] for d, a in data.items()}))
    done, status = run_units(data, fold_of_pos, prog, DESIGN["cpu_seconds_cap"], DESIGN["rss_mb_cap"], log)
    cells = [{"drug": d, "fold": k, **({"flag": done[(d, k)]["info"]["flag"], "train_rows": done[(d, k)]["info"]["train_rows"], "test_rows": done[(d, k)]["info"]["test_rows"], "alpha": done[(d, k)]["info"]["alpha"]} if (d, k) in done else {"flag": "not_run"})}
             for d in data for k in range(DESIGN["outer_blocks"])]
    failed = [c for c in cells if c["flag"] not in ("ok", "not_run")]
    out = {"schema": "topofold.met-eval-results/1", "status": status, "freeze": fz, "design": DESIGN, "cells": cells, "non_ok_cells": failed, "row_counts": {d: a["counts"] for d, a in data.items()}}
    complete = status == "complete" and not any(c["flag"] in ("no_training_rows", "not_run") for c in cells)
    if complete:
        preds = {d: {} for d in data}
        for d, a in data.items():
            idx = {r: i for i, r in enumerate(a["row_id"])}
            for key in [f"{s}|{m}" for s in DESIGN["feature_sets"] for m in DESIGN["models"]]:
                arr = np.full(len(a["y"]), np.nan)
                for k in range(DESIGN["outer_blocks"]):
                    rec = done[(d, k)]
                    if key in rec["preds"]:
                        for rid, v in zip(rec["row_ids"], rec["preds"][key]):
                            arr[idx[rid]] = v
                preds[d][tuple(key.split("|"))] = arr
            if any(np.isnan(v).any() for v in preds[d].values()):
                raise RuntimeError(f"rows without an out-of-fold prediction for {d}")
        out["evaluation"] = evaluate_predictions(data, preds, fold_of_pos, n)
        med = np.nanmedian(np.concatenate([a["sd"] for a in data.values()]))
        mask = {d: np.nan_to_num(a["sd"], nan=np.inf) <= med for d, a in data.items()}
        out["sensitivity_low_sd_rows"] = {"median_sd_gamma_approx": float(med), "rows_kept": {d: int(m.sum()) for d, m in mask.items()}, **evaluate_predictions(data, preds, fold_of_pos, n, mask)}
        out["secondary_correlations_per_drug_pooled_rows"] = correlations(data, preds)
        allsd = np.concatenate([a["sd"][~np.isnan(a["sd"])] for a in data.values()])
        out["noise_context"] = {"mean_approx_sd_gamma_squared": float(np.mean(allsd ** 2)), "caveat": "sd_gamma^2 = sd_inh^2 + sd_DMSO^2 assumes zero covariance between the two estimates, which share controls; context only, not a measurement-error model and not a detectable-improvement bound"}
        with gzip.open(HERE / "results/oof_predictions.csv.gz", "wt", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["drug", "row_id", "outer_fold", "feature_set", "model", "prediction"])
            for d, a in data.items():
                for (s, m), arr in preds[d].items():
                    for rid, p, v in zip(a["row_id"], a["pos"], arr):
                        w.writerow([d, rid, int(fold_of_pos[p - 1]), s, m, f"{v:.6g}"])
        with open(HERE / "results/folds.csv", "w") as f:
            f.write("position,outer_fold\n" + "\n".join(f"{p + 1},{int(k)}" for p, k in enumerate(fold_of_pos)) + "\n")
    out["resources"] = {"cpu_seconds_total_process_and_children": cpu_seconds(), "wall_seconds": time.time() - t0, "max_rss_mb": rss_mb(), "cpu_cap_seconds": DESIGN["cpu_seconds_cap"], "rss_cap_mb": DESIGN["rss_mb_cap"], "threads": "bounded to 1 via environment before numpy import", "child_processes_started": 0}
    (HERE / "results/results.json").write_text(json.dumps(out, indent=1))
    log(f"finished status={out['status']} non_ok_cells={len(failed)}")


if __name__ == "__main__":
    {"freeze": freeze, "run": run}[sys.argv[1]]()

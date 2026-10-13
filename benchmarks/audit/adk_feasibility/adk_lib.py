"""Pure functions for WP-B G0/G1 (SPEC.md + AMENDMENTS.md, Amendment 1). Units: Å, ps."""
from pathlib import Path

import numpy as np

THREE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
         'GLY': 'G', 'HIS': 'H', 'HSD': 'H', 'HSE': 'H', 'HSP': 'H', 'ILE': 'I', 'LEU': 'L',
         'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
         'TYR': 'Y', 'VAL': 'V'}


# --- reference parsing and mapping (A1) -------------------------------------------------------------

def parse_pdb_ca(path):
    """All CA records per chain and resSeq (every altloc kept), hetero names, REVDAT dates."""
    chains, het, revdat = {}, set(), []
    for line in Path(path).read_text().splitlines():
        rec = line[:6]
        if rec == 'REVDAT':
            revdat.append(line[13:22].strip())
        elif rec == 'HETATM':
            het.add(line[17:20].strip())
        elif rec == 'ATOM  ' and line[12:16].strip() == 'CA':
            chains.setdefault(line[21], {}).setdefault(int(line[22:26]), []).append(
                (line[16], line[17:20].strip(),
                 (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return chains, sorted(het), revdat


def select_altloc(records):
    """A1 policy: blank altloc, else 'A', else None (mapping failure)."""
    for want in (' ', 'A'):
        for alt, name, xyz in records:
            if alt == want:
                return alt, name, xyz
    return None


def map_reference(chain_records, declared_seq, resids):
    """Enforce residue correspondence; return (xyz or None, report)."""
    xyz, missing, mismatched, altlocs = [], [], [], {}
    for i, r in enumerate(resids):
        recs = chain_records.get(r)
        sel = select_altloc(recs) if recs else None
        if sel is None:
            missing.append(r)
            xyz.append((np.nan,) * 3)
            continue
        alt, name, p = sel
        if len(recs) > 1 or alt != ' ':
            altlocs[r] = {'available': sorted({a for a, _, _ in recs}), 'used': alt}
        if THREE.get(name) != declared_seq[i]:
            mismatched.append({'resseq': r, 'resname': name, 'declared': declared_seq[i]})
        xyz.append(p)
    valid = not missing and not mismatched
    return (np.array(xyz, float) if valid else None,
            {'valid': valid, 'missing': missing, 'mismatched': mismatched, 'altlocs': altlocs})


# --- geometry --------------------------------------------------------------------------------------

def kabsch(P, Q):
    """Proper rotation R and centroids with (P - cp) @ R.T + cq ≈ Q (least squares, no reflection)."""
    cp, cq = P.mean(0), Q.mean(0)
    H = (P - cp).T @ (Q - cq)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T, cp, cq


def fit_rmsd(X, Ref, fit_idx, eval_idx_list):
    R, cx, cr = kabsch(X[fit_idx], Ref[fit_idx])
    Y = (X - cx) @ R.T + cr
    return [float(np.sqrt(np.mean(np.sum((Y[i] - Ref[i]) ** 2, axis=1)))) for i in eval_idx_list]


def derive_core(A, B, cutoff, max_iter, min_size, trim_frac=0.1, max_trim=200):
    """A2 as revised by Amendment 2.

    Phase 1 (trim): fit on the current set and drop the worst ceil(trim_frac·n) residues, until
    every deviation in the set is < cutoff (at most max_trim steps).
    Phase 2 (refine): fit on the set, take all residues with deviation < cutoff, and repeat until
    the set is unchanged (at most max_iter steps).
    Converged only if phase 1 met its condition and phase 2 reached a fixed point.
    """
    idx = np.arange(len(A))
    history, trimmed, converged = [], False, False
    for _ in range(max_trim):
        if len(idx) < 3:
            break
        R, cb, ca = kabsch(B[idx], A[idx])
        dev = np.linalg.norm((B - cb) @ R.T + ca - A, axis=1)
        history.append(['trim', int(len(idx)), float(dev[idx].max())])
        if dev[idx].max() < cutoff:
            trimmed = True
            break
        k = max(1, int(np.ceil(trim_frac * len(idx))))
        idx = np.sort(idx[np.argsort(dev[idx], kind='stable')[:-k]])
    if trimmed:
        for _ in range(max_iter):
            R, cb, ca = kabsch(B[idx], A[idx])
            dev = np.linalg.norm((B - cb) @ R.T + ca - A, axis=1)
            new = np.where(dev < cutoff)[0]
            history.append(['refine', int(len(new))])
            if np.array_equal(new, idx):
                converged = True
                break
            idx = new
            if len(idx) < 3:
                break
    valid = converged and len(idx) >= min_size
    reason = None if valid else ('not converged' if not converged else f'size {len(idx)} < {min_size}')
    return idx, {'size': int(len(idx)), 'iterations': history, 'converged': converged,
                 'valid': valid, 'invalid_reason': reason}


def make_whole(X, box):
    """Sequential minimum image along the chain, orthorhombic box lengths box[:3]."""
    Y = np.array(X, float, copy=True)
    L = np.asarray(box[:3], float)
    for i in range(1, len(Y)):
        Y[i] -= L * np.round((Y[i] - Y[i - 1]) / L)
    return Y


def max_consecutive(X):
    return float(np.linalg.norm(np.diff(X, axis=0), axis=1).max())


# --- statistics ------------------------------------------------------------------------------------

def tau_int(x, c=5.0):
    """Integrated autocorrelation time (frames), Sokal windowing; an estimate."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    if acf[0] <= 0:
        return None
    rho = acf / acf[0]
    tau = 1.0
    for m in range(1, n):
        tau = 1.0 + 2.0 * rho[1:m + 1].sum()
        if m >= c * tau:
            break
    return float(tau)


# --- classification and runs (A3, A5) ----------------------------------------------------------

def classify_frames(r_open, r_closed, radius):
    """near X: rmsd_X <= radius and strictly closer to X than to the other; ties -> ambiguous."""
    r_open, r_closed = np.asarray(r_open), np.asarray(r_closed)
    near_o = (r_open <= radius) & (r_open < r_closed)
    near_c = (r_closed <= radius) & (r_closed < r_open)
    return np.where(near_o, 'open', np.where(near_c, 'closed', 'ambiguous'))


def runs_with_time(labels, times, nominal_dt, tol):
    """Maximal same-class runs; never spanning an interval that differs from nominal_dt by > tol."""
    out, start = [], 0
    for i in range(1, len(labels) + 1):
        boundary = (i == len(labels) or labels[i] != labels[i - 1]
                    or abs((times[i] - times[i - 1]) - nominal_dt) > tol)
        if boundary:
            out.append({'class': str(labels[start]), 'start': int(start), 'end': int(i - 1),
                        'n_frames': int(i - start), 'span_ps': float(times[i - 1] - times[start])})
            start = i
    return out


def summarize_runs(runs, sustained_ps):
    s = {}
    for c in ('open', 'closed', 'ambiguous'):
        rr = [r for r in runs if r['class'] == c]
        sus = [r for r in rr if r['span_ps'] >= sustained_ps]
        s[c] = {'frames': int(sum(r['n_frames'] for r in rr)), 'sustained_runs': len(sus),
                'brief_runs': len(rr) - len(sus),
                'longest_span_ns': max((r['span_ps'] for r in rr), default=0.0) / 1000.0,
                'longest_run_frames': max((r['n_frames'] for r in rr), default=0)}
    return s


def starts_near(labels, times, start_ref, window_ps):
    sel = (np.asarray(times) - times[0]) <= window_ps
    return bool(np.all(np.asarray(labels)[sel] == start_ref))


def variant_category(trajs, dominance):
    """A3 per-variant rules, first match wins. trajs: {'open_start': {...}, 'closed_start': {...}},
    each with 'summary' (summarize_runs), 'n_frames', 'start_ref' ('open'/'closed')."""
    other = {'open': 'closed', 'closed': 'open'}
    if any(t['summary']['ambiguous']['frames'] > dominance * t['n_frames'] for t in trajs.values()):
        return 'unresolved'
    if any(t['summary']['open']['sustained_runs'] > 0 and t['summary']['closed']['sustained_runs'] > 0
           for t in trajs.values()):
        return 'suitable'
    if not any(t['summary'][t['start_ref']]['sustained_runs'] > 0 for t in trajs.values()):
        return 'unsuitable'
    if not any(t['summary'][other[t['start_ref']]]['sustained_runs'] > 0 for t in trajs.values()):
        return 'restricted'
    return 'unresolved'


def final_outcome(categories):
    """A3: identical category across all valid variants, otherwise unresolved."""
    vals = set(categories.values())
    return (vals.pop() if len(vals) == 1 else 'unresolved'), len(vals) <= 1


OUTCOME_TEXT = {'suitable': 'suitable for the proposed development task',
                'restricted': 'suitable only for a restricted task',
                'unsuitable': 'unsuitable',
                'unresolved': 'unresolved pending specified evidence'}

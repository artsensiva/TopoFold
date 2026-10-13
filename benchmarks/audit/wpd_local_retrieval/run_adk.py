#!/usr/bin/env python3
"""WP-D AdK development tier: Cα-only local retrieval vs full-backbone φ/ψ target (descriptive).

Protocol: docs/WPD_LOCAL_RETRIEVAL_PROTOCOL.md (§2–§6, §10), frozen in results/wpd_protocol_freeze.json.
Inputs: WP-B raw files (Zenodo 10.5281/zenodo.5583119, CC BY 4.0), bound by SHA-256 to the WP-B G0
record. Outputs are third-party-MD derivatives (not SYNTHETIC, not TopoFold-generated data).

Stop rule (§9, §10 A2/A3): any hash, binding-identity, reference-check or target-oracle failure
raises and nothing is written as a result.

Run: build/venv/bin/python -B run_adk.py [--workers 4]
     build/venv/bin/python -B run_adk.py --smoke   (2 windows per trajectory/W, scratch output only)
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import argparse  # noqa: E402
import datetime  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import platform  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
from wpd_lib import (R_ILL, SIN_DEGENERATE, circular_mean_R, ca_pair_vectors, d_chord,  # noqa: E402
                     d_circ, d_rms, d_scalar, dihedral, horn_rmsd, kabsch_rmsd, query_metrics,
                     sin_angle, spearman, turning_angle, wrap)

WPB = HERE.parent / 'adk_feasibility'
WPA = HERE.parent / 'angular_characterization'
RAW = WPB / 'data' / 'raw'
RES = HERE / 'results'
PROTOCOL = REPO / 'docs' / 'WPD_LOCAL_RETRIEVAL_PROTOCOL.md'
FREEZE = RES / 'wpd_protocol_freeze.json'
ELIG = RES / 'wpd_adk_backbone_eligibility.json'
G0 = WPB / 'results' / 'adk_g0_provenance.json'
WPA_PROV = WPA / 'results' / 'SYNTHETIC_angular_characterization_provenance.json'
OUT = RES / 'wpd_adk_THIRD_PARTY_MD_DERIVATIVE_retrieval.json'
OUT_PROV = RES / 'wpd_adk_THIRD_PARTY_MD_DERIVATIVE_provenance.json'

PAIRS = {'open_start': ('ORADD_open_state.gro', 'ORADD_open_state_traj.xtc'),
         'closed_start': ('ORADD_closed_state.gro', 'ORADD_closed_state_traj.xtc')}
N_RES = 214
WINDOWS = (5, 8)
KS = (1, 10)
METHODS = ('M1', 'M2', 'B1', 'B2', 'B3', 'C1')
BOND_MAX_A = 2.0
REF_TOL = 1e-9
HORN_KABSCH_TOL = 1e-6
N_SPEARMAN = 20000
SEED_BASE = 20261005
BATCH = 128
N_HORN_CHECK = 10


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def split_ranges(n):
    """§5: query first ⌊0.45 n⌋, gap ⌊0.10 n⌋, database the rest; half-open, 0-based."""
    nq, ng = int(np.floor(0.45 * n)), int(np.floor(0.10 * n))
    return (0, nq), (nq, nq + ng), (nq + ng, n)


class Stop(RuntimeError):
    pass


# --- input binding (§9 stop rule, §10 A2) -------------------------------------------------------

def verify_inputs():
    freeze = json.loads(FREEZE.read_text())
    out = {'protocol_sha256': sha256(PROTOCOL), 'eligibility_sha256': sha256(ELIG),
           'freeze_sha256': sha256(FREEZE), 'wpb_g0_sha256': sha256(G0)}
    if out['protocol_sha256'] != freeze['protocol_sha256']:
        raise Stop('protocol changed after freeze')
    if out['eligibility_sha256'] != freeze['eligibility_artifact_sha256']:
        raise Stop('eligibility artifact changed after freeze')
    elig = json.loads(ELIG.read_text())
    if not elig.get('eligible') or elig['wpb_g0_sha256'] != out['wpb_g0_sha256']:
        raise Stop('eligibility artifact not eligible or bound to a different WP-B G0 record')
    g0 = json.loads(G0.read_text())
    if not g0.get('g0_passed'):
        raise Stop('WP-B G0 not passed')
    out['raw_sha256'] = {}
    for pair in PAIRS.values():
        for n in pair:
            h = sha256(RAW / n)
            if h != g0['files'][n]['sha256']:
                raise Stop(f'raw file {n} does not match WP-B G0')
            out['raw_sha256'][n] = h
    # binding identity
    import topofold
    import topofold.topofold as ext
    wpa = json.loads(WPA_PROV.read_text())
    so = Path(ext.__file__).resolve()
    so_sha = sha256(so)
    if so_sha != wpa['extension']['loaded_extension_sha256']:
        raise Stop('loaded extension differs from the WP-A extension')
    src = {}
    for rel, h in wpa['source_sha256'].items():
        if not (rel.startswith('crates/topofold-core/src/') or rel.startswith('crates/topofold-python/src/')
                or rel == 'Cargo.lock'):
            continue
        src[rel] = sha256(REPO / rel)
        if src[rel] != h:
            raise Stop(f'production source {rel} differs from the WP-A record')
    out['binding'] = {'wheel': wpa['extension']['wheel'], 'package_file': str(Path(topofold.__file__).resolve()),
                      'extension_file': str(so), 'extension_sha256': so_sha,
                      'equals_wpa_extension': True, 'production_sources_sha256': src,
                      'production_sources_equal_wpa': True}
    return out


# --- per-trajectory geometry -------------------------------------------------------------------

def _make_whole(X, box):
    """WP-B adk_lib.make_whole (sequential minimum image), imported by path to keep it identical."""
    return ADK_LIB.make_whole(X, box)


def load_traj(label):
    import MDAnalysis as mda
    gro, xtc = PAIRS[label]
    u = mda.Universe(str(RAW / gro), str(RAW / xtc))
    prot = u.select_atoms('protein')
    sel = {}
    for name in ('N', 'CA', 'C'):
        s = prot.select_atoms(f'name {name}')
        s = s[np.argsort(s.resids, kind='stable')]
        if s.n_atoms != N_RES or list(s.resids) != list(range(1, N_RES + 1)):
            raise Stop(f'{label}: atom mapping for {name} differs from the eligibility artifact')
        sel[name] = s
    F = len(u.trajectory)
    BB = np.empty((F, N_RES, 3, 3))          # N, CA, C
    CA = np.empty((F, N_RES, 3))
    times = np.empty(F)
    stats = {'ca_whole_changed_frames': 0, 'bb_whole_changed_frames': 0, 'ca_whole_max_shift_A': 0.0,
             'bb_whole_max_shift_A': 0.0, 'float32_roundtrip_changed_frames': 0, 'nonfinite_frames': 0}
    for f, ts in enumerate(u.trajectory):
        times[f] = ts.time
        box = ts.dimensions
        n, ca, c = (sel[k].positions.astype(np.float64) for k in ('N', 'CA', 'C'))
        if not (np.isfinite(n).all() and np.isfinite(ca).all() and np.isfinite(c).all()):
            raise Stop(f'{label}: non-finite coordinates in frame {f}')
        ca_w = _make_whole(ca, box)
        bb_raw = np.stack([n, ca, c], axis=1).reshape(-1, 3)   # N1 CA1 C1 N2 ...
        bb_w = _make_whole(bb_raw, box)
        sh_ca, sh_bb = float(np.abs(ca_w - ca).max()), float(np.abs(bb_w - bb_raw).max())
        stats['ca_whole_changed_frames'] += sh_ca > 0
        stats['bb_whole_changed_frames'] += sh_bb > 0
        stats['ca_whole_max_shift_A'] = max(stats['ca_whole_max_shift_A'], sh_ca)
        stats['bb_whole_max_shift_A'] = max(stats['bb_whole_max_shift_A'], sh_bb)
        ca32 = ca_w.astype(np.float32).astype(np.float64)       # A2: identical values for all methods
        stats['float32_roundtrip_changed_frames'] += bool(np.any(ca32 != ca_w))
        CA[f] = ca32
        BB[f] = bb_w.reshape(N_RES, 3, 3)
    return CA, BB, times, stats


def frame_geometry(label, CA, BB):
    """Binding τ/κ with per-frame reference check, degeneracy sines, target dihedrals, bonds, ω."""
    import topofold as tf
    F = CA.shape[0]
    tau = np.full((F, N_RES - 3), np.nan)
    kap = np.full((F, N_RES - 2), np.nan)
    binfail = np.zeros(F, bool)
    binfail_msgs = []
    sinCA = np.full((F, N_RES), np.nan)
    sinCA[:, 1:-1] = sin_angle(CA[:, :-2], CA[:, 1:-1], CA[:, 2:])
    tau_ref = dihedral(CA[:, :-3], CA[:, 1:-2], CA[:, 2:-1], CA[:, 3:])
    kap_ref = turning_angle(CA[:, :-2], CA[:, 1:-1], CA[:, 2:])
    tau_degenerate = (sinCA[:, 1:-2] < SIN_DEGENERATE) | (sinCA[:, 2:-1] < SIN_DEGENERATE)
    ref = {'max_tau_circ_diff': 0.0, 'max_kappa_diff': 0.0, 'tau_positions_skipped_degenerate':
           int(tau_degenerate.sum()), 'frames_checked': 0}
    for f in range(F):
        try:
            k_, t_, _ = tf.compute_invariants(CA[f].astype(np.float32))
        except Exception as e:  # noqa: BLE001  (A1: binding error -> frame excluded, counted)
            binfail[f] = True
            if len(binfail_msgs) < 10:
                binfail_msgs.append({'frame': f, 'error': repr(e)})
            continue
        t_, k_ = np.asarray(t_, np.float64), np.asarray(k_, np.float64)
        if t_.shape != (N_RES - 3,) or k_.shape != (N_RES - 2,):
            raise Stop(f'{label}: binding output shape {t_.shape}, {k_.shape} in frame {f}')
        dt = np.abs(wrap(t_ - tau_ref[f]))[~tau_degenerate[f]]
        dk = np.abs(k_ - kap_ref[f])
        if not (np.all(np.isfinite(dt)) and np.all(np.isfinite(dk))):
            raise Stop(f'{label}: non-finite binding/reference value in frame {f}')
        ref['max_tau_circ_diff'] = max(ref['max_tau_circ_diff'], float(dt.max(initial=0.0)))
        ref['max_kappa_diff'] = max(ref['max_kappa_diff'], float(dk.max()))
        if dt.max(initial=0.0) > REF_TOL or dk.max() > REF_TOL:
            raise Stop(f'{label}: reference check violated in frame {f}: {dt.max()}, {dk.max()}')
        tau[f], kap[f] = t_, k_
        ref['frames_checked'] += 1
    N, A, C = BB[:, :, 0], BB[:, :, 1], BB[:, :, 2]
    phi = np.full((F, N_RES), np.nan)
    psi = np.full((F, N_RES), np.nan)
    phi[:, 1:] = dihedral(C[:, :-1], N[:, 1:], A[:, 1:], C[:, 1:])
    psi[:, :-1] = dihedral(N[:, :-1], A[:, :-1], C[:, :-1], N[:, 1:])
    omega = dihedral(A[:, :-1], C[:, :-1], N[:, 1:], A[:, 1:])          # bond r -> r+1
    s_cna = np.full((F, N_RES), np.nan)
    s_cna[:, 1:] = sin_angle(C[:, :-1], N[:, 1:], A[:, 1:])
    s_nac = sin_angle(N, A, C)
    s_acn = np.full((F, N_RES), np.nan)
    s_acn[:, :-1] = sin_angle(A[:, :-1], C[:, :-1], N[:, 1:])
    b_nca = np.linalg.norm(A - N, axis=-1)
    b_cac = np.linalg.norm(C - A, axis=-1)
    b_cn = np.linalg.norm(N[:, 1:] - C[:, :-1], axis=-1)
    cis = np.abs(omega) < np.pi / 2
    geo = {'tau': tau, 'kap': kap, 'binfail': binfail, 'sinCA': sinCA, 'phi': phi, 'psi': psi,
           's_cna': s_cna, 's_nac': s_nac, 's_acn': s_acn, 'b_nca': b_nca, 'b_cac': b_cac, 'b_cn': b_cn,
           'cis': cis}
    diag = {'reference_check': ref, 'binding_failure_frames': int(binfail.sum()),
            'binding_failure_examples': binfail_msgs,
            'cis_peptide_frame_bonds': int(cis.sum()),
            'cis_peptide_bonds_residue_r_to_r_plus_1': [int(r) + 1 for r in np.where(cis.any(0))[0]],
            'min_sin_ca_vertex': float(np.nanmin(sinCA)),
            'min_sin_backbone': float(min(np.nanmin(s_cna), np.nanmin(s_nac), np.nanmin(s_acn))),
            'backbone_bond_max_A_after_whole': float(max(b_nca.max(), b_cac.max(), b_cn.max()))}
    return geo, diag


# --- one window (worker) -------------------------------------------------------------------------

G = {}


def window_validity(geo, s, W):
    """Frames excluded for this window (all methods and the target alike, A1 / §3)."""
    v = np.arange(s + 1, s + W - 1)                    # interior residue indices (0-based)
    deg_ca = (geo['sinCA'][:, v] < SIN_DEGENERATE).any(1)
    deg_t = ((geo['s_cna'][:, v] < SIN_DEGENERATE) | (geo['s_nac'][:, v] < SIN_DEGENERATE)
             | (geo['s_acn'][:, v] < SIN_DEGENERATE)).any(1)
    brk = ((geo['b_nca'][:, s:s + W] > BOND_MAX_A).any(1) | (geo['b_cac'][:, s:s + W] > BOND_MAX_A).any(1)
           | (geo['b_cn'][:, s:s + W - 1] > BOND_MAX_A).any(1))
    bf = geo['binfail']
    return ~(bf | deg_ca | deg_t | brk), {'binding_failure': bf, 'degenerate_ca': deg_ca & ~bf,
                                          'degenerate_target': deg_t & ~bf & ~deg_ca,
                                          'chain_break': brk & ~bf & ~deg_ca & ~deg_t}


def run_window(task):
    label, W, i = task
    geo, CA, (q0, q1), (d0, d1) = G[label]['geo'], G[label]['CA'], G[label]['q'], G[label]['d']
    s = i - 1
    valid, why = window_validity(geo, s, W)
    qidx = np.arange(q0, q1)[valid[q0:q1]]
    didx = np.arange(d0, d1)[valid[d0:d1]]
    nQ, nD = len(qidx), len(didx)
    rec = {'trajectory': label, 'W': W, 'i': i, 'n_query_valid': nQ, 'n_db_valid': nD,
           'excluded_query': {k: int(v[q0:q1].sum()) for k, v in why.items()},
           'excluded_db': {k: int(v[d0:d1].sum()) for k, v in why.items()}}
    if nQ == 0 or nD < max(KS):
        rec['undefined'] = 'too few valid frames'
        return rec
    tau = geo['tau'][:, s:s + W - 3]
    kap = geo['kap'][:, s:s + W - 2]
    X = CA[:, s:s + W]
    v = np.arange(s + 1, s + W - 1)
    tgt = np.concatenate([geo['phi'][:, v], geo['psi'][:, v]], axis=1)
    pairs = ca_pair_vectors(X)
    mean, Rlen = circular_mean_R(tau)
    mean_c1 = np.where(Rlen < R_ILL, np.nan, mean)
    feats = {'M1': (tau, d_chord), 'M2': (tau, d_circ), 'B1': (pairs, d_rms), 'B2': (X, horn_rmsd),
             'B3': (kap, d_rms), 'C1': (mean_c1, d_scalar)}
    Rv = Rlen[np.concatenate([qidx, didx])]
    rec['c1_conditioning'] = {
        'R_min': float(Rv.min()), 'R_p05': float(np.quantile(Rv, 0.05)), 'R_median': float(np.median(Rv)),
        'n_frames': int(len(Rv)), 'n_ill': int(np.sum(Rv < R_ILL)),
        'n_query_ill': int(np.sum(Rlen[qidx] < R_ILL)), 'n_db_ill': int(np.sum(Rlen[didx] < R_ILL))}
    rng = np.random.default_rng(SEED_BASE + 1000 * W + i)
    sq = rng.integers(0, nQ, N_SPEARMAN)
    sd = rng.integers(0, nD, N_SPEARMAN)
    samp = {m: np.empty(N_SPEARMAN) for m in METHODS + ('T',)}
    acc = {(m, k): {key: [] for key in ('valid', 'overlap', 'regret', 'norm_regret', 'norm_defined',
                                        'random_regret', 'ties')} for m in METHODS for k in KS}
    cnt = {(m, k): {'rounded': 0, 'failures': 0} for m in METHODS for k in KS}
    t_all = np.empty((nQ, nD))
    oracle_failures = 0
    for b0 in range(0, nQ, BATCH):
        qb = qidx[b0:b0 + BATCH]
        T = d_circ(tgt[qb], tgt[didx])
        t_all[b0:b0 + len(qb)] = T
        for k in KS:
            oracle_failures += query_metrics(T, T, k)['n_regret_failures']
        insamp = (sq >= b0) & (sq < b0 + len(qb))
        samp['T'][insamp] = T[sq[insamp] - b0, sd[insamp]]
        for m, (feat, fn) in feats.items():
            R = fn(feat[qb], feat[didx])
            samp[m][insamp] = R[sq[insamp] - b0, sd[insamp]]
            for k in KS:
                r = query_metrics(R, T, k)
                for key in acc[(m, k)]:
                    acc[(m, k)][key].append(r[key])
                cnt[(m, k)]['rounded'] += r['n_regret_rounded']
                cnt[(m, k)]['failures'] += r['n_regret_failures']
    if oracle_failures:
        raise Stop(f'{label} W={W} i={i}: target-oracle regret failure')
    # Horn vs Kabsch on the first sampled pairs (A2)
    hk = max(abs(samp['B2'][j] - kabsch_rmsd(X[qidx[sq[j]]], X[didx[sd[j]]])) for j in range(N_HORN_CHECK))
    if hk > HORN_KABSCH_TOL:
        raise Stop(f'{label} W={W} i={i}: Horn vs Kabsch {hk}')
    rec['horn_kabsch_max_abs_diff_A'] = float(hk)
    rec['target_spread_median'] = float(np.median(t_all))
    rec['target_quantiles'] = [float(x) for x in np.quantile(t_all, [0.05, 0.25, 0.5, 0.75, 0.95])]
    rec['random'] = {f'k{k}': {'overlap': k / nD, 'norm_regret': 1.0} for k in KS}
    rec['methods'] = {}
    for m in METHODS:
        mr = {}
        for k in KS:
            a = {key: np.concatenate(v) for key, v in acc[(m, k)].items()}
            ok = a['valid']
            nd = ok & a['norm_defined']
            mr[f'k{k}'] = {
                'n_queries_valid': int(ok.sum()), 'n_queries_undefined': int((~ok).sum()),
                'mean_overlap': float(a['overlap'][ok].mean()) if ok.any() else None,
                'mean_regret': float(a['regret'][ok].mean()) if ok.any() else None,
                'mean_norm_regret': float(a['norm_regret'][nd].mean()) if nd.any() else None,
                'n_norm_undefined': int((ok & ~a['norm_defined']).sum()),
                'mean_random_regret': float(a['random_regret'][ok].mean()) if ok.any() else None,
                'n_queries_with_boundary_ties': int(a['ties'].sum()),
                'n_regret_rounded_to_zero': cnt[(m, k)]['rounded'],
                'n_regret_numerical_failures': cnt[(m, k)]['failures']}
        fin = np.isfinite(samp[m])
        mr['spearman'] = spearman(samp[m][fin], samp['T'][fin]) if fin.sum() > 2 else None
        mr['spearman_pairs_used'] = int(fin.sum())
        mr['spearman_pairs_dropped'] = int((~fin).sum())
        rec['methods'][m] = mr
    rec['spearman_seed'] = SEED_BASE + 1000 * W + i
    return rec


def init_worker(g):
    G.update(g)


# --- aggregation ---------------------------------------------------------------------------------

def _mi(vals):
    v = np.array([x for x in vals if x is not None], float)
    if len(v) == 0:
        return {'n_windows': 0, 'median': None, 'q25': None, 'q75': None}
    q = np.quantile(v, [0.25, 0.5, 0.75])
    return {'n_windows': int(len(v)), 'median': float(q[1]), 'q25': float(q[0]), 'q75': float(q[2])}


def aggregate(windows):
    out = {}
    for label in PAIRS:
        for W in WINDOWS:
            ws = [w for w in windows if w['trajectory'] == label and w['W'] == W and 'methods' in w]
            undefined = [w['i'] for w in windows if w['trajectory'] == label and w['W'] == W
                         and 'methods' not in w]
            spread = np.array([w['target_spread_median'] for w in ws])
            order = np.argsort(spread, kind='stable')
            terc = np.empty(len(ws), int)
            for t, part in enumerate(np.array_split(order, 3)):
                terc[part] = t
            cell = {'n_windows': len(ws), 'undefined_windows': undefined,
                    'target_spread_median_across_windows': _mi(spread.tolist()),
                    'tercile_spread_bounds': [[float(spread[p].min()), float(spread[p].max())]
                                              if len(p) else None for p in np.array_split(order, 3)],
                    'random_reference': {f'k{k}': {'overlap': _mi([w['random'][f'k{k}']['overlap'] for w in ws]),
                                                   'norm_regret': 1.0,
                                                   'raw_regret': _mi([w['methods']['M1'][f'k{k}']['mean_random_regret']
                                                                      for w in ws])} for k in KS},
                    'methods': {}}
            for m in METHODS:
                cm = {}
                for k in KS:
                    sub = {}
                    for key in ('mean_overlap', 'mean_regret', 'mean_norm_regret'):
                        sub[key] = _mi([w['methods'][m][f'k{k}'][key] for w in ws])
                        sub[key + '_by_spread_tercile'] = [
                            _mi([ws[j]['methods'][m][f'k{k}'][key] for j in range(len(ws)) if terc[j] == t])
                            for t in range(3)]
                    for key in ('n_queries_valid', 'n_queries_undefined', 'n_norm_undefined',
                                'n_queries_with_boundary_ties', 'n_regret_rounded_to_zero',
                                'n_regret_numerical_failures'):
                        sub['total_' + key] = int(sum(w['methods'][m][f'k{k}'][key] for w in ws))
                    sub['windows_without_defined_norm_regret'] = int(
                        sum(w['methods'][m][f'k{k}']['mean_norm_regret'] is None for w in ws))
                    cm[f'k{k}'] = sub
                cm['spearman'] = _mi([w['methods'][m]['spearman'] for w in ws])
                cm['spearman_by_spread_tercile'] = [
                    _mi([ws[j]['methods'][m]['spearman'] for j in range(len(ws)) if terc[j] == t])
                    for t in range(3)]
                cm['spearman_pairs_dropped_total'] = int(sum(w['methods'][m]['spearman_pairs_dropped'] for w in ws))
                cell['methods'][m] = cm
            cc = [w['c1_conditioning'] for w in ws]
            cell['c1_conditioning'] = {
                'R_min': float(min(c['R_min'] for c in cc)),
                'R_p05_median_across_windows': float(np.median([c['R_p05'] for c in cc])),
                'R_median_median_across_windows': float(np.median([c['R_median'] for c in cc])),
                'ill_fraction_pooled': sum(c['n_ill'] for c in cc) / sum(c['n_frames'] for c in cc),
                'n_query_ill_total': int(sum(c['n_query_ill'] for c in cc)),
                'n_db_ill_total': int(sum(c['n_db_ill'] for c in cc)),
                'windows_with_any_ill': int(sum(c['n_ill'] > 0 for c in cc))}
            cell['coverage'] = {
                'query_frame_windows_valid': int(sum(w['n_query_valid'] for w in ws)),
                'db_frame_windows_valid': int(sum(w['n_db_valid'] for w in ws)),
                'excluded_query': {k: int(sum(w['excluded_query'][k] for w in windows
                                              if w['trajectory'] == label and w['W'] == W))
                                   for k in ws[0]['excluded_query']} if ws else {},
                'excluded_db': {k: int(sum(w['excluded_db'][k] for w in windows
                                           if w['trajectory'] == label and w['W'] == W))
                                for k in ws[0]['excluded_db']} if ws else {}}
            cell['horn_kabsch_max_abs_diff_A'] = float(max(w['horn_kabsch_max_abs_diff_A'] for w in ws))
            out[f'{label}/W{W}'] = cell
    return out


# --- main ----------------------------------------------------------------------------------------

def run_tests():
    p = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-p', 'no:cacheprovider', '-q',
                        str(HERE / 'test_wpd.py')], capture_output=True, text=True, cwd=HERE)
    return {'exit_code': p.returncode, 'summary': p.stdout.strip().splitlines()[-1] if p.stdout else ''}


def git_state():
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip()
    st = subprocess.run(['git', 'status', '--porcelain'], capture_output=True, text=True, cwd=REPO).stdout
    return {'head': head, 'status_porcelain': st.rstrip('\n').splitlines()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--smoke-out', default=None)
    args = ap.parse_args()
    t0 = time.time()
    started = datetime.datetime.now(datetime.UTC).isoformat()
    tests = run_tests()
    if tests['exit_code'] != 0:
        raise Stop(f'harness tests failed: {tests}')
    inputs = verify_inputs()
    import MDAnalysis
    import scipy
    import topofold  # noqa: F401
    g, traj_diag = {}, {}
    for label in PAIRS:
        CA, BB, times, wstats = load_traj(label)
        geo, diag = frame_geometry(label, CA, BB)
        F = CA.shape[0]
        q, gap, d = split_ranges(F)
        dt = np.diff(times)
        g[label] = {'geo': geo, 'CA': CA, 'q': q, 'd': d}
        traj_diag[label] = {'n_frames': F, 'splits_half_open': {'query': q, 'gap': gap, 'database': d},
                            'split_times_ps': {'query': [float(times[q[0]]), float(times[q[1] - 1])],
                                               'database': [float(times[d[0]]), float(times[d[1] - 1])]},
                            'cadence_ps_unique': sorted({round(float(x), 3) for x in dt}),
                            'whole_reconstruction': wstats, **diag}
        print(label, F, 'frames; geometry done', f'{time.time() - t0:.0f}s', flush=True)
    tasks = [(label, W, i) for label in PAIRS for W in WINDOWS for i in range(1, N_RES + 2 - W)]
    if args.smoke:
        tasks = [t for t in tasks if t[2] in (1, 100)]
    windows = []
    ctx = mp.get_context('fork')
    with ctx.Pool(args.workers, initializer=init_worker, initargs=(g,)) as pool:
        for n, rec in enumerate(pool.imap_unordered(run_window, tasks, chunksize=1), 1):
            windows.append(rec)
            if n % 25 == 0 or n == len(tasks):
                print(f'{n}/{len(tasks)} windows {time.time() - t0:.0f}s', flush=True)
    windows.sort(key=lambda w: (list(PAIRS).index(w['trajectory']), w['W'], w['i']))
    summary = aggregate(windows)
    meta = {'producer': 'run_adk.py', 'label': 'THIRD-PARTY MD DERIVATIVE (not SYNTHETIC)',
            'data': 'Orädd et al. 2021 AdK MD, Zenodo 10.5281/zenodo.5583119, CC BY 4.0; derived quantities only',
            'tier': 'AdK development tier: descriptive only, no decision (protocol §7); not held-out validation',
            'protocol': 'docs/WPD_LOCAL_RETRIEVAL_PROTOCOL.md', 'protocol_sha256': inputs['protocol_sha256'],
            'smoke': bool(args.smoke), 'started_utc': started,
            'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(),
            'wall_seconds': round(time.time() - t0, 1),
            'parameters': {'windows_W': WINDOWS, 'k': KS, 'batch_queries': BATCH, 'n_spearman': N_SPEARMAN,
                           'seed_rule': 'SEED_BASE + 1000*W + i', 'SEED_BASE': SEED_BASE,
                           'sin_degenerate': SIN_DEGENERATE, 'R_ill': R_ILL, 'bond_max_A': BOND_MAX_A,
                           'reference_tol_rad': REF_TOL, 'horn_kabsch_tol_A': HORN_KABSCH_TOL,
                           'horn_kabsch_pairs_per_window': N_HORN_CHECK}}
    result = {'_meta': meta, 'trajectories': traj_diag, 'summary': summary, 'windows': windows}
    prov = {'_meta': {'producer': 'run_adk.py', 'label': 'THIRD-PARTY MD DERIVATIVE provenance',
                      'smoke': bool(args.smoke)},
            'command': ' '.join([Path(sys.executable).name, '-B'] + sys.argv),
            'tests': tests, 'inputs': inputs, 'git': git_state(),
            'harness_sha256': {n: sha256(HERE / n) for n in ('wpd_lib.py', 'run_adk.py', 'test_wpd.py')},
            'wpb_adk_lib_sha256': sha256(WPB / 'adk_lib.py'),
            'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                         'mdanalysis': MDAnalysis.__version__, 'scipy': scipy.__version__,
                         'workers': args.workers},
            'result_file': None}
    if args.smoke:
        out = Path(args.smoke_out)
        out.write_text(json.dumps(result, indent=1) + '\n')
        print('smoke output', out)
        return 0
    OUT.write_text(json.dumps(result, indent=1) + '\n')
    prov['result_file'] = {'path': str(OUT.relative_to(REPO)), 'sha256': sha256(OUT)}
    OUT_PROV.write_text(json.dumps(prov, indent=1) + '\n')
    print('wrote', OUT.name, OUT_PROV.name, f'{time.time() - t0:.0f}s')
    return 0


def _load_adk_lib():
    import importlib.util
    spec = importlib.util.spec_from_file_location('wpb_adk_lib', WPB / 'adk_lib.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ADK_LIB = _load_adk_lib()

if __name__ == '__main__':
    sys.exit(main())

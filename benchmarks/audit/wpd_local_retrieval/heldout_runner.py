#!/usr/bin/env python3
"""WP-D held-out evaluation runner (protocol + Amendments 1–3). Third-party MD derivatives (mdCATH, CC BY 4.0).

Stages:
  eligibility    structural eligibility (original §4.3 with the Amendment-3 cap correction) of the downloaded
                 PRIMARY files → results/wpd_mdcath_eligibility_ledger.json. Reads no target, feature or score.
  validate-dev   full domain stage on the schema sample (development data) → results/wpd_heldout_runner_dev_validation.json
  score          authorized cohort scoring with explicit paths (validated authorization; checkpoints)
  score-heldout  `score` with the official cohort, authorization and output paths
Check classes: MANDATORY (integrity, binding identity, reference check, reader invariants) → exit status 2 on failure;
ELIGIBILITY (per-domain §4.3 criteria) → ledger entries, never silently replaced; DESCRIPTIVE → recorded only.
Run: build/venv/bin/python -B heldout_runner.py {eligibility|validate-dev}
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import datetime  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'legacy_engine'))
import heldout_decision_a2 as HD2  # noqa: E402
from heldout_dependence_r1 import classify_series  # noqa: E402
from legacy_lib import d_aligned  # noqa: E402
from mdcath_reader_r1 import MdcathDomainR1, ReaderError, make_whole_r1  # noqa: E402
from wpd_lib import (R_ILL, SIN_DEGENERATE, ca_pair_vectors, circular_mean_R, d_chord, d_circ,  # noqa: E402
                     d_rms, d_scalar, dihedral, horn_rmsd, kabsch_rmsd, query_metrics, sin_angle, spearman,
                     turning_angle, wrap)

RES = HERE / 'results'
SEL = RES / 'wpd_mdcath_selection.json'
MANIFEST = RES / 'wpd_mdcath_download_manifest.json'
LEDGER = RES / 'wpd_mdcath_eligibility_ledger.json'
AUTH = RES / 'wpd_heldout_scoring_authorization.json'
TEMPS, REPS, WINDOWS, KS = ('320', '348'), ('0', '1'), (5, 8), (1, 10)
METHODS = ('M1', 'M2', 'B1', 'B2', 'B3', 'A', 'C1')
CANDIDATES = HD2.BEST_ORDER
TARGETS = ('phipsi', 'pp')
MIN_FRAMES, MIN_RESIDUES, MIN_WINDOWS, BURN = 200, 50, 10, 0.10
REF_TOL, BOND_MAX = 1e-9, 2.0
HORN_KABSCH_TOL, N_HORN_CHECK = 1e-6, 10          # frozen: AdK protocol §10 A2
N_SPEARMAN, SEED_SPEARMAN = 20000, 20261005
WPA_EXT_SHA = '732909cc99ec88835f0b8a1f60152b7f134ae4a353d152be7bd9f048e232ba01'


class MandatoryCheckError(RuntimeError):
    pass


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def binding_identity():
    import topofold.topofold as ext
    s = sha256(ext.__file__)
    if s != WPA_EXT_SHA:
        raise MandatoryCheckError(f'binding extension differs from the WP-A extension: {s}')
    return s


# --- eligibility (original §4.3; cap correction Amendment 3 §5) -------------------------------------------------

def structural_eligibility(path, expected_size=None, expected_sha=None):
    """Integrity mismatch raises (mandatory). Other criteria give ineligibility reasons (ledger)."""
    rec = {'file': str(path)}
    if expected_size is not None:
        size, sha = Path(path).stat().st_size, sha256(path)
        rec.update({'size': size, 'sha256': sha})
        if size != expected_size or sha != expected_sha:
            raise MandatoryCheckError(f'integrity failure for {path}')
    reasons = []
    try:
        d = MdcathDomainR1(path)
    except (ReaderError, KeyError, OSError) as e:
        rec.update({'eligible': False, 'reasons': [f'reader: {e}']})
        return rec
    g = d.g
    if int(g.attrs.get('numChains', -1)) != 1:
        reasons.append('numChains != 1')
    reps = d.replicas()
    for t in TEMPS:
        for r in REPS:
            if t not in reps or r not in reps[t]:
                reasons.append(f'missing trajectory {t}/{r}')
                continue
            nF = reps[t][r]
            ds = g[t][r]['coords']
            if nF < MIN_FRAMES:
                reasons.append(f'{t}/{r}: numFrames {nF} < {MIN_FRAMES}')
            if ds.shape != (nF, int(g.attrs['numProteinAtoms']), 3):
                reasons.append(f'{t}/{r}: coords shape {ds.shape}')
            try:
                d.box_A(t, r)
            except ReaderError as e:
                reasons.append(f'{t}/{r}: {e}')
            nonfin = 0
            for f0 in range(0, ds.shape[0], 64):
                nonfin += int(np.sum(~np.isfinite(ds[f0:f0 + 64])))
            if nonfin:
                reasons.append(f'{t}/{r}: {nonfin} non-finite coordinate values')
    m = d.mapping
    if not d.mandatory_ok():
        reasons.append('missing or duplicated N/CA/C: ' + json.dumps({a: (m['missing'][a], m['duplicated'][a])
                                                                       for a in ('N', 'CA', 'C')}))
    if not m['backbone_contiguous']:
        reasons.append('backbone residues not one consecutive run')
    if len(m['backbone_resids']) < MIN_RESIDUES:
        reasons.append(f"backbone run {len(m['backbone_resids'])} < {MIN_RESIDUES} residues")
    rec.update({'eligible': not reasons, 'reasons': reasons, 'n_backbone_residues': len(m['backbone_resids']),
                'numFrames': {t: {r: reps.get(t, {}).get(r) for r in REPS} for t in TEMPS},
                'descriptive': {'residues_without_unique_O': m['missing']['O'] + m['duplicated']['O'],
                                'separate_cap_residues': sorted({x for x in d.resname if x in ('ACE', 'NME', 'NMA')})}})
    d.close()
    return rec


# --- geometry per trajectory --------------------------------------------------------------------------------------

def trajectory_geometry(dom, temp, rep):
    import topofold as tf
    bb, O, o_present = dom.backbone_r1(temp, rep)
    box = dom.box_A(temp, rep)
    whole, Ow = make_whole_r1(bb, O, box)
    nF = whole.shape[0]
    b0 = int(np.floor(BURN * nF))
    whole, Ow = whole[b0:], Ow[b0:]
    CA32 = whole[:, :, 1].astype(np.float32)               # Amendment 3 §3: round AFTER periodic reconstruction
    CA = CA32.astype(np.float64)                            # identical values for every Cα method and the binding
    F, L = CA.shape[:2]
    tau = np.full((F, L - 3), np.nan)
    kap = np.full((F, L - 2), np.nan)
    binfail = np.zeros(F, bool)
    tref = dihedral(CA[:, :-3], CA[:, 1:-2], CA[:, 2:-1], CA[:, 3:])
    kref = turning_angle(CA[:, :-2], CA[:, 1:-1], CA[:, 2:])
    sinCA = np.full((F, L), np.nan)
    sinCA[:, 1:-1] = sin_angle(CA[:, :-2], CA[:, 1:-1], CA[:, 2:])
    tdeg = (sinCA[:, 1:-2] < SIN_DEGENERATE) | (sinCA[:, 2:-1] < SIN_DEGENERATE)
    maxdt = maxdk = 0.0
    for f in range(F):
        try:
            k_, t_, _ = tf.compute_invariants(CA32[f])
        except Exception:  # noqa: BLE001  (frame excluded for all windows; counted)
            binfail[f] = True
            continue
        t_, k_ = np.asarray(t_, float), np.asarray(k_, float)
        dt = np.abs(wrap(t_ - tref[f]))[~tdeg[f]]
        dk = np.abs(k_ - kref[f])
        if not (np.all(np.isfinite(dt)) and np.all(np.isfinite(dk))):
            raise MandatoryCheckError(f'{temp}/{rep} frame {f}: non-finite binding or reference value')
        maxdt, maxdk = max(maxdt, float(dt.max(initial=0))), max(maxdk, float(dk.max()))
        tau[f], kap[f] = t_, k_
    if maxdt > REF_TOL or maxdk > REF_TOL:
        raise MandatoryCheckError(f'{temp}/{rep}: binding vs NumPy reference {maxdt}, {maxdk} > {REF_TOL}')
    N, A, C = whole[:, :, 0], whole[:, :, 1], whole[:, :, 2]
    phi = np.full((F, L), np.nan)
    psi = np.full((F, L), np.nan)
    phi[:, 1:] = dihedral(C[:, :-1], N[:, 1:], A[:, 1:], C[:, 1:])
    psi[:, :-1] = dihedral(N[:, :-1], A[:, :-1], C[:, :-1], N[:, 1:])
    s_cna = np.full((F, L), np.nan)
    s_cna[:, 1:] = sin_angle(C[:, :-1], N[:, 1:], A[:, 1:])
    s_acn = np.full((F, L), np.nan)
    s_acn[:, :-1] = sin_angle(A[:, :-1], C[:, :-1], N[:, 1:])
    return {'whole': whole, 'O': Ow, 'o_present': o_present, 'CA': CA, 'tau': tau, 'kap': kap, 'binfail': binfail,
            'sinCA': sinCA, 's_cna': s_cna, 's_nac': sin_angle(N, A, C), 's_acn': s_acn,
            'b_nca': np.linalg.norm(A - N, axis=-1), 'b_cac': np.linalg.norm(C - A, axis=-1),
            'b_cn': np.linalg.norm(N[:, 1:] - C[:, :-1], axis=-1), 'phi': phi, 'psi': psi,
            'burn_in_frames': b0, 'reference_max': [maxdt, maxdk],
            'rounding_changed_values': int(np.sum(CA32.astype(np.float64) != whole[:, :, 1]))}


def window_valid(geo, s, W):
    v = np.arange(s + 1, s + W - 1)
    deg = ((geo['sinCA'][:, v] < SIN_DEGENERATE).any(1) | (geo['s_cna'][:, v] < SIN_DEGENERATE).any(1)
           | (geo['s_nac'][:, v] < SIN_DEGENERATE).any(1) | (geo['s_acn'][:, v] < SIN_DEGENERATE).any(1))
    brk = ((geo['b_nca'][:, s:s + W] > BOND_MAX).any(1) | (geo['b_cac'][:, s:s + W] > BOND_MAX).any(1)
           | (geo['b_cn'][:, s:s + W - 1] > BOND_MAX).any(1))
    return ~(geo['binfail'] | deg | brk)


def pp_atoms(geo, s, W):
    r = np.arange(s, s + W - 1)
    w = geo['whole']
    return np.stack([w[:, r, 2], geo['O'][:, r], w[:, r + 1, 0]], axis=2).reshape(w.shape[0], 3 * (W - 1), 3)


# --- domain stage -----------------------------------------------------------------------------------------------

def score_domain(path, domain_index):
    dom = MdcathDomainR1(path)
    out = {'domain': dom.domain, 'domain_index': domain_index, 'cells': {}}
    for ti, t in enumerate(TEMPS):
        gq, gd = trajectory_geometry(dom, t, '0'), trajectory_geometry(dom, t, '1')
        L = gq['CA'].shape[1]
        out.setdefault('geometry', {})[t] = {k: {'burn_in_frames': g['burn_in_frames'], 'n_frames_used': g['CA'].shape[0],
                                                 'binding_failures': int(g['binfail'].sum()),
                                                 'reference_max_tau_kappa': g['reference_max'],
                                                 'rounding_changed_values': g['rounding_changed_values']}
                                             for k, g in (('query_rep0', gq), ('database_rep1', gd))}
        for W in WINDOWS:
            for tgt in TARGETS:
                out['cells'][f'{t}/W{W}/{tgt}'] = {'windows': [], 'per_query': []}
            for s in range(0, L - W + 1):
                vq, vd = window_valid(gq, s, W), window_valid(gd, s, W)
                qi, di = np.where(vq)[0], np.where(vd)[0]
                if len(qi) == 0 or len(di) < max(KS):
                    for tgt in TARGETS:
                        out['cells'][f'{t}/W{W}/{tgt}']['windows'].append({'s': s, 'undefined': 'too few valid frames'})
                    continue
                F = {}
                for name, g, ix in (('q', gq, qi), ('d', gd, di)):
                    tau = g['tau'][ix, s:s + W - 3]
                    kap = g['kap'][ix, s:s + W - 2]
                    X = g['CA'][ix, s:s + W]
                    mean, Rl = circular_mean_R(tau)
                    v = np.arange(s + 1, s + W - 1)
                    F[name] = {'tau': tau, 'kap': kap, 'X': X, 'pairs': ca_pair_vectors(X),
                               'c1': np.where(Rl < R_ILL, np.nan, mean),
                               'tgt': np.concatenate([g['phi'][ix][:, v], g['psi'][ix][:, v]], 1)}
                q, d = F['q'], F['d']
                R = {'M1': d_chord(q['tau'], d['tau']), 'M2': d_circ(q['tau'], d['tau']),
                     'B1': d_rms(q['pairs'], d['pairs']), 'B2': horn_rmsd(q['X'], d['X']),
                     'B3': d_rms(q['kap'], d['kap']), 'A': d_aligned(q['kap'], q['tau'], d['kap'], d['tau']),
                     'C1': d_scalar(q['c1'], d['c1'])}
                pp_ok = bool(gq['o_present'][s:s + W - 1].all() and gd['o_present'][s:s + W - 1].all())
                T = {'phipsi': d_circ(q['tgt'], d['tgt'])}
                if pp_ok:
                    PPq, PPd = pp_atoms(gq, s, W)[qi], pp_atoms(gd, s, W)[di]
                    T['pp'] = horn_rmsd(PPq, PPd)
                rng = np.random.default_rng(SEED_SPEARMAN + 1000 * W + (s + 1) + 100000 * domain_index + 10000000 * ti)
                sq, sd = rng.integers(0, len(qi), N_SPEARMAN), rng.integers(0, len(di), N_SPEARMAN)
                hk = 0.0                             # execution correction EC2: sampled Horn vs Kabsch (mandatory)
                for j in range(N_HORN_CHECK):
                    hk = max(hk, abs(R['B2'][sq[j], sd[j]] - kabsch_rmsd(q['X'][sq[j]], d['X'][sd[j]])))
                    if 'pp' in T:
                        hk = max(hk, abs(T['pp'][sq[j], sd[j]] - kabsch_rmsd(PPq[sq[j]], PPd[sd[j]])))
                if hk > HORN_KABSCH_TOL:
                    raise MandatoryCheckError(f'Horn vs Kabsch {hk} Å > {HORN_KABSCH_TOL} at {t}/W{W} s={s}')
                out.setdefault('horn_kabsch_max_A', 0.0)
                out['horn_kabsch_max_A'] = max(out['horn_kabsch_max_A'], float(hk))
                for tgt in TARGETS:
                    cell = out['cells'][f'{t}/W{W}/{tgt}']
                    if tgt not in T:
                        cell['windows'].append({'s': s, 'undefined': 'peptide-plane O unavailable'})
                        continue
                    for k in KS:                     # AdK §10 A3: oracle-side failures stop the run
                        if query_metrics(T[tgt], T[tgt], k)['n_regret_failures']:
                            raise MandatoryCheckError(f'target-oracle regret failure {t}/W{W}/{tgt} s={s}')
                    wrec = {'s': s, 'n_query': int(len(qi)), 'n_db': int(len(di)), 'methods': {},
                            'target_spread_median': float(np.median(T[tgt])),
                            'random_overlap': {f'k{k}': k / len(di) for k in KS}}
                    nrq = {}
                    for m in METHODS:
                        e = {}
                        for k in KS:
                            r = query_metrics(R[m], T[tgt], k)
                            ok = r['valid']
                            nd = ok & r['norm_defined']
                            e[f'k{k}'] = {'mean_overlap': float(r['overlap'][ok].mean()) if ok.any() else None,
                                          'mean_regret': float(r['regret'][ok].mean()) if ok.any() else None,
                                          'mean_norm_regret': float(r['norm_regret'][nd].mean()) if nd.any() else None,
                                          'n_valid': int(ok.sum()), 'n_norm_defined': int(nd.sum()),
                                          'n_regret_failures': r['n_regret_failures']}
                            if k == 10:
                                full = np.full(len(vq), np.nan)
                                full[qi] = np.where(nd, r['norm_regret'], np.nan)
                                nrq[m] = full
                        a, b = R[m][sq, sd], T[tgt][sq, sd]
                        fin = np.isfinite(a)
                        e['spearman'] = spearman(a[fin], b[fin]) if fin.sum() > 2 else None
                        wrec['methods'][m] = e
                    c1ok = np.isfinite(nrq['C1'])
                    wrec['c1_matched_means'] = ({m: float(np.nanmean(nrq[m][c1ok])) for m in METHODS}
                                                if c1ok.any() and all(np.isfinite(nrq[m][c1ok]).all() for m in METHODS)
                                                else None)
                    cell['windows'].append(wrec)
                    cell['per_query'].append(nrq)
    dom.close()
    return summarize_domain(out)


def summarize_domain(out):
    """Per cell: per-method median over windows (≥ MIN_WINDOWS defined for every non-C1 method), coverage-matched
    medians, and dependence diagnostics on per-query window-averaged contrast series (NaN kept in place)."""
    for key, cell in out['cells'].items():
        ws = [w for w in cell['windows'] if 'methods' in w]
        vals = {m: [w['methods'][m]['k10']['mean_norm_regret'] for w in ws] for m in METHODS}
        nd = {m: sum(v is not None for v in vals[m]) for m in METHODS}
        defined = all(nd[m] >= MIN_WINDOWS for m in METHODS if m != 'C1')
        cell['n_windows_defined'] = nd
        cell['domain_nr10'] = ({m: float(np.median([v for v in vals[m] if v is not None])) if nd[m] else None
                                for m in METHODS} if defined else None)
        mt = [w['c1_matched_means'] for w in ws if w.get('c1_matched_means')]
        cell['n_matched_windows'] = len(mt)
        cell['domain_nr10_c1_matched'] = ({m: float(np.median([x[m] for x in mt])) for m in METHODS}
                                          if len(mt) >= MIN_WINDOWS and defined else None)
        pq = cell.pop('per_query')
        dep = {}
        if pq:
            for i, x in enumerate(CANDIDATES):
                for y in CANDIDATES[i + 1:]:
                    S = np.stack([p[x] - p[y] for p in pq])
                    with np.errstate(all='ignore'):
                        cnt = np.isfinite(S).sum(0)
                        ser = np.where(cnt > 0, np.nansum(np.where(np.isfinite(S), S, 0.0), 0) / np.maximum(cnt, 1),
                                       np.nan)
                    dep[f'{x}-{y}'] = classify_series(ser)
        cell['dependence'] = dep
        cell['windows_undefined'] = [w for w in cell['windows'] if 'methods' not in w]
    return out


# --- cell stage (never run on held-out data in this task) --------------------------------------------------------

def analyse_cells(domains, target, n_boot=2000):
    res = {}
    for t in TEMPS:
        for W in WINDOWS:
            key = f'{t}/W{W}/{target}'
            rows = [d['cells'][key] for d in domains if d['cells'][key]['domain_nr10'] is not None]
            if len(rows) < HD2._mod.MIN_DOMAINS:
                res[key] = {'status': 'insufficient_evidence', 'n_defined_domains': len(rows)}
                continue
            nr = {m: np.array([r['domain_nr10'][m] if r['domain_nr10'][m] is not None else np.nan for r in rows])
                  for m in METHODS}
            nr_cm = {m: np.array([r['domain_nr10_c1_matched'][m] if r['domain_nr10_c1_matched'] else np.nan
                                  for r in rows]) for m in METHODS}
            res[key] = HD2.analyse_cell(nr, nr_cm, n_boot=n_boot)
            res[key]['n_defined_domains'] = len(rows)
    overall = {}
    for m in METHODS:
        cells = [res[k] for k in res]
        if any('methods' not in c for c in cells):
            overall[m] = 'insufficient_evidence'
            continue
        overall[m] = HD2.overall_outcome([(c['methods'][m]['loss_class'], c['methods'][m]['gain_class']) for c in cells])
    return {'cells': res, 'overall': overall}


# --- execution layer (execution corrections EC1–EC4) ----------------------------------------------------------

PROTOCOL_DOCS = ('docs/WPD_LOCAL_RETRIEVAL_PROTOCOL.md', 'docs/WPD_HELDOUT_MDCATH_PROTOCOL.md',
                 'docs/WPD_HELDOUT_MDCATH_AMENDMENT_1.md', 'docs/WPD_HELDOUT_MDCATH_AMENDMENT_2.md',
                 'docs/WPD_HELDOUT_MDCATH_AMENDMENT_3.md', 'docs/WPD_HELDOUT_MDCATH_EXECUTION_CORRECTIONS.md')
IMPLEMENTATION = (HERE / 'heldout_runner.py', HERE / 'mdcath_reader_r1.py', HERE / 'mdcath_reader.py',
                  HERE / 'heldout_dependence_r1.py', HERE / 'heldout_dependence.py', HERE / 'heldout_decision.py',
                  HERE / 'heldout_decision_a2.py', HERE / 'wpd_lib.py', HERE.parent / 'legacy_engine' / 'legacy_lib.py',
                  HERE.parent / 'adk_feasibility' / 'adk_lib.py')
WPA_PROV = HERE.parent / 'angular_characterization' / 'results' / 'SYNTHETIC_angular_characterization_provenance.json'
AUTH_REQUIRED = ('authorized_by', 'scope', 'protocol_sha256', 'selection_sha256', 'cohort_sha256',
                 'implementation_sha256', 'cohort_domains')


class AuthorizationError(RuntimeError):
    pass


def _rel(p):
    return str(Path(p).resolve().relative_to(REPO))


def implementation_hashes():
    return {_rel(p): sha256(p) for p in IMPLEMENTATION}


def protocol_hashes():
    return {d: sha256(REPO / d) if (REPO / d).exists() else None for d in PROTOCOL_DOCS}


def production_sources():
    """Mandatory: production sources compiled into the binding equal the WP-A record (as run_adk.verify_inputs)."""
    wpa = json.loads(WPA_PROV.read_text())
    out = {}
    for rel, h in wpa['source_sha256'].items():
        if rel.startswith(('crates/topofold-core/src/', 'crates/topofold-python/src/')) or rel == 'Cargo.lock':
            out[rel] = sha256(REPO / rel)
            if out[rel] != h:
                raise MandatoryCheckError(f'production source {rel} differs from the WP-A record')
    return out


def expected_authorization(selection_path, cohort_path):
    cohort = json.loads(Path(cohort_path).read_text())
    return {'protocol_sha256': protocol_hashes(), 'selection_sha256': sha256(selection_path),
            'cohort_sha256': sha256(cohort_path), 'implementation_sha256': implementation_hashes(),
            'cohort_domains': [d['domain'] for d in cohort['domains']]}


def validate_authorization(auth_path, selection_path, cohort_path):
    auth_path = Path(auth_path)
    if not auth_path.exists():
        raise AuthorizationError('authorization record absent')
    try:
        auth = json.loads(auth_path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        raise AuthorizationError(f'authorization record malformed: {e}') from e
    if not isinstance(auth, dict) or any(k not in auth for k in AUTH_REQUIRED):
        raise AuthorizationError('authorization record malformed: missing required fields')
    exp = expected_authorization(selection_path, cohort_path)
    bad = [k for k, v in exp.items() if auth.get(k) != v]
    if bad:
        raise AuthorizationError(f'authorization record does not match current files: {bad}')
    return auth


def _checkpoint_key(entry):
    return {'domain': entry['domain'], 'domain_index': entry['order'], 'input_sha256': entry['sha256'],
            'implementation_sha256': implementation_hashes(), 'protocol_sha256': protocol_hashes()}


def _score_entry(args):
    entry, path, ckpt_dir = args
    key = _checkpoint_key(entry)
    ck = Path(ckpt_dir) / f"{entry['domain']}.json"
    if ck.exists():
        try:
            prev = json.loads(ck.read_text())
            if prev.get('key') == key:
                return entry['domain'], prev['result'], True
        except json.JSONDecodeError:
            pass
    res = json.loads(json.dumps(score_domain(path, entry['order']), default=float))
    tmp = ck.with_suffix('.tmp')
    tmp.write_text(json.dumps({'key': key, 'written_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                               'result': res}) + '\n')
    os.replace(tmp, ck)
    return entry['domain'], res, False


def _median(vals):
    v = [x for x in vals if x is not None]
    return float(np.median(v)) if v else None


def secondary(domains):
    """Reported, not deciding (original §7.6): domain medians of secondary metrics, coverage, target-spread
    terciles (windows pooled across domains), dependence-class counts."""
    out = {}
    keys = list(domains[0]['cells']) if domains else []
    for key in keys:
        cells = [d['cells'][key] for d in domains]
        per_dom = {m: {k: [] for k in ('overlap10', 'nr1', 'regret10', 'spearman')} for m in METHODS}
        rnd, wins = [], []
        for c in cells:
            ws = [w for w in c['windows'] if 'methods' in w]
            if not ws or c['domain_nr10'] is None:
                continue
            rnd.append(float(np.median([w['random_overlap']['k10'] for w in ws])))
            for m in METHODS:
                per_dom[m]['overlap10'].append(_median([w['methods'][m]['k10']['mean_overlap'] for w in ws]))
                per_dom[m]['nr1'].append(_median([w['methods'][m]['k1']['mean_norm_regret'] for w in ws]))
                per_dom[m]['regret10'].append(_median([w['methods'][m]['k10']['mean_regret'] for w in ws]))
                per_dom[m]['spearman'].append(_median([w['methods'][m]['spearman'] for w in ws]))
            wins += ws
        terc = None
        if len(wins) >= 3:
            sp = np.array([w['target_spread_median'] for w in wins])
            parts = np.array_split(np.argsort(sp, kind='stable'), 3)
            terc = {'bounds': [[float(sp[p].min()), float(sp[p].max())] for p in parts],
                    'mean_window_nr10': {m: [_mean([wins[j]['methods'][m]['k10']['mean_norm_regret'] for j in p])
                                             for p in parts] for m in METHODS}}
        dep = {}
        for c in cells:
            for pair, e in c['dependence'].items():
                dep.setdefault(pair, {}).setdefault(e['class'], 0)
                dep[pair][e['class']] += 1
        undefined = {}
        for c in cells:
            for w in c['windows_undefined']:
                undefined[w['undefined']] = undefined.get(w['undefined'], 0) + 1
        out[key] = {'n_domains': len(cells), 'n_defined_domains': sum(c['domain_nr10'] is not None for c in cells),
                    'n_c1_matched_domains': sum(c['domain_nr10_c1_matched'] is not None for c in cells),
                    'windows_defined_total': sum(len([w for w in c['windows'] if 'methods' in w]) for c in cells),
                    'windows_undefined_by_reason': undefined,
                    'random_overlap10_median': _median(rnd),
                    'domain_medians': {m: {k: _median(v) for k, v in per_dom[m].items()} for m in METHODS},
                    'target_spread_terciles': terc, 'dependence_class_counts': dep}
    return out


def _mean(vals):
    v = [x for x in vals if x is not None]
    return float(np.mean(v)) if v else None


def run_cohort(cohort_path, auth_path, selection_path, data_dir, out_path, ckpt_dir, workers=4, n_boot=2000):
    """Authorized cohort scoring. Raises AuthorizationError (exit 3) or MandatoryCheckError (exit 2)."""
    import time
    t0 = time.time()
    auth = validate_authorization(auth_path, selection_path, cohort_path)
    ext = binding_identity()
    prod = production_sources()
    cohort = json.loads(Path(cohort_path).read_text())
    entries = cohort['domains']
    for e in entries:
        p = Path(data_dir) / f"mdcath_dataset_{e['domain']}.h5"
        if not p.exists() or p.stat().st_size != e['size_bytes'] or sha256(p) != e['sha256']:
            raise MandatoryCheckError(f"input integrity failure for {e['domain']}")
    Path(ckpt_dir).mkdir(parents=True, exist_ok=True)
    tasks = [(e, Path(data_dir) / f"mdcath_dataset_{e['domain']}.h5", ckpt_dir) for e in entries]
    results, reused = {}, []
    if workers <= 1:
        for t in tasks:
            dom, res, r = _score_entry(t)
            results[dom] = res
            if r:
                reused.append(dom)
            print('domain', dom, 'reused' if r else 'scored', f'{time.time() - t0:.0f}s', flush=True)
    else:
        import multiprocessing as mp
        pool = mp.get_context('fork').Pool(workers)
        try:
            for dom, res, r in pool.imap_unordered(_score_entry, tasks):
                results[dom] = res
                if r:
                    reused.append(dom)
                print('domain', dom, 'reused' if r else 'scored', f'{time.time() - t0:.0f}s', flush=True)
        except BaseException:
            pool.terminate()
            raise
        finally:
            pool.close()
            pool.join()
    domains = [results[e['domain']] for e in entries]
    agg = {tgt: analyse_cells(domains, tgt, n_boot=n_boot) for tgt in TARGETS}
    out = {'_meta': {'producer': 'heldout_runner.py score', 'label': 'THIRD-PARTY MD DERIVATIVE (mdCATH, CC BY 4.0)',
                     'authorization_sha256': sha256(auth_path), 'authorized_by': auth['authorized_by'],
                     'scope': auth['scope'], 'binding_extension_sha256': ext, 'production_sources_sha256': prod,
                     'implementation_sha256': implementation_hashes(), 'protocol_sha256': protocol_hashes(),
                     'cohort_sha256': sha256(cohort_path), 'n_domains': len(entries),
                     'reused_checkpoints': sorted(reused), 'workers': workers, 'n_boot': n_boot,
                     'wall_seconds': round(time.time() - t0, 1),
                     'finished_utc': datetime.datetime.now(datetime.UTC).isoformat()},
           'aggregate': agg, 'secondary': secondary(domains), 'domains': domains}
    Path(out_path).write_text(json.dumps(out, indent=1, default=float) + '\n')
    return out


def cmd_score(argv):
    import argparse
    ap = argparse.ArgumentParser(prog='heldout_runner.py score')
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--authorization', required=True)
    ap.add_argument('--selection', default=str(SEL))
    ap.add_argument('--data', default=str(HERE / 'data' / 'mdcath'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--checkpoints', required=True)
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args(argv)
    try:
        run_cohort(a.cohort, a.authorization, a.selection, a.data, a.out, a.checkpoints, a.workers)
    except AuthorizationError as e:
        print(f'AUTHORIZATION REJECTED: {e}')
        return 3
    except MandatoryCheckError as e:
        fail = Path(a.out).with_suffix('.FAILED.json')
        fail.write_text(json.dumps({'mandatory_failure': str(e), 'utc': datetime.datetime.now(datetime.UTC).isoformat(),
                                    'checkpoints_preserved_in': a.checkpoints}, indent=1) + '\n')
        print(f'MANDATORY CHECK FAILED: {e}')
        return 2
    print('score written', a.out)
    return 0


# --- CLI --------------------------------------------------------------------------------------------------------

def cmd_eligibility():
    sel = json.loads(SEL.read_text())
    man = json.loads(MANIFEST.read_text())
    verified = {r['domain']: r for r in man['files']}
    ledger = {'_meta': {'producer': 'heldout_runner.py eligibility', 'run_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                        'criteria': 'original §4.3 with Amendment 3 §5 cap correction', 'reads': 'structure only; no '
                        'target, feature, ranking or score', 'binding_extension_sha256': binding_identity()},
              'domains': []}
    for f in [x for x in sel['files'] if x['role'] == 'primary']:
        p = HERE / 'data' / 'mdcath' / f"mdcath_dataset_{f['domain']}.h5"
        v = verified.get(f['domain'])
        if not p.exists() or not v or not v.get('verified'):
            ledger['domains'].append({'order': f['order'], 'domain': f['domain'], 'eligible': None,
                                      'status': 'not downloaded or not verified'})
            continue
        rec = structural_eligibility(p, f['size_bytes'], f['lfs_sha256'])
        rec.update({'order': f['order'], 'domain': f['domain'], 'superfamily': f['superfamily']})
        rec['file'] = str(p.relative_to(HERE))
        ledger['domains'].append(rec)
        print(f['order'], f['domain'], 'eligible' if rec['eligible'] else f"INELIGIBLE {rec['reasons']}", flush=True)
    el = [d for d in ledger['domains'] if d.get('eligible')]
    ledger['_meta'].update({'n_primary': 30, 'n_checked': sum(d.get('eligible') is not None for d in ledger['domains']),
                            'n_eligible': len(el), 'ineligible': [d['domain'] for d in ledger['domains']
                                                                  if d.get('eligible') is False],
                            'replacements': 'none performed; any replacement requires documented failure and '
                                            'authorization'})
    LEDGER.write_text(json.dumps(ledger, indent=1) + '\n')
    return 0


def cmd_validate_dev():
    sel = json.loads(SEL.read_text())
    f = sel['files'][0]
    assert f['role'] == 'schema_sample_development'
    p = HERE / 'data' / 'mdcath' / f"mdcath_dataset_{f['domain']}.h5"
    ext = binding_identity()
    elig = structural_eligibility(p, f['size_bytes'], f['lfs_sha256'])
    import time
    t0 = time.time()
    dres = score_domain(p, domain_index=0) if elig['eligible'] else None
    out = {'_meta': {'producer': 'heldout_runner.py validate-dev', 'label': 'DEVELOPMENT DATA (schema sample, '
                     'superfamily excluded from held-out); third-party MD derivative', 'binding_extension_sha256': ext,
                     'run_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'wall_seconds': round(time.time() - t0, 1)},
           'eligibility': elig, 'domain_result': dres}
    (RES / 'wpd_heldout_runner_dev_validation.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('dev validation written; eligible', elig['eligible'])
    return 0


def main(argv):
    try:
        if argv[:1] == ['eligibility']:
            return cmd_eligibility()
        if argv[:1] == ['validate-dev']:
            return cmd_validate_dev()
        if argv[:1] == ['score']:
            return cmd_score(argv[1:])
        if argv[:1] == ['score-heldout']:
            return cmd_score(['--cohort', str(RES / 'wpd_mdcath_cohort.json'), '--authorization', str(AUTH),
                              '--out', str(RES / 'wpd_heldout_THIRD_PARTY_MD_DERIVATIVE_results.json'),
                              '--checkpoints', str(HERE / 'build' / 'heldout_checkpoints')] + argv[1:])
        print(__doc__)
        return 1
    except MandatoryCheckError as e:
        print(f'MANDATORY CHECK FAILED: {e}')
        return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

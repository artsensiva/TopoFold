"""WP-D held-out decision statistics, as amended (docs/WPD_HELDOUT_MDCATH_AMENDMENT_1.md).

Input per cell (temperature × W): per-domain NR@10 values (median over windows, protocol §6.3).
- nr[m]: primary values for m in M1, M2, B1, B2, B3 (and C1, reported only), shape (n_domains,).
- nr_cm[m]: coverage-matched values (C1 query population, amendment §3) for every method including C1;
  NaN for domains outside the common C1 population. The NaN pattern must be identical across methods.
Harness code, not production. No mdCATH data has been read when this module was written.
"""
import numpy as np

BEST_ORDER = ('B2', 'B1', 'B3', 'M1', 'M2')     # candidates for best; order breaks exact ties only
C1 = 'C1'
RHO = 0.10                 # acceptable-loss margin and informativeness threshold (original §7.3, unchanged)
EPS_G = 1e-6               # amendment §2: best-method gain at or below this is not a valid denominator
TIE = 1e-12
N_BOOT = 2000
SEED = 20261008
ALPHA = 0.05
MIN_DOMAINS = 24
BURN = 0.10


def post_burn_frames(n):
    """Frames left after excluding the first ⌊0.10 n⌋ (protocol §5)."""
    return n - int(np.floor(BURN * n))


def select_best(med):
    """argmin of median NR over BEST_ORDER; |Δ| ≤ TIE resolved by BEST_ORDER. C1 is never a candidate."""
    best = BEST_ORDER[0]
    for m in BEST_ORDER[1:]:
        if med[m] < med[best] - TIE:
            best = m
    return best


def cell_statistics(nr, nr_cm, idx):
    """Point statistics for the domain multiset idx (observed sample: idx = arange(n))."""
    med = {m: float(np.median(nr[m][idx])) for m in BEST_ORDER}
    b = select_best(med)
    gb = 1.0 - med[b]
    valid = gb > EPS_G
    out = {'best': b, 'median_nr': med, 'gain_best': gb, 'loss_valid': bool(valid),
           'G': {m: 1.0 - med[m] for m in BEST_ORDER},
           'L': {m: ((med[m] - med[b]) / gb if valid else None) for m in BEST_ORDER}}
    sub = idx[~np.isnan(nr_cm[C1][idx])]
    out['c1_n_domains'] = int(len(sub))
    if len(sub) == 0:
        out.update({'c1_loss_valid': False, 'L_C1': None, 'G_C1': None, 'gain_best_cm': None})
        return out
    mc = float(np.median(nr_cm[C1][sub]))
    mb = float(np.median(nr_cm[b][sub]))
    gbc = 1.0 - mb
    out.update({'gain_best_cm': gbc, 'G_C1': 1.0 - mc, 'c1_loss_valid': bool(gbc > EPS_G),
                'L_C1': (mc - mb) / gbc if gbc > EPS_G else None})
    return out


def conservative_interval(values, invalid, alpha=ALPHA):
    """Percentile interval (order statistics, numpy 'inverted_cdf') in which every invalid resample is
    placed at −∞ for the lower bound and at +∞ for the upper bound. Nothing is discarded; if more than
    alpha/2 of resamples are invalid, the corresponding bound is infinite."""
    v = np.asarray(values, float)
    inv = np.asarray(invalid, bool)
    lo = np.quantile(np.where(inv, -np.inf, v), alpha / 2, method='inverted_cdf')
    hi = np.quantile(np.where(inv, np.inf, v), 1 - alpha / 2, method='inverted_cdf')
    return float(lo), float(hi), float(inv.mean())


def check_inputs(nr, nr_cm):
    n = len(nr[BEST_ORDER[0]])
    for m in BEST_ORDER + (C1,):
        if len(nr_cm[m]) != n or (m in nr and len(nr[m]) != n):
            raise ValueError('length mismatch')
    pattern = np.isnan(nr_cm[C1])
    for m in BEST_ORDER:
        if not np.array_equal(np.isnan(nr_cm[m]), pattern):
            raise ValueError('coverage-matched population differs between methods')
        if np.isnan(nr[m]).any():
            raise ValueError('primary values must be defined for every domain in the cell')
    return n


def bootstrap(nr, nr_cm, n_boot=N_BOOT, seed=SEED):
    n = check_inputs(nr, nr_cm)
    rng = np.random.default_rng(seed)
    L = {m: np.empty(n_boot) for m in BEST_ORDER + (C1,)}
    G = {m: np.empty(n_boot) for m in BEST_ORDER + (C1,)}
    inv_L = np.zeros(n_boot, bool)
    inv_L1 = np.zeros(n_boot, bool)
    inv_G1 = np.zeros(n_boot, bool)
    best = []
    for r in range(n_boot):
        idx = rng.integers(0, n, n)                      # paired: same domains for every method
        s = cell_statistics(nr, nr_cm, idx)
        best.append(s['best'])
        inv_L[r] = not s['loss_valid']
        for m in BEST_ORDER:
            L[m][r] = s['L'][m] if s['loss_valid'] else np.nan
            G[m][r] = s['G'][m]
        inv_L1[r] = not s['c1_loss_valid']
        inv_G1[r] = s['G_C1'] is None
        L[C1][r] = s['L_C1'] if s['c1_loss_valid'] else np.nan
        G[C1][r] = s['G_C1'] if s['G_C1'] is not None else np.nan
    return {'L': L, 'G': G, 'invalid_L': inv_L, 'invalid_L_C1': inv_L1, 'invalid_G_C1': inv_G1,
            'best_counts': {m: best.count(m) for m in BEST_ORDER}}


def classify_loss(observed_valid, lo, hi, rho=RHO):
    if not observed_valid:
        return 'not_assessable'
    if hi <= rho:
        return 'within_margin'
    if lo > rho:
        return 'beyond_margin'
    return 'inconclusive'


def classify_gain(observed_defined, lo, hi, rho=RHO):
    if not observed_defined:
        return 'not_assessable'
    if lo > rho:
        return 'informative'
    if hi <= rho:
        return 'uninformative'
    return 'inconclusive'


def analyse_cell(nr, nr_cm, n_boot=N_BOOT, seed=SEED):
    n = check_inputs(nr, nr_cm)
    obs = cell_statistics(nr, nr_cm, np.arange(n))
    bs = bootstrap(nr, nr_cm, n_boot, seed)
    res = {'n_domains': n, 'observed': obs, 'best_counts': bs['best_counts'], 'methods': {}}
    for m in BEST_ORDER:
        lo, hi, fi = conservative_interval(bs['L'][m], bs['invalid_L'])
        glo, ghi, _ = conservative_interval(bs['G'][m], np.zeros(n_boot, bool))
        res['methods'][m] = {'L': obs['L'][m], 'L_interval': [lo, hi], 'invalid_fraction': fi,
                             'G': obs['G'][m], 'G_interval': [glo, ghi],
                             'loss_class': classify_loss(obs['loss_valid'], lo, hi),
                             'gain_class': classify_gain(True, glo, ghi)}
    c1_ok = obs['c1_n_domains'] >= MIN_DOMAINS
    lo, hi, fi = conservative_interval(bs['L'][C1], bs['invalid_L_C1'])
    glo, ghi, gfi = conservative_interval(bs['G'][C1], bs['invalid_G_C1'])
    res['methods'][C1] = {'population': 'coverage-matched', 'n_domains': obs['c1_n_domains'],
                          'L': obs['L_C1'], 'L_interval': [lo, hi], 'invalid_fraction': fi,
                          'G': obs['G_C1'], 'G_interval': [glo, ghi], 'G_invalid_fraction': gfi,
                          'loss_class': classify_loss(c1_ok and obs['c1_loss_valid'], lo, hi),
                          'gain_class': classify_gain(c1_ok and obs['G_C1'] is not None, glo, ghi)}
    return res


def overall_outcome(cells):
    """cells: list of (loss_class, gain_class), one per cell, for one method (amendment §4)."""
    ok = [lc == 'within_margin' and gc == 'informative' for lc, gc in cells]
    if all(ok):
        return 'meets_criterion_in_all_cells'
    if all(lc == 'beyond_margin' for lc, _ in cells) or all(gc == 'uninformative' for _, gc in cells):
        return 'fails_criterion_in_all_cells'
    if any(ok) and any(lc == 'beyond_margin' for lc, _ in cells):
        return 'cell_dependent'
    return 'inconclusive'

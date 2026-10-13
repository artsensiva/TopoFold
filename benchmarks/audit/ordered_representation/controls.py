"""New SYNTHETIC controls N1 (paired ±ε) and N2 (mirror), protocol §4.

N1: one permutation and one per-frame noise draw, from default_rng(20261104), shared by all six
ε values, so ±ε cases differ only in ε. N2: F1a_s0 with every torsion negated.
"""
import numpy as np

from inputs import f1_labels, load_case, wpa_build, wpa_common

N1_SEED = 20261104
N1_EPS = [1e-6, -1e-6, 1e-3, -1e-3, 1e-1, -1e-1]
N1_F, N1_P, N1_SIGMA = 1000, 0.49, np.deg2rad(1.0)


def efmt(x):
    return f'{x:+.0e}'


def build_controls():
    """Return {case_id: dict(tors, coords, labels, params, draw_order)}."""
    out = {}
    rng = np.random.default_rng(N1_SEED)
    perm = rng.permutation(N1_F)
    n2 = int(round(N1_P * N1_F))
    state = np.zeros(N1_F, dtype=bool)
    state[perm[:n2]] = True
    nu = rng.normal(0.0, N1_SIGMA, size=N1_F)
    draws = [f'default_rng({N1_SEED})', f'permutation({N1_F}); first {n2} -> phi2',
             f'normal(0, {N1_SIGMA!r}, size={N1_F}) per frame, added to both torsions',
             'same permutation and noise reused for every eps']
    for eps in N1_EPS:
        x = np.where(state, np.pi - eps, 0.0) + nu
        tors = wpa_common.wrap(np.stack([x, x], axis=1))
        out[f'N1_e{efmt(eps)}'] = dict(tors=tors, labels=state.astype(int), W=5,
                                       params=dict(p=N1_P, eps_rad=eps, sigma_rad=N1_SIGMA,
                                                   F=N1_F, n_phi2=n2),
                                       draw_order=draws)
    a_tors, _, _ = load_case('F1a_s0')
    out['N2_mirror_F1a_s0'] = dict(tors=wpa_common.wrap(-a_tors), labels=f1_labels('F1a_s0'), W=5,
                                   params=dict(derived_from='F1a_s0', transform='tau -> -tau'),
                                   draw_order=['deterministic: negation of F1a_s0 torsions'])
    for c in out.values():
        F, W = c['tors'].shape[0], c['W']
        c['coords'] = wpa_build.build_traces(np.full((F, W - 2), wpa_common.KAPPA), c['tors'])
    return out

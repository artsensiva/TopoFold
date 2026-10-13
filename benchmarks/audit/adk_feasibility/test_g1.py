"""Focused tests for WP-B G1 logic (Amendment 1, A4). Pure functions only; no trajectory data.

Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_g1.py
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adk_lib import (classify_frames, derive_core, final_outcome, fit_rmsd, kabsch,  # noqa: E402
                     make_whole, map_reference, max_consecutive, runs_with_time, select_altloc,
                     starts_near, summarize_runs, variant_category)


def rot(axis, ang):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K


# --- alignment --------------------------------------------------------------------------------

def test_kabsch_recovers_rigid_motion():
    rng = np.random.default_rng(0)
    P = rng.normal(size=(50, 3)) * 10
    Q = P @ rot([1, 2, 3], 0.9).T + np.array([5.0, -3.0, 12.0])
    assert fit_rmsd(P, Q, np.arange(50), [np.arange(50)])[0] < 1e-12
    R, _, _ = kabsch(P, Q)
    assert abs(np.linalg.det(R) - 1) < 1e-12


def test_kabsch_does_not_reflect():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(30, 3)) * 5
    assert fit_rmsd(P, P * np.array([1, 1, -1]), np.arange(30), [np.arange(30)])[0] > 0.5


# --- periodic reconstruction ---------------------------------------------------------------------

def chain(n=20):
    t = np.arange(n)
    return np.stack([3.8 * t * 0.9, 1.5 * np.sin(t), 1.5 * np.cos(t)], 1) + 10.0


def test_split_chain_is_made_whole():
    X = chain()
    box = np.array([30.0, 30.0, 30.0, 90, 90, 90])
    split = X.copy()
    split[10:] -= np.array([30.0, 0, 0])            # second half wrapped by one box length
    assert max_consecutive(split) > 4.2
    Y = make_whole(split, box)
    assert max_consecutive(Y) < 4.2
    d = Y - X
    assert np.allclose(d, d[0])                     # whole chain = original up to a lattice shift


def test_whole_chain_unchanged():
    X = chain()
    assert np.array_equal(make_whole(X, np.array([100.0, 100.0, 100.0, 90, 90, 90])), X)


# --- core convergence -----------------------------------------------------------------------------

def test_core_converges_to_rigid_part():
    rng = np.random.default_rng(2)
    A = rng.normal(size=(60, 3)) * 8
    B = A @ rot([0, 0, 1], 0.4).T + 3.0
    B[40:] += np.array([6.0, 0, 0])                  # moved segment
    idx, meta = derive_core(A, B, 1.0, 50, 30)
    assert meta['converged'] and meta['valid'] and set(idx) == set(range(40))


def test_core_not_converged_is_invalid():
    rng = np.random.default_rng(2)
    A = rng.normal(size=(60, 3)) * 8
    B = A + 0.0
    B[40:] += 6.0
    _, meta = derive_core(A, B, 1.0, 50, 30, max_trim=1)   # trimming limit reached first
    assert not meta['converged'] and not meta['valid'] and meta['invalid_reason'] == 'not converged'


def test_core_too_small_is_invalid():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(40, 3)) * 8
    B = A + rng.normal(size=(40, 3)) * 3.0          # no rigid subset
    _, meta = derive_core(A, B, 1.0, 50, 30)
    assert not meta['valid']


# --- reference mapping and altlocs ----------------------------------------------------------------

def test_altloc_policy_and_mapping():
    assert select_altloc([('B', 'ALA', (0, 0, 0)), ('A', 'ALA', (1, 1, 1))])[0] == 'A'
    assert select_altloc([('A', 'ALA', (1, 1, 1)), (' ', 'ALA', (2, 2, 2))])[0] == ' '
    assert select_altloc([('B', 'ALA', (1, 1, 1))]) is None
    recs = {1: [(' ', 'MET', (0, 0, 0))], 2: [('B', 'ARG', (1, 0, 0))], 3: [(' ', 'ILE', (2, 0, 0))]}
    xyz, rep = map_reference(recs, 'MRI', [1, 2, 3])
    assert xyz is None and rep['missing'] == [2]
    recs[2] = [(' ', 'LYS', (1, 0, 0))]
    xyz, rep = map_reference(recs, 'MRI', [1, 2, 3])
    assert xyz is None and rep['mismatched'][0]['resseq'] == 2
    recs[2] = [(' ', 'ARG', (1, 0, 0))]
    xyz, rep = map_reference(recs, 'MRI', [1, 2, 3])
    assert rep['valid'] and xyz.shape == (3, 3)


# --- classification and timestamped runs -----------------------------------------------------------

def test_ties_and_boundary():
    lab = classify_frames([1.0, 2.0, 2.0, 3.1], [2.0, 2.0, 1.5, 9.0], radius=2.0)
    assert list(lab) == ['open', 'ambiguous', 'closed', 'ambiguous']


def test_sustained_uses_elapsed_time_not_count():
    t = np.arange(26) * 200.0
    lab25 = np.array(['open'] * 25 + ['ambiguous'])
    lab26 = np.array(['open'] * 26)
    s25 = summarize_runs(runs_with_time(lab25, t, 200.0, 1e-3), 5000.0)
    s26 = summarize_runs(runs_with_time(lab26, t, 200.0, 1e-3), 5000.0)
    assert s25['open']['sustained_runs'] == 0 and abs(s25['open']['longest_span_ns'] - 4.8) < 1e-12
    assert s26['open']['sustained_runs'] == 1


def test_cadence_gap_splits_run():
    t = np.r_[np.arange(20) * 200.0, 3800.0 + 400.0 + np.arange(20) * 200.0]
    runs = runs_with_time(np.array(['open'] * 40), t, 200.0, 1e-3)
    assert len(runs) == 2


def test_start_window():
    t = np.arange(40) * 200.0
    lab = np.array(['open'] * 26 + ['closed'] * 14)
    assert starts_near(lab, t, 'open', 5000.0)
    lab[3] = 'ambiguous'
    assert not starts_near(lab, t, 'open', 5000.0)


# --- outcome precedence ----------------------------------------------------------------------

def traj(start, o_sus, c_sus, amb_frames, n=100):
    s = {'open': {'sustained_runs': o_sus, 'frames': 0}, 'closed': {'sustained_runs': c_sus, 'frames': 0},
         'ambiguous': {'sustained_runs': 0, 'frames': amb_frames}}
    return {'summary': s, 'n_frames': n, 'start_ref': start}


def test_precedence_ambiguity_over_suitable():
    trajs = {'open_start': traj('open', 1, 1, 60), 'closed_start': traj('closed', 0, 1, 0)}
    assert variant_category(trajs, 0.5) == 'unresolved'


def test_categories():
    assert variant_category({'open_start': traj('open', 1, 1, 0),
                             'closed_start': traj('closed', 0, 1, 0)}, 0.5) == 'suitable'
    assert variant_category({'open_start': traj('open', 1, 0, 0),
                             'closed_start': traj('closed', 0, 1, 0)}, 0.5) == 'restricted'
    assert variant_category({'open_start': traj('open', 0, 0, 10),
                             'closed_start': traj('closed', 0, 0, 10)}, 0.5) == 'unsuitable'


def test_variant_disagreement_is_unresolved():
    out, agree = final_outcome({'v1': 'suitable', 'v2': 'restricted'})
    assert out == 'unresolved' and not agree
    assert final_outcome({'v1': 'restricted', 'v2': 'restricted'}) == ('restricted', True)

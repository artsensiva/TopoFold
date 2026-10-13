"""SYNTHETIC tests for the held-out runner, revised reader and revised dependence diagnostics (Amendment 3).
Fixtures are constructed backbones wrapped into a periodic box (not proteins, not MD).

Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_heldout_runner.py
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import heldout_runner as HR  # noqa: E402
from heldout_dependence_r1 import IdentityViolation, classify_series  # noqa: E402
from mdcath_reader_r1 import MdcathDomainR1, make_whole_r1  # noqa: E402
from synthetic_backbone import backbone, trajectory, write_domain  # noqa: E402
from wpd_lib import dihedral  # noqa: E402

L, F, BOX_NM = 55, 205, (3.0, 3.1, 3.2)


@pytest.fixture(scope='module')
def trajs():
    out, base = {}, None
    for i, (t, r) in enumerate([(t, r) for t in ('320', '348') for r in ('0', '1')]):
        X, base = trajectory(L, F, seed=10 + i, base=base, shift=(7.0, -3.0, 11.0))
        out[(t, r)] = X
    return out


def test_generator_reproduces_requested_dihedrals():
    r = np.random.default_rng(0)
    phi, psi = r.uniform(-3, 3, 12), r.uniform(-3, 3, 12)
    X = backbone(phi, psi)
    N, A, C = X[:, 0], X[:, 1], X[:, 2]
    assert np.allclose(np.angle(np.exp(1j * (dihedral(C[:-1], N[1:], A[1:], C[1:]) - phi[1:]))), 0, atol=1e-9)
    assert np.allclose(np.angle(np.exp(1j * (dihedral(N[:-1], A[:-1], C[:-1], N[1:]) - psi[:-1]))), 0, atol=1e-9)


def test_boundary_crossing_fixture_whole_and_rust_numpy_identical_inputs(tmp_path, trajs):
    p = tmp_path / 'wrap.h5'
    write_domain(p, trajs, BOX_NM)
    d = MdcathDomainR1(p)
    bb, O, op = d.backbone_r1('320', '0')
    raw_bond = np.linalg.norm(bb[:, 1:, 0] - bb[:, :-1, 2], axis=-1)
    assert raw_bond.max() > 10.0                                   # the stored molecule really crosses the box
    whole, Ow = make_whole_r1(bb, O, d.box_A('320', '0'))
    ref = trajs[('320', '0')]
    shift = whole[:, :1, 0] - ref[:, :1, 0]                        # whole chain differs by one lattice vector
    assert np.abs(whole - ref[:, :, :3] - shift[:, :, None]).max() < 1e-4
    assert np.abs(Ow - ref[:, :, 3] - shift).max() < 1e-4
    g = HR.trajectory_geometry(d, '320', '0')
    assert max(g['reference_max']) <= 1e-9
    assert np.array_equal(g['CA'], g['CA'].astype(np.float32).astype(np.float64))   # identical rounded inputs
    assert g['burn_in_frames'] == int(np.floor(0.1 * F))


def test_reference_check_failure_is_mandatory(tmp_path, trajs, monkeypatch):
    import topofold as tf
    p = tmp_path / 'ref.h5'
    write_domain(p, trajs, BOX_NM)
    orig = tf.compute_invariants

    def perturbed(x):
        k, t, w = orig(x)
        return k, np.asarray(t) + 1e-6, w
    monkeypatch.setattr(tf, 'compute_invariants', perturbed)
    with pytest.raises(HR.MandatoryCheckError):
        HR.trajectory_geometry(MdcathDomainR1(p), '320', '0')


def test_exit_codes(monkeypatch):
    def boom():
        raise HR.MandatoryCheckError('synthetic')
    monkeypatch.setattr(HR, 'cmd_eligibility', boom)
    assert HR.main(['eligibility']) == 2
    monkeypatch.setattr(HR, 'AUTH', Path('/nonexistent/authorization.json'))
    assert HR.main(['score-heldout']) == 3


def test_eligibility_ledger_cases(tmp_path, trajs):
    good = tmp_path / 'good.h5'
    write_domain(good, trajs, BOX_NM, missing_o=(20,))
    rec = HR.structural_eligibility(good)
    assert rec['eligible'] and rec['descriptive']['residues_without_unique_O'] == [20]
    bad = tmp_path / 'nan.h5'
    write_domain(bad, trajs, BOX_NM, nonfinite=(('348', '1'), 7))
    rec = HR.structural_eligibility(bad)
    assert not rec['eligible'] and any('non-finite' in x for x in rec['reasons'])
    short = tmp_path / 'short.h5'
    write_domain(short, {k: v[:150] for k, v in trajs.items()}, BOX_NM)
    assert any('numFrames 150' in x for x in HR.structural_eligibility(short)['reasons'])
    with pytest.raises(HR.MandatoryCheckError):
        HR.structural_eligibility(good, expected_size=1, expected_sha='0' * 64)


@pytest.fixture(scope='module')
def scored(tmp_path_factory, trajs):
    tmp = tmp_path_factory.mktemp('score')
    a, b = tmp / 'withO.h5', tmp / 'missingO.h5'
    write_domain(a, trajs, BOX_NM)
    write_domain(b, trajs, BOX_NM, missing_o=(20,))
    return HR.score_domain(a, 1), HR.score_domain(b, 1)


def test_optional_oxygen_affects_only_peptide_plane_windows(scored):
    full, miss = scored
    for key in full['cells']:
        if key.endswith('phipsi'):
            assert full['cells'][key]['domain_nr10'] == miss['cells'][key]['domain_nr10']
    for W in (5, 8):
        und = [w['s'] for w in miss['cells'][f'320/W{W}/pp']['windows_undefined']]
        expect = [s for s in range(L - W + 1) if s <= 19 <= s + W - 2]   # residue 20 = index 19 in r = s…s+W−2
        assert und == expect
        assert not full['cells'][f'320/W{W}/pp']['windows_undefined']


def test_domain_stage_outputs(scored):
    full, _ = scored
    c = full['cells']['348/W8/phipsi']
    assert set(c['domain_nr10']) == set(HR.METHODS) and all(np.isfinite(v) for v in c['domain_nr10'].values())
    assert c['n_matched_windows'] >= HR.MIN_WINDOWS or c['domain_nr10_c1_matched'] is None
    assert set(c['dependence']) == {f'{x}-{y}' for i, x in enumerate(HR.CANDIDATES) for y in HR.CANDIDATES[i + 1:]}
    assert full['geometry']['320']['query_rep0']['reference_max_tau_kappa'][0] <= 1e-9


def test_cell_stage_minimum_domains_and_candidate_A():
    r = np.random.default_rng(3)

    def fake(n):
        doms = []
        for _ in range(n):
            vals = {m: float(r.uniform(0.6, 0.9)) for m in HR.METHODS}
            vals['A'] = float(r.uniform(0.4, 0.5))
            cell = {'domain_nr10': vals, 'domain_nr10_c1_matched': dict(vals)}
            doms.append({'cells': {f'{t}/W{W}/phipsi': cell for t in HR.TEMPS for W in HR.WINDOWS}})
        return doms
    small = HR.analyse_cells(fake(20), 'phipsi', n_boot=50)
    assert all(v['status'] == 'insufficient_evidence' for v in small['cells'].values())
    big = HR.analyse_cells(fake(30), 'phipsi', n_boot=100)
    assert all(c['observed']['best'] == 'A' for c in big['cells'].values())
    assert big['overall']['A'] in ('meets_criterion_in_all_cells', 'inconclusive')


# --- dependence diagnostics r1 --------------------------------------------------------------------------------

def test_no_time_compression_and_unavailable_vs_violation():
    r = np.random.default_rng(1)
    s = r.normal(size=400)
    s[150:160] = np.nan
    c = classify_series(s)
    assert c['segment'] == [160, 400] and c['n_gaps'] == 1 and c['n_unavailable'] == 10
    z = np.zeros(300)
    z[5] = np.nan
    ci = classify_series(z, identity=True)
    assert ci['class'] == 'structural_identity' and ci['n_unavailable'] == 1
    z[7] = 1e-6
    with pytest.raises(IdentityViolation):
        classify_series(z, identity=True)
    assert classify_series(np.r_[r.normal(size=40), np.nan, r.normal(size=40)])['class'] == 'insufficient_observations'
    assert classify_series(np.full(300, 0.02))['class'] == 'constant_observed'
    white = classify_series(r.normal(size=5000))
    assert white['class'] == 'estimated_below_heuristic_cutoff'

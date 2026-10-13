"""SYNTHETIC tests of the actual scoring command path (execution corrections EC2–EC4).
Constructed backbones in mdCATH layout; no mdCATH data. Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_heldout_exec.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import heldout_runner as HR  # noqa: E402
from synthetic_backbone import trajectory, write_domain  # noqa: E402

L, F, BOX = 50, 200, (3.0, 3.1, 3.2)


def make_cohort(tmp, n=2):
    data = tmp / 'data'
    data.mkdir()
    entries = []
    for k in range(n):
        name = f'syn{k:02d}A00'
        trajs, base = {}, None
        for i, (t, r) in enumerate([(t, r) for t in ('320', '348') for r in ('0', '1')]):
            trajs[(t, r)], base = trajectory(L, F, seed=100 * k + i, base=base, shift=(5.0, 1.0, -2.0))
        p = data / f'mdcath_dataset_{name}.h5'
        write_domain(p, trajs, BOX, name=name)
        entries.append({'order': k + 1, 'domain': name, 'role': 'synthetic', 'sha256': HR.sha256(p),
                        'size_bytes': p.stat().st_size})
    sel = tmp / 'selection.json'
    sel.write_text(json.dumps({'synthetic': True}))
    cohort = tmp / 'cohort.json'
    cohort.write_text(json.dumps({'_meta': {'synthetic': True}, 'domains': entries}))
    return data, sel, cohort


def authorize(tmp, sel, cohort, **override):
    a = {'authorized_by': 'SYNTHETIC test', 'scope': 'test', **HR.expected_authorization(sel, cohort), **override}
    p = tmp / 'auth.json'
    p.write_text(json.dumps(a))
    return p


@pytest.fixture(scope='module')
def env(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('exec')
    data, sel, cohort = make_cohort(tmp)
    return tmp, data, sel, cohort


def args(tmp, data, sel, cohort, auth, out='out.json', ck='ck'):
    return ['--cohort', str(cohort), '--authorization', str(auth), '--selection', str(sel), '--data', str(data),
            '--out', str(tmp / out), '--checkpoints', str(tmp / ck), '--workers', '1']


def test_authorization_absent_malformed_mismatched(env):
    tmp, data, sel, cohort = env
    assert HR.cmd_score(args(tmp, data, sel, cohort, tmp / 'none.json')) == 3
    bad = tmp / 'bad.json'
    bad.write_text('{not json')
    assert HR.cmd_score(args(tmp, data, sel, cohort, bad)) == 3
    missing = tmp / 'missing.json'
    missing.write_text(json.dumps({'authorized_by': 'x'}))
    assert HR.cmd_score(args(tmp, data, sel, cohort, missing)) == 3
    wrong = authorize(tmp, sel, cohort, cohort_sha256='0' * 64)
    assert HR.cmd_score(args(tmp, data, sel, cohort, wrong)) == 3
    impl = HR.implementation_hashes()
    impl[next(iter(impl))] = 'f' * 64
    wrong2 = authorize(tmp, sel, cohort, implementation_sha256=impl)
    assert HR.cmd_score(args(tmp, data, sel, cohort, wrong2)) == 3


def test_integrity_and_numerical_failures_exit_2(env, monkeypatch, tmp_path):
    tmp, data, sel, cohort = env
    c = json.loads(cohort.read_text())
    c['domains'][0]['sha256'] = '0' * 64
    bad_cohort = tmp_path / 'cohort_bad.json'
    bad_cohort.write_text(json.dumps(c))
    auth = authorize(tmp_path, sel, bad_cohort)
    assert HR.cmd_score(args(tmp_path, data, sel, bad_cohort, auth)) == 2
    assert (tmp_path / 'out.FAILED.json').exists()
    import topofold as tf
    orig = tf.compute_invariants
    monkeypatch.setattr(tf, 'compute_invariants', lambda x: (lambda k, t, w: (k, np.asarray(t) + 1e-6, w))(*orig(x)))
    auth = authorize(tmp_path, sel, cohort)
    assert HR.cmd_score(args(tmp_path, data, sel, cohort, auth, out='o2.json', ck='ck2')) == 2
    assert 'reference' in json.loads((tmp_path / 'o2.FAILED.json').read_text())['mandatory_failure']


def test_horn_kabsch_failure_is_mandatory(env, monkeypatch, tmp_path):
    tmp, data, sel, cohort = env
    monkeypatch.setattr(HR, 'horn_rmsd', lambda a, b: HR.__dict__['d_rms'](a.reshape(len(a), -1),
                                                                         b.reshape(len(b), -1)) * 1.01)
    with pytest.raises(HR.MandatoryCheckError, match='Horn vs Kabsch'):
        HR.score_domain(data / 'mdcath_dataset_syn00A00.h5', 1)


@pytest.fixture(scope='module')
def full_run(env):
    tmp, data, sel, cohort = env
    mp = pytest.MonkeyPatch()
    mp.setattr(HR.HD2._mod, 'MIN_DOMAINS', 2)
    auth = authorize(tmp, sel, cohort)
    out = HR.run_cohort(cohort, auth, sel, data, tmp / 'full.json', tmp / 'ck_full', workers=1, n_boot=200)
    yield out, auth
    mp.undo()


def test_final_aggregation_both_targets(full_run):
    out, _ = full_run
    for tgt in ('phipsi', 'pp'):
        agg = out['aggregate'][tgt]
        assert set(agg['overall']) == set(HR.METHODS)
        c = agg['cells']['320/W8/' + tgt]
        assert c['n_defined_domains'] == 2 and 'A' in c['methods']
    assert out['secondary']['348/W5/phipsi']['n_defined_domains'] == 2
    assert all(d['horn_kabsch_max_A'] <= 1e-6 for d in out['domains'])


def test_resumed_equals_uninterrupted_and_mismatched_checkpoint_recomputed(env, full_run, monkeypatch):
    tmp, data, sel, cohort = env
    full, auth = full_run
    monkeypatch.setattr(HR.HD2._mod, 'MIN_DOMAINS', 2)
    real = HR.score_domain

    def stop_on_second(path, idx):
        if idx == 2:
            raise HR.MandatoryCheckError('SYNTHETIC interruption')
        return real(path, idx)
    monkeypatch.setattr(HR, 'score_domain', stop_on_second)
    with pytest.raises(HR.MandatoryCheckError):
        HR.run_cohort(cohort, auth, sel, data, tmp / 'res.json', tmp / 'ck_res', workers=1, n_boot=200)
    assert (tmp / 'ck_res' / 'syn00A00.json').exists() and not (tmp / 'ck_res' / 'syn01A00.json').exists()
    monkeypatch.setattr(HR, 'score_domain', real)
    res = HR.run_cohort(cohort, auth, sel, data, tmp / 'res.json', tmp / 'ck_res', workers=1, n_boot=200)
    assert res['_meta']['reused_checkpoints'] == ['syn00A00']
    assert res['aggregate'] == full['aggregate'] and res['domains'] == full['domains']
    ck = tmp / 'ck_res' / 'syn01A00.json'
    d = json.loads(ck.read_text())
    d['key']['input_sha256'] = '1' * 64
    ck.write_text(json.dumps(d))
    res2 = HR.run_cohort(cohort, auth, sel, data, tmp / 'res2.json', tmp / 'ck_res', workers=1, n_boot=200)
    assert res2['_meta']['reused_checkpoints'] == ['syn00A00'] and res2['aggregate'] == full['aggregate']

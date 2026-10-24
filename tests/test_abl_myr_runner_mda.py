"""SYNTHETIC tests of the ABL-MYR-1 worker functions that need MDAnalysis (topology parsing, residue mapping, XTC time
stamps and atom counts, φ/ψ extraction, QCP-based RMSF). Engineering checks only; no Dryad file is used."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1"))

# With TOPOFOLD_REQUIRE_MDANALYSIS_TESTS=1 (the CI adapter job) a missing MDAnalysis fails instead of skipping.
if os.environ.get("TOPOFOLD_REQUIRE_MDANALYSIS_TESTS") == "1":
    import MDAnalysis  # noqa: F401
else:
    pytest.importorskip("MDAnalysis")
import run_abl_myr as R  # noqa: E402
from test_abl_myr_runner import (AUTHOR_SHIFT, FIRST, LAST, N_FRAMES, base_backbone, frames_for, pdb_text,  # noqa: E402
                                 residues, write_xtc)

ANCH = R.CFG["anchors"]
FASTA = str(ROOT / R.CFG["reference"]["fasta"])
N = LAST - FIRST + 1


def topo(pdb, anchors=ANCH):
    return R._w_topo({"topology": str(pdb), "fasta": FASTA, "anchors": anchors, "accession": "P00519",
                      "isoform": "P00519-2"})


def make_pdb(tmp_path, drop=()):
    res, X0 = residues(), base_backbone()
    keep = [i for i in range(N) if i not in drop]
    p = tmp_path / "ref.pdb"
    p.write_text(pdb_text([res[i] for i in keep], X0[keep]))
    return p, keep


def test_topology_facts_and_mapping_use_the_declared_isoform_numbering(tmp_path):
    pdb, _ = make_pdb(tmp_path)
    t = topo(pdb)
    assert t["chains"] == ["A"] and t["n_atoms"] == 4 * N and set(t["names"]) == {"N", "CA", "C", "O"}
    assert [r[3] for r in t["residues"]] == [None] + [True] * (N - 1)               # bonded flags from geometry
    assert all(r[4:7] == [True, True, True] for r in t["residues"])
    m = t["mapping"]
    assert m["accepted"] and m["report"]["roles"]["DFG_Asp"]["reference_index"] == 400
    assert m["report"]["roles"]["DFG_Asp"]["resid"] == 400 - AUTHOR_SHIFT           # author numbering is 1b − 19
    assert m["report"]["alignment"]["author_numbering_diagnostic"]["reference_index_minus_author_resid"] == {
        str(AUTHOR_SHIFT): N}
    R.check_g2(t, {"x": {"n_atoms": 4 * N}})                                        # atom identity topology vs trajectory
    with pytest.raises(R.GateError, match="topology has"):
        R.check_g2(t, {"x": {"n_atoms": 4 * N - 4}})


def test_mapping_refuses_wrong_declared_numbering_and_reports_gaps(tmp_path):
    pdb, _ = make_pdb(tmp_path)
    shifted = {k: v - 1 for k, v in ANCH.items()}                                   # literature-style offset, not applied
    m = topo(pdb, shifted)["mapping"]
    assert not m["accepted"] and any("expected" in p for p in m["problems"])
    assert m["report"]["diagnostic_offsets_not_applied"]["d"] == [1]
    drop = set(range(402 - FIRST, 408 - FIRST))                                     # six residues after DFG-Gly are absent
    pdb2, keep = make_pdb(tmp_path / "gap" if (tmp_path / "gap").mkdir() is None else tmp_path, drop=drop)
    t = topo(pdb2)
    assert t["mapping"]["accepted"] and t["mapping"]["report"]["roles"]["DFGxxX"]["status"] == "unobserved"
    assert [i for i, r in enumerate(t["residues"]) if r[3] is False] == [keep.index(408 - FIRST)]
    assert t["mapping"]["report"]["alignment"]["unobserved_reference_ranges"] == [[402, 407]]


def test_xtc_facts_report_atom_count_frames_and_time_stamps(tmp_path):
    pdb, _ = make_pdb(tmp_path)
    X = frames_for(np.random.default_rng(1), np.full(N, 0.1), 12)
    write_xtc(tmp_path / "a.xtc", pdb, X, t0_ps=0.0)
    write_xtc(tmp_path / "b.xtc", pdb, X[:10], t0_ps=1.0e6, times=[1.0e6 + 2000.0 * i + (50.0 if i > 5 else 0.0)
                                                                    for i in range(10)])
    f = R._w_xtcs({"files": {"a": str(tmp_path / "a.xtc"), "b": str(tmp_path / "b.xtc")}})
    assert f["a"]["n_atoms"] == 4 * N and f["a"]["n_frames"] == 12 and f["a"]["times_ps"][:3] == [0.0, 2000.0, 4000.0]
    assert R.frame_gate(f["a"]["times_ps"], {**R.CFG, "min_frames": 1})["n_analysis"] == 0
    with pytest.raises(R.GateError, match="interval"):                               # 2050 ps step inside the analysis frames
        R.frame_gate(f["b"]["times_ps"], R.CFG)
    assert R._w_count({"topology": str(pdb)}) == {"n_atoms": 4 * N}


def numpy_phi_psi(X, bonded):
    def dih(p0, p1, p2, p3):
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1)
        v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
        return np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))
    n = len(X)
    phi, psi = np.full(n, np.nan), np.full(n, np.nan)
    for i in range(1, n):
        if bonded[i]:
            phi[i] = dih(X[i - 1, 2], X[i, 0], X[i, 1], X[i, 2])
    for i in range(n - 1):
        if bonded[i + 1]:
            psi[i] = dih(X[i, 0], X[i, 1], X[i, 2], X[i + 1, 0])
    return phi, psi


def test_phi_psi_extraction_matches_an_independent_formula_and_flags_invalid_residues(tmp_path):
    drop = {60, 61}
    pdb, keep = make_pdb(tmp_path, drop=drop)
    t = topo(pdb)
    X = frames_for(np.random.default_rng(3), np.full(N, 0.05), 6)[:, keep]
    X[2:, 30] += 15.0                                                               # SYNTHETIC displaced residue, frames 2+
    write_xtc(tmp_path / "r.xtc", pdb, X, t0_ps=0.0)
    R._w_phipsi({"topology": str(pdb), "chain": t["chain"], "residues": t["residues"], "bond_max": 2.0,
                 "runs": {"r": {"xtc": str(tmp_path / "r.xtc"), "start": 1}}, "npz": str(tmp_path / "pp.npz")})
    pp = np.load(tmp_path / "pp.npz")
    phi, psi = pp["r__phi"], pp["r__psi"]
    assert phi.shape == (5, len(keep))                                              # frames 1..5 (start = 1)
    bonded = [bool(r[3]) if r[3] is not None else False for r in t["residues"]]
    for f in (0, 1):                                                                # frames 1 and 2 of the file: 1 is intact
        ref_phi, ref_psi = numpy_phi_psi(X[f + 1], bonded)
        cols = [i for i in range(len(keep)) if i not in (29, 30, 31)]
        ok = np.isfinite(ref_phi[cols])
        if f == 0:
            for got, want in ((phi[f], ref_phi), (psi[f], ref_psi)):
                m = [i for i in cols if np.isfinite(want[i])]
                assert np.isfinite(got[m]).all() and len(m) > 100
                assert np.abs(np.angle(np.exp(1j * (got[m] - want[m])))).max() < 0.03      # circular difference
                assert np.array_equal(np.isfinite(got[cols]), np.isfinite(want[cols]))
    gap = keep.index(62) if 62 in keep else None
    first_after_gap = [i for i, r in enumerate(t["residues"]) if r[3] is False][0]
    assert np.isnan(phi[:, 0]).all() and np.isnan(psi[:, -1]).all()                 # termini have no φ / ψ
    assert np.isnan(phi[:, first_after_gap]).all() and np.isnan(psi[:, first_after_gap - 1]).all()   # chain break
    assert gap is not None
    assert np.isfinite(phi[0, 30]) and np.isnan(phi[1:, 30]).all() and np.isnan(psi[1:, 29]).all()    # displaced residue
    assert np.isnan(phi[1:, 31]).all() and np.isfinite(psi[0, 30]) and np.isnan(psi[1:, 30]).all()


def kabsch_rmsf(X):
    ref = X[0] - X[0].mean(0)
    al = []
    for f in X:
        m = f - f.mean(0)
        u, _, vt = np.linalg.svd(m.T @ ref)
        d = np.sign(np.linalg.det(u @ vt))
        r = u @ np.diag([1, 1, d]) @ vt
        al.append(m @ r)
    al = np.array(al)
    return np.sqrt(((al - al.mean(0)) ** 2).sum(-1).mean(0))


def random_rotation(rng):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_rmsf_uses_the_qcp_rotation_convention_and_is_invariant_to_rigid_motion(tmp_path):
    rng = np.random.default_rng(5)
    ca = frames_for(rng, np.full(N, 0.4), 30)[:, :, 1]
    moved = np.array([f @ random_rotation(rng).T + rng.normal(scale=20, size=3) for f in ca])
    np.savez(tmp_path / "in.npz", plain=ca, moved=moved)
    R._w_rmsf({"npz_in": str(tmp_path / "in.npz"), "npz_out": str(tmp_path / "out.npz")})
    out = np.load(tmp_path / "out.npz")
    np.testing.assert_allclose(out["plain"], kabsch_rmsf(ca), atol=1e-6)            # independent SVD superposition
    np.testing.assert_allclose(out["moved"], out["plain"], atol=1e-6)               # rigid motion changes nothing
    assert out["plain"].min() > 0.2 and out["plain"].shape == (N,)
    bad = ca.copy()
    bad[:, 7] = np.nan
    np.savez(tmp_path / "nan.npz", bad=bad)
    R._w_rmsf({"npz_in": str(tmp_path / "nan.npz"), "npz_out": str(tmp_path / "nan_out.npz")})
    r = np.load(tmp_path / "nan_out.npz")["bad"]
    assert np.isnan(r[7]) and np.isfinite(np.delete(r, 7)).all()                    # non-finite rows are excluded, not zeroed


def test_worker_dispatch_writes_a_json_result(tmp_path):
    pdb, _ = make_pdb(tmp_path)
    pj = tmp_path / "p.json"
    import json
    pj.write_text(json.dumps({"topology": str(pdb), "out": str(tmp_path / "o.json")}))
    assert R.main(["_worker", "count", str(pj)]) == 0
    assert json.loads((tmp_path / "o.json").read_text()) == {"n_atoms": 4 * N}
    assert N_FRAMES == 40

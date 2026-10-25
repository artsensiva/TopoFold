"""Tests of the KRAS K1 crystal endpoint (step S1) on SYNTHETIC structures: analytic references, invariance, ties, bounds."""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
K1 = ROOT / "benchmarks/audit/kras_k1"
sys.path.insert(0, str(K1))
import k1_endpoint as K  # noqa: E402

P = json.loads((K1 / "protocol_params.json").read_text())
EP = P["endpoint"]
AA = "ACDEFGHIKLMNPQRSTVWY"


def _rot(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    return q if np.linalg.det(q) > 0 else q * np.array([1, 1, -1])


def test_kabsch_recovers_a_known_rigid_motion_and_is_a_proper_rotation():
    rng = np.random.default_rng(1)
    ref = rng.normal(size=(40, 3)) * 8
    r_true, t_true = _rot(rng), np.array([5.0, -3.0, 2.0])
    mov = (ref - t_true) @ r_true            # ref = mov @ r_true.T + t_true
    r, t = K.kabsch(mov, ref)
    assert np.linalg.det(r) == pytest.approx(1.0)
    assert np.abs(ref - (mov @ r.T + t)).max() < 1e-9


def test_a_mirror_image_cannot_be_fitted_by_a_proper_rotation():
    rng = np.random.default_rng(2)
    ref = rng.normal(size=(40, 3)) * 8
    mirrored = ref * np.array([-1.0, 1.0, 1.0])
    r, t = K.kabsch(mirrored, ref)
    assert np.linalg.det(r) == pytest.approx(1.0)
    assert np.linalg.norm(ref - (mirrored @ r.T + t), axis=1).mean() > 1.0


def test_trimmed_fit_returns_the_injected_local_displacement_exactly():
    rng = np.random.default_rng(3)
    n = 50
    ref = rng.normal(size=(n, 3)) * 8
    mov = ref.copy()
    mov[10:15] += np.array([3.0, 0.0, 0.0])           # 5 residues displaced by exactly 3 Å
    mov = mov @ _rot(rng).T + np.array([1.0, 2.0, 3.0])
    d = K.pair_displacement(ref, mov, list(range(n)), 0.2)   # trims floor(0.2 * 50) = 10 residues
    assert d[10:15] == pytest.approx(3.0, abs=1e-9)
    assert np.delete(d, range(10, 15)).max() < 1e-9


def test_pair_displacement_is_invariant_to_rigid_motion_of_either_chain_and_symmetric():
    rng = np.random.default_rng(4)
    a = rng.normal(size=(30, 3)) * 8
    b = a + rng.normal(size=(30, 3)) * 0.3
    base = K.pair_displacement(a, b, list(range(30)), 0.2)
    moved = K.pair_displacement(a @ _rot(rng).T + 7.0, b @ _rot(rng).T - 4.0, list(range(30)), 0.2)
    swapped = K.pair_displacement(b, a, list(range(30)), 0.2)
    assert moved == pytest.approx(base, abs=1e-9)
    assert swapped == pytest.approx(base, abs=1e-9)


def test_trim_ties_drop_the_lower_residue_number_first():
    disp = np.array([1.0, 1.0, 1.0, 0.5])
    assert K.trim_indices(disp, [30, 10, 20, 40], 2) == [1, 2]         # residues 10 and 20 dropped; 30 kept
    assert K.trim_indices(np.array([2.0, 3.0, 1.0]), [1, 2, 3], 1) == [1]


def test_positive_rule_is_inclusive_within_tolerance_and_takes_the_larger_threshold():
    b = np.array([1.0, 1.0 - 1e-6, 2.0, 1.9999999995, 0.99])
    w = np.array([0.1, 0.1, 1.0, 1.0, 0.4])           # thresholds: 1.0, 1.0, 2.0, 2.0, 1.0
    thr, pos, bnd = K.compute_positive(b, w, 1.0, 2.0)
    assert thr.tolist() == [1.0, 1.0, 2.0, 2.0, 1.0]
    assert pos.tolist() == [True, False, True, True, False]
    assert bnd.tolist() == [True, False, True, True, False]


def test_segments_merge_rule_and_length_counts_positives():
    assert K.segments([10, 11, 13], 2) == [[10, 11, 13]]              # gap of 1 residue
    assert K.segments([10, 13], 2) == [[10, 13]]                      # gap of 2 residues: merged
    assert K.segments([10, 14], 2) == [[10], [14]]                    # gap of 3 residues: separate
    assert K.segments([58, 59, 68, 69], 2) == [[58, 59], [68, 69]]    # the 60-67 gap separates segments
    assert [len(s) for s in K.segments([10, 13, 40], 2)] == [2, 1]


G = EP["gate"]


def _gate(res, n=153):
    return K.gate(res, n, G, EP["segments"]["merge_gap"])


def test_gate_bounds_are_inclusive():
    eight = [10, 11, 12, 13] + [50, 51, 52, 53]
    assert _gate(eight)["passes"] and not _gate(eight[:-1])["passes"]
    assert _gate(list(range(10, 29)) + list(range(60, 79)))["n_positive"] == 38 and _gate(list(range(10, 29)) + list(range(60, 79)))["passes"]
    assert not _gate(list(range(10, 29)) + list(range(60, 80)))["passes"]          # 39 > floor(0.25 * 153) = 38
    assert _gate(eight)["max_positive"] == math.floor(0.25 * 153) == 38


def test_gate_needs_two_segments_of_at_least_three_residues():
    assert _gate([10, 11, 12, 13, 14, 15, 50, 51, 52])["passes"]               # 6 and 3
    assert not _gate([10, 11, 12, 13, 14, 15, 16, 50, 51])["passes"]           # 7 and 2: only one long segment
    assert not _gate(list(range(10, 20)))["passes"]                            # one segment of 10


def _pdb(chains: dict) -> str:
    lines, serial = [], 1
    for chain, residues in chains.items():
        for resid, (resname, xyz) in sorted(residues.items()):
            lines.append(f"ATOM  {serial:5d} {' CA ':4s} {resname:>3s} {chain}{resid:4d}    "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{20.0:6.2f}           C")
            serial += 1
    lines.append("HETATM    1  CA  CA  Z 900       0.000   0.000   0.000  1.00  0.00          CA")
    return "\n".join(lines) + "\n"


def _synthetic(seed=7, displaced=(range(30, 39), range(100, 111)), shift=3.0, noise=0.05):
    """SYNTHETIC 8-chain G12D-like structures: a random-walk Calpha trace, GNP chains displaced by `shift` Å in `displaced`."""
    rng = np.random.default_rng(seed)
    topo = "".join(rng.choice(list(AA), size=161))                      # residues 4-164
    names = {r: next(k for k, v in K.ONE.items() if v == topo[r - 4]) for r in range(4, 165)}
    base = np.cumsum(rng.normal(size=(166, 3)) * 2.2, axis=0)
    missing = {("GDP", "A"): range(60, 65), ("GDP", "B"): range(60, 67), ("GNP", "A"): range(61, 64), ("GNP", "D"): range(61, 68)}
    out = {}
    for state, chains in (("GDP", "AB"), ("GNP", "ABCDEF")):
        text = {}
        for c in chains:
            xyz = base.copy()
            if state == "GNP":
                for rg in displaced:
                    for r in rg:
                        xyz[r - 1] += np.array([shift, 0.0, 0.0])
            xyz = xyz + rng.normal(size=xyz.shape) * noise
            xyz = xyz @ _rot(rng).T + rng.normal(size=3) * 10
            gone = set(missing.get((state, c), ()))
            text[c] = {r: ("ASP" if r == 12 else names.get(r, "GLY"), xyz[r - 1]) for r in range(1, 167) if r not in gone}
        out[state] = _pdb(text)
    return out["GDP"], out["GNP"], topo


def test_endpoint_pipeline_recovers_the_injected_segments_exactly():
    gdp, gnp, topo = _synthetic()
    res = K.analyse(gdp, gnp, P, topo)
    assert len(res["scored"]) == 153 and 60 not in res["scored"] and 67 not in res["scored"]
    assert res["E"] == list(range(30, 39)) + list(range(100, 111))
    assert res["gate"]["passes"] and [s["length"] for s in res["gate"]["segments"]] == [9, 11]
    assert res["n_between_pairs"] == 12 and res["n_within_pairs"] == 16
    b = res["b_mean"]
    pos = np.array(res["positive"])
    assert b[pos] == pytest.approx(3.0, abs=0.2) and b[~pos].max() < 0.5
    assert res["sensitivity_sets"]["primary_pair_only"]["E"] == res["E"]
    assert res["sensitivity_sets"]["median_over_pairs"]["jaccard_with_E"] == 1.0


def test_endpoint_is_deterministic_and_unchanged_by_moving_one_chain():
    gdp, gnp, topo = _synthetic(seed=11)
    r1, r2 = K.analyse(gdp, gnp, P, topo), K.analyse(gdp, gnp, P, topo)
    assert r1["E"] == r2["E"] and np.array_equal(r1["b_mean"], r2["b_mean"])
    ca, names = K.parse_ca(gnp)
    rng = np.random.default_rng(5)
    rot, shift = _rot(rng), np.array([30.0, -20.0, 10.0])
    moved = {c: {r: (names[c][r], v @ rot.T + shift) for r, v in ca[c].items()} if c == "C" else {r: (names[c][r], v) for r, v in ca[c].items()} for c in ca}
    r3 = K.analyse(gdp, _pdb(moved), P, topo)
    assert r3["E"] == r1["E"] and r3["b_mean"] == pytest.approx(r1["b_mean"], abs=2e-3)   # PDB text rounds to 1e-3 Å


def test_a_failed_gate_is_reported_not_repaired():
    gdp, gnp, topo = _synthetic(displaced=(range(30, 39),))
    res = K.analyse(gdp, gnp, P, topo)
    assert res["E"] == list(range(30, 39)) and not res["gate"]["passes"]
    assert res["gate"]["checks"]["segments_of_min_length_at_least_min"] is False
    gdp2, gnp2, topo2 = _synthetic(displaced=(), shift=0.0)
    assert not K.analyse(gdp2, gnp2, P, topo2)["gate"]["passes"]            # nothing differs: |E| = 0 < 8


def test_identity_checks_stop_on_any_inconsistency():
    gdp, gnp, topo = _synthetic()
    with pytest.raises(ValueError, match="not ASP"):
        K.analyse(gdp.replace(" ASP A  12", " GLY A  12"), gnp, P, topo)
    with pytest.raises(ValueError, match="topology sequence"):
        other = "".join("A" if c != "A" else "C" for c in topo)
        K.analyse(gdp, gnp, P, other)
    with pytest.raises(ValueError, match="chains found"):
        K.analyse(gdp, "\n".join(ln for ln in gnp.splitlines() if " F " not in ln[16:23]) + "\n", P, topo)
    with pytest.raises(ValueError, match="frozen KRAS-STRUCT-2 set"):
        K.analyse(gdp, "\n".join(ln for ln in gnp.splitlines() if not (ln.startswith("ATOM") and ln[21] in "ABCDEF" and int(ln[22:26]) == 20)) + "\n", P, topo)


def test_parse_ca_rejects_altlocs_insertion_codes_duplicates_and_zero_occupancy_and_ignores_hetatm():
    good = _pdb({"A": {5: ("ALA", (1.0, 2.0, 3.0))}})
    ca, names = K.parse_ca(good)
    assert list(ca) == ["A"] and names["A"][5] == "ALA"                      # the HETATM CA of chain Z is ignored
    line = good.splitlines()[0]
    with pytest.raises(ValueError, match="alternate"):
        K.parse_ca(line[:16] + "A" + line[17:])
    with pytest.raises(ValueError, match="insertion"):
        K.parse_ca(line[:26] + "B" + line[27:])
    with pytest.raises(ValueError, match="duplicate"):
        K.parse_ca(line + "\n" + line)
    with pytest.raises(ValueError, match="occupancy"):
        K.parse_ca(line[:54] + "  0.00" + line[60:])

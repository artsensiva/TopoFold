"""SYNTHETIC tests of the ABL-MYR-1b comparability check (benchmarks/audit/abl_myr_1b/run_comp_check.py).

Engineering checks only; archives and run-input records are generated here and nothing is read from Dryad. Real md.tpr files
cannot be generated without GROMACS, so run_check is exercised with a loader returning records of the same shape as
load_tpr. Needs MDAnalysis (minimum-image utilities)."""

import hashlib
import json
import os
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1b"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1"))
if os.environ.get("TOPOFOLD_REQUIRE_MDANALYSIS_TESTS") == "1":
    import MDAnalysis  # noqa: F401
else:
    pytest.importorskip("MDAnalysis")
import run_comp_check as C  # noqa: E402
from test_abl_myr_runner import base_backbone, residues  # noqa: E402

SPEC = json.loads((ROOT / "benchmarks/audit/abl_myr_1b/check_spec.json").read_text())
BUDGET = {"id": "X", "budget": {"extracted_bytes_max": 10 ** 9, "cpu_seconds_max": 1800}}
BOX_NM = 9.0
MDP = ("integrator = md\ntcoupl = V-rescale\ntc-grps = {grp} Water_NA_CL\ntau_t = 1.0 1.0\nref_t = 300 300\n"
       "pcoupl = Parrinello-Rahman\ngen_vel = no\ncontinuation = yes\ndefine = -DPOSRES\n")


# --- SYNTHETIC system ----------------------------------------------------------------------------------------------------

def system_atoms(myr: bool, n_sol=30, n_na=3, n_cl=2):
    """[(resid, resname, atom name, xyz nm)] for a SYNTHETIC protein (N, CA, C, O) plus ligands, water and ions."""
    out = []
    X = base_backbone() / 10.0 + 1.0
    for (resid, rn), xyz in zip(residues(), X):
        for nm, p in zip(("N", "CA", "C", "O"), xyz):
            out.append((resid, rn, nm, p))
    nxt = 1000
    if myr:
        out += [(nxt, "MYRP", "C1", np.array([1.0, 1.0, 5.0])), (nxt, "MYRP", "C2", np.array([1.1, 1.0, 5.0]))]
    out += [(nxt + 1, "ATP", "PG", np.array([2.0, 2.0, 5.0])), (nxt + 2, "MG", "MG", np.array([2.1, 2.0, 5.0]))]
    for k in range(n_sol):
        o = np.array([0.3 + 0.25 * (k % 10), 7.0, 0.5 + 0.3 * (k // 10)])
        out += [(nxt + 3 + k, "SOL", "OW", o), (nxt + 3 + k, "SOL", "HW1", o + 0.09), (nxt + 3 + k, "SOL", "HW2", o - 0.09)]
    for k in range(n_na):
        out.append((nxt + 100 + k, "NA", "NA", np.array([8.0, 1.0 + k, 1.0])))
    for k in range(n_cl):
        out.append((nxt + 200 + k, "CL", "CL", np.array([8.0, 1.0 + k, 3.0])))
    return out


def gro(atoms, box_line=f"   {BOX_NM:.5f}   {BOX_NM:.5f}   {BOX_NM:.5f}"):
    lines = ["SYNTHETIC", str(len(atoms))]
    for i, (resid, rn, nm, p) in enumerate(atoms, start=1):
        lines.append(f"{resid % 100000:5d}{rn:<5s}{nm:>5s}{i % 100000:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
    return "\n".join(lines) + "\n" + box_line + "\n"


def ndx(atoms, myr):
    solute = [i + 1 for i, a in enumerate(atoms) if a[1] not in ("SOL", "NA", "CL")]
    solvent = [i + 1 for i, a in enumerate(atoms) if a[1] in ("SOL", "NA", "CL")]
    g = "Protein_MYRP_ATP_MG" if myr else "Protein_ATP_MG"
    return f"[ {g} ]\n" + " ".join(map(str, solute)) + "\n[ Water_NA_CL ]\n" + " ".join(map(str, solvent)) + "\n"


def make_record(atoms, rng, *, noise=0.03, vel=True, vel_rng=None):
    """A load_tpr-shaped SYNTHETIC record (Å) for one run input."""
    pos = np.array([a[3] for a in atoms]) * 10.0
    prot = [i for i, a in enumerate(atoms) if a[1] not in ("MYRP", "ATP", "MG", "SOL", "NA", "CL")]
    pos[prot] += rng.normal(scale=noise, size=(len(prot), 3))
    resid = np.array([a[0] for a in atoms])
    _, resindex = np.unique(resid, return_inverse=True)
    names = [a[2] for a in atoms]
    resn = np.array([a[1] for a in atoms])
    bonds = [(i, i + 1) for i in prot[:-1] if i + 1 in prot and resid[i] == resid[i + 1]]       # N-CA, CA-C, C-O
    by = {(resid[i], names[i]): i for i in prot}
    rids = sorted({resid[i] for i in prot})
    bonds += [(by[(a, "C")], by[(b, "N")]) for a, b in zip(rids, rids[1:])]                       # peptide bonds
    ca = [i for i in prot if names[i] == "CA"]
    v = (vel_rng or rng).normal(size=pos.shape) if vel else None
    counts = {}
    for r in np.unique(resid):
        nm = resn[np.flatnonzero(resid == r)[0]]
        counts[nm] = counts.get(nm, 0) + 1
    return {"n_atoms": len(atoms), "names": names, "resnames": resn, "resindex": resindex,
            "protein_resnames": set(resn[prot]), "protein_residues_outside_gro_set": [], "positions": pos,
            "dims": np.array([BOX_NM * 10] * 3 + [90.0] * 3), "velocities": v, "protein_bonds": np.array(bonds, np.int64),
            "ca_positions": pos[ca], "ca_resnames": list(resn[ca]), "heavy_whole": pos[prot],
            "protein_vel": None if v is None else v[prot], "residue_counts": dict(sorted(counts.items())),
            "net_charge": 0.0}


def build(tmp_path, *, mdp_override=None, ndx_override=None, na_nomyr=3, drop=(), bad_mdp=False):
    stage = tmp_path / "stage"
    for c, model in C.CONDITIONS.items():
        myr = c == "Myr"
        atoms = system_atoms(myr, n_na=3 if myr else na_nomyr)
        for s, reps in C.SETS.items():
            d = stage / model / s
            d.mkdir(parents=True, exist_ok=True)
            grp = "Protein_MYRP_ATP_MG" if myr else "Protein_ATP_MG"
            for f in C.MDP_FILES:
                (d / f).write_text((mdp_override or {}).get((s, c, f), MDP.format(grp=grp)))
            (d / "index.ndx").write_text((ndx_override or {}).get((s, c), ndx(atoms, myr)))
            (d / "topol.top").write_text("; SYNTHETIC\n")
            (d / "system_em.gro").write_text(gro(atoms))
            for r in reps:
                (d / str(r)).mkdir(exist_ok=True)
                (d / str(r) / "md.tpr").write_bytes(f"SYNTHETIC {c} {s} {r}".encode())
    if bad_mdp:
        (stage / "kinase_kinked_Myr/Set1/index.ndx").write_text("1 2 3\n[ X ]\n")
    data = tmp_path / "data"
    data.mkdir()
    orig = {"inputs": {}}
    for model in C.CONDITIONS.values():
        with zipfile.ZipFile(data / f"{model}.zip", "w") as zf:
            for p in sorted((stage / model).rglob("*")):
                rel = str(p.relative_to(stage / model))
                if p.is_file() and not any(rel.endswith(x) for x in drop):
                    zf.write(p, f"{model}/{rel}")
        b = (data / f"{model}.zip").read_bytes()
        orig["inputs"][f"{model}.zip"] = {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
    return orig, data


def loader_factory(*, fail=(), shift=None, reorder=(), copies=(), same_vel=False, charge=None, na_nomyr=3):
    """A SYNTHETIC loader keyed on the md.tpr bytes written by build()."""
    def load(path):
        tag = Path(path).read_bytes().decode().split()
        c, s, r = tag[1], tag[2], int(tag[3])
        key = f"{s}:{c}:{r}"
        if key in fail:
            raise OSError("SYNTHETIC unreadable")
        atoms = system_atoms(c == "Myr", n_na=3 if c == "Myr" else na_nomyr)
        run_seed = 1000 * int(s[-1]) + 100 * (c == "Myr") + r
        seed = 7 if (s, c) in copies else run_seed                     # copies share positions, not velocities
        rec = make_record(atoms, np.random.default_rng(seed), vel_rng=np.random.default_rng(run_seed + 50000))
        if same_vel and (s, c) in copies:
            rec["protein_vel"] = np.zeros_like(rec["protein_vel"])
        if shift and key in shift:
            rec["positions"] = rec["positions"].copy()
            rec["ca_positions"] = rec["ca_positions"] + np.linspace(0, shift[key], len(rec["ca_positions"]))[:, None]
        if key in reorder:
            rec["names"] = rec["names"][::-1]
        if charge and c in charge:
            rec["net_charge"] = charge[c]
        return rec
    return load


def run(tmp_path, loader=None, **kw):
    orig, data = build(tmp_path, **kw)
    res, ca = C.run_check(orig, data, tmp_path / "ext", BUDGET, loader=loader or loader_factory(), log=lambda *_: None)
    return res


# --- pure helpers --------------------------------------------------------------------------------------------------------

def test_constants_and_member_plan_match_the_spec():
    k = SPEC["constants"]
    assert (C.BOND_MIN, C.BOND_MAX, C.CA_RMSD_MAX, C.IDENT_RMSD) == (k["bond_min_A"], k["bond_max_A"], k["ca_rmsd_to_em_max_A"],
                                                                   k["identical_start_rmsd_A"])
    assert (C.SOL_REL_MAX, C.CHARGE_DIFF_MAX, C.REF_T) == (k["sol_relative_difference_max"], k["net_charge_difference_max_e"], k["ref_t_K"])
    plan = C.member_plan()
    assert all(len(v) == 2 * (len(C.MDP_FILES) + 3) + 20 for v in plan.values())
    assert any(p.endswith("abl_myr_1a/results/ABL-MYR-1a_audit_results.json") for p in C.PRESERVED)


def test_mdp_comparison_classifies_myr_group_names_and_other_differences():
    a = MDP.format(grp="Protein_MYRP_ATP_MG") + "compressed-x-grps = Protein_MYRP_ATP_MG\n"
    b = MDP.format(grp="Protein_ATP_MG") + "compressed-x-grps = Protein_ATP_MG\n"
    m = C.compare_mdp(a, b)
    assert {k: d["class"] for k, d in m["differences"].items()} == {"tc_grps": "inherent_group_name",
                                                                    "compressed_x_grps": "inherent_group_name"}
    assert m["potential_confounds"] == []
    m = C.compare_mdp(a, b.replace("ref_t = 300 300", "ref_t = 300 310") + "nstcalcenergy = 50\n")
    assert m["potential_confounds"] == ["nstcalcenergy", "ref_t"]
    assert C.strip_myr("Protein_MYRP_ATP_MG Water_NA_CL") == "Protein_ATP_MG Water_NA_CL"
    assert C.compare_mdp("x = Protein_ABC\n", "x = Protein\n")["potential_confounds"] == ["x"]   # not a MYRP-only difference


def test_ref_t_is_compared_numerically():
    assert C.ref_t_ok("300 300") and C.ref_t_ok("300.0 300")
    assert not C.ref_t_ok("300 310") and not C.ref_t_ok(None) and not C.ref_t_ok("") and not C.ref_t_ok("warm")


def test_index_parsing_and_errors():
    g = C.parse_ndx("[ A ]\n1 2 3\n4\n[ B ]\n5\n")
    assert g["A"].tolist() == [0, 1, 2, 3] and g["B"].tolist() == [4]
    with pytest.raises(ValueError, match="duplicate"):
        C.parse_ndx("[ A ]\n1\n[ A ]\n2\n")
    with pytest.raises(ValueError, match="before the first group"):
        C.parse_ndx("1 2\n[ A ]\n3\n")


def _sys(myr):
    atoms = system_atoms(myr)
    rec = make_record(atoms, np.random.default_rng(0))
    return atoms, rec


def test_tc_group_content_rule():
    for myr, cond in ((True, "Myr"), (False, "noMyr")):
        atoms, rec = _sys(myr)
        g = C.parse_ndx(ndx(atoms, myr))
        names = list(g)
        r = C.tc_group_check(g, names, rec["resnames"], rec["resindex"], rec["protein_resnames"], cond)
        assert r["ok"], r["problems"]
        assert r["groups"][names[1]]["residue_counts"] == {"CL": 2, "NA": 3, "SOL": 30}
        assert ("MYRP" in r["groups"][names[0]]["residue_counts"]) == myr
    atoms, rec = _sys(True)
    g = C.parse_ndx(ndx(atoms, True))
    r = C.tc_group_check(g, ["Protein_MYRP_ATP_MG", "Missing"], rec["resnames"], rec["resindex"], rec["protein_resnames"], "Myr")
    assert not r["ok"] and "not in index.ndx" in r["problems"][0]
    g2 = dict(g)
    g2["Water_NA_CL"] = g["Water_NA_CL"][:-1]                                       # one ion in no group
    r = C.tc_group_check(g2, list(g), rec["resnames"], rec["resindex"], rec["protein_resnames"], "Myr")
    assert not r["ok"] and "partition" in r["problems"][0]
    r = C.tc_group_check(g, list(g), rec["resnames"], rec["resindex"], rec["protein_resnames"], "noMyr")
    assert not r["ok"] and any("MYRP" in p for p in r["problems"])                  # caught by the content rule
    mixed = {"A": np.concatenate([g["Protein_MYRP_ATP_MG"], g["Water_NA_CL"][:3]]), "B": g["Water_NA_CL"][3:]}
    r = C.tc_group_check(mixed, ["A", "B"], rec["resnames"], rec["resindex"], rec["protein_resnames"], "Myr")
    assert not r["ok"] and any("SOL" in p for p in r["problems"])


def test_composition_rule():
    base = {"ATP": 1, "CL": 87, "MG": 1, "NA": 100, "SOL": 29938, "ALA": 10}
    ok = C.composition_check({**base, "MYRP": 1, "SOL": 29949}, base, 0.0, 0.0)
    assert ok["ok"] and len(ok["inherent"]) == 2 and ok["potential_confounds"] == []
    assert not C.composition_check({**base, "MYRP": 1, "NA": 101}, base, 0.0, 0.0)["ok"]
    assert not C.composition_check({**base, "SOL": 29938 + 400}, base, 0.0, 0.0)["ok"]                # 1.3 % > 1 %
    bad = C.composition_check({**base, "MYRP": 1}, base, -1.0, 0.0)
    assert not bad["ok"] and "net charge" in bad["potential_confounds"][0]


def test_kabsch_rmsd_is_rigid_motion_invariant_and_detects_mirrors():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(50, 3))
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.linalg.det(q))
    assert C.kabsch_rmsd(P, P @ q.T + 5.0) < 1e-9
    assert C.kabsch_rmsd(P, P * np.array([1, 1, -1])) > 0.1                          # a mirror image is not superposable
    Q = P.copy()
    Q[0] += np.array([3.0, 0, 0])
    assert 0 < C.kabsch_rmsd(P, Q) < 3.0 / np.sqrt(50) + 1e-9
    assert C.radius_of_gyration(np.array([[1.0, 0, 0], [-1.0, 0, 0]])) == pytest.approx(1.0)


def test_whole_chain_and_bonds_under_periodic_boundaries():
    dims = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    chain = np.array([[18.0 + 3.8 * i, 5.0, 5.0] for i in range(5)])
    wrapped = chain.copy()
    wrapped[:, 0] %= 20.0
    whole = C.whole_chain(wrapped, dims)
    assert np.allclose(np.linalg.norm(np.diff(whole, axis=0), axis=1), 3.8, atol=1e-4)
    b = C.bond_stats(wrapped, np.array([[i, i + 1] for i in range(4)]), dims)
    assert b["outside"] == 4 and b["min_A"] == pytest.approx(3.8, abs=1e-4)               # 3.8 Å exceeds the 2.3 Å limit
    pair = np.array([[19.5, 1.0, 1.0], [0.9, 1.0, 1.0]])                                 # 1.4 Å across the boundary
    assert C.bond_stats(pair, np.array([[0, 1]]), dims)["ok"]
    tri = np.array([50.0, 50.0, 50.0, 60.0, 60.0, 90.0])
    w = C.whole_chain(np.array([[1.0, 1.0, 1.0], [49.0, 1.0, 1.0]]), tri)
    assert np.linalg.norm(w[1] - w[0]) == pytest.approx(2.0, abs=1e-4)


def test_gro_box_parsing(tmp_path):
    (tmp_path / "a.gro").write_text("x\n0\n   5.00000   6.00000   7.00000\n")
    assert C.gro_box_dims(tmp_path / "a.gro").tolist() == [50.0, 60.0, 70.0, 90.0, 90.0, 90.0]
    (tmp_path / "b.gro").write_text("x\n0\n   5.0 5.0 3.53553 0.0 0.0 0.0 0.0 2.5 2.5\n")      # rhombic dodecahedron (xy-square)
    d = C.gro_box_dims(tmp_path / "b.gro")
    assert d[:3] == pytest.approx([50.0, 50.0, 50.0], abs=1e-3) and d[3:] == pytest.approx([60.0, 60.0, 90.0], abs=1e-3)
    (tmp_path / "c.gro").write_text("x\n0\n 1 2\n")
    with pytest.raises(ValueError):
        C.gro_box_dims(tmp_path / "c.gro")


def test_distinct_groups_and_decision():
    assert C.groups_distinct(4, []) == 4 and C.groups_distinct(4, [(0, 1), (1, 2), (0, 2), (3, 3)]) == 2
    ok = {k: True for k in ("K1_usable_coordinates", "K2_thermostat_comparable", "K3_tc_groups", "K4_composition",
                            "K5_replicates_not_copies")}
    assert C.decide(ok)["outcome"] == "comparable_for_1b"
    for k in ok:
        d = C.decide({**ok, k: False})
        assert d["outcome"] == "not_comparable" and d["failed"] == [k]


# --- run_check on SYNTHETIC archives ---------------------------------------------------------------------------------------

def test_comparable_case(tmp_path):
    r = run(tmp_path)
    assert r["outcome"] == "comparable_for_1b", r.get("decision")
    assert r["R1_usable_coordinates"]["runs_ok"] == 40
    g = r["R2_production_start_geometry"]["Set1"]
    assert g["within_Myr"]["ca_rmsd_A"]["n"] == 45 and g["between_Myr_noMyr"]["ca_rmsd_A"]["n"] == 100
    assert g["within_Myr"]["distinct_starts"] == 10 and g["within_Myr"]["distinct_velocity_sets"] == 10
    assert r["R2_production_start_geometry"]["cross_set"]["Myr"]["ca_rmsd_A"]["n"] == 100
    cls = r["R6_classification"]
    assert any("tc_grps" in x for x in cls["inherent_to_Myr_model"]) and any("MYRP present" in x for x in cls["inherent_to_Myr_model"])
    assert cls["potential_confound"] == []
    assert r["resources"]["download_bytes"] == 0
    md = C.render_report(r)
    assert "comparable_for_1b" in md and "40 of 40" in md and "no causal Myr claim" in md


@pytest.mark.parametrize("kw,failed", [
    (dict(loader=loader_factory(fail=("Set2:noMyr:25",))), "K1_usable_coordinates"),
    (dict(loader=loader_factory(shift={"Set1:Myr:11": 8.0})), "K1_usable_coordinates"),           # Cα far from the em reference
    (dict(loader=loader_factory(reorder=("Set1:noMyr:12",))), "K1_usable_coordinates"),
    (dict(loader=loader_factory(copies={("Set1", "Myr")}, same_vel=True)), "K5_replicates_not_copies"),
    (dict(loader=loader_factory(charge={"Myr": -1.0})), "K4_composition"),
    (dict(na_nomyr=4, loader=loader_factory(na_nomyr=4)), "K4_composition"),
    (dict(mdp_override={("Set2", "noMyr", "nvt.mdp"): MDP.format(grp="Protein_ATP_MG").replace("300 300", "310 310")}),
     "K2_thermostat_comparable"),
    (dict(mdp_override={("Set1", "Myr", "npt.mdp"): MDP.format(grp="Protein_MYRP_ATP_MG").replace("tau_t = 1.0 1.0", "tau_t = 0.5 0.5")}),
     "K2_thermostat_comparable"),
    (dict(ndx_override={("Set1", "Myr"): "[ Protein_MYRP_ATP_MG ]\n1 2\n[ Water_NA_CL ]\n3\n"}), "K3_tc_groups"),
])
def test_each_criterion_fails_on_its_synthetic_defect(tmp_path, kw, failed):
    r = run(tmp_path, **kw)
    assert r["outcome"] == "not_comparable" and failed in r["decision"]["failed"], r["decision"]


def test_identical_positions_with_distinct_velocities_are_not_copies(tmp_path):
    r = run(tmp_path, loader=loader_factory(copies={("Set1", "Myr")}))
    w = r["R2_production_start_geometry"]["Set1"]["within_Myr"]
    assert w["distinct_starts"] == 1 and w["distinct_velocity_sets"] == 10 and w["pairs_below_0.01A"] == 45
    assert r["outcome"] == "comparable_for_1b"
    assert r["decision"]["criteria"]["K5_replicates_not_copies"]


def test_stop_conditions(tmp_path):
    orig, data = build(tmp_path / "s1")
    orig["inputs"]["kinase_kinked_Myr.zip"]["sha256"] = "0" * 64
    r, _ = C.run_check(orig, data, tmp_path / "x", BUDGET, loader=loader_factory())
    assert r["stop"]["code"] == "S1" and r["outcome"] == "check_incomplete"
    assert run(tmp_path / "s2", drop=("Set2/index.ndx",))["stop"]["code"] == "S2"
    assert run(tmp_path / "s3", bad_mdp=True)["stop"]["code"] == "S3"
    orig, data = build(tmp_path / "s4")
    r, _ = C.run_check(orig, data, tmp_path / "y", {"id": "X", "budget": {"extracted_bytes_max": 10, "cpu_seconds_max": 1800}},
                       loader=loader_factory())
    assert r["stop"]["code"] == "S4"


def test_preserved_artifacts_unchanged_when_a_freeze_exists():
    frz = ROOT / "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_freeze.json"
    if not frz.exists():
        pytest.skip("check freeze not yet written")
    for p, h in json.loads(frz.read_text())["preserved_sha256"].items():
        assert C.R0.sha256_file(ROOT / p) == h, p

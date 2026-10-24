"""SYNTHETIC tests for the ABL-MYR-1 runner (benchmarks/audit/abl_myr_1/run_abl_myr.py). Engineering checks only.

No Dryad file is used or created, and nothing here is a scientific result. Pure tests (gates, statistics, parsing, archive
handling, the frozen-record checks) run everywhere. The end-to-end tests need both the compiled TopoFold core and an
MDAnalysis interpreter (set TOPOFOLD_ADAPTER_PYTHON) and are skipped otherwise; the MDAnalysis-only worker tests are in
tests/test_abl_myr_runner_mda.py.
"""

import hashlib
import importlib.util
import json
import os
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1"))
import run_abl_myr as R  # noqa: E402
from topofold_inputs.kinase_roles import STANDARD  # noqa: E402

THREE = {v: k for k, v in STANDARD.items()}
FIRST, LAST, AUTHOR_SHIFT = 285, 410, 19          # SYNTHETIC chain = isoform-1b residues 285-410; author = 1b - 19
SMALL = {**R.CFG, "sets": {"Set1": [11, 12, 13, 14], "Set2": [21, 22, 23, 24]}, "min_runs": 3, "min_frames": 30}
N_FRAMES = 40
REGION = range(40, 50)                              # residue indices given a 3x larger SYNTHETIC fluctuation in Myr runs
ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON")
HAVE_CORE = importlib.util.find_spec("topofold") is not None
e2e = pytest.mark.skipif(not (ADAPTER and HAVE_CORE), reason="needs the TopoFold core and TOPOFOLD_ADAPTER_PYTHON")


# --- SYNTHETIC structure and archive builders ---------------------------------------------------------------------------

def ref_sequence() -> str:
    text = (ROOT / R.CFG["reference"]["fasta"]).read_text()
    return "".join(ln.strip() for ln in text.splitlines() if not ln.startswith(">"))


def residues():
    seq = ref_sequence()
    return [(i - AUTHOR_SHIFT, THREE[seq[i - 1]]) for i in range(FIRST, LAST + 1)]


def base_backbone():
    """(n, 4, 3) SYNTHETIC N, CA, C, O on a helix (CA–CA 3.83 Å); each N lies near the C(i−1)–CA(i) segment so that the
    peptide bond C(i−1)–N(i) is below 1.6 Å (the geometry is not stereochemically realistic)."""
    n = LAST - FIRST + 1
    k = np.arange(n)
    ca = np.stack([2.3 * np.cos(np.radians(100) * k), 2.3 * np.sin(np.radians(100) * k), 1.5 * k], axis=1)
    u = np.zeros_like(ca)
    u[:-1] = ca[1:] - ca[:-1]
    u[-1] = u[-2]
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    c = ca + 1.52 * u
    nn = np.empty_like(ca)
    w0 = np.cross(u[0], np.array([0.0, 0.0, 1.0]))
    nn[0] = ca[0] - 1.0 * u[0] + 0.8 * w0 / np.linalg.norm(w0)
    for i in range(1, n):
        seg = ca[i] - c[i - 1]
        w = np.cross(seg, np.array([0.0, 0.0, 1.0]))
        w /= np.linalg.norm(w)
        nn[i] = c[i - 1] + (1.33 / np.linalg.norm(seg)) * seg + 0.8 * w          # off the line: φ stays well conditioned
    perp = np.cross(u, np.array([0.0, 0.0, 1.0]))
    perp /= np.linalg.norm(perp, axis=1, keepdims=True)
    return np.stack([nn, ca, c, c + 1.23 * perp], axis=1)


def frames_for(rng, sigma_res, n_frames=N_FRAMES):
    base = base_backbone()
    return base[None] + rng.normal(size=(n_frames,) + base.shape) * sigma_res[None, :, None, None]


def pdb_text(res, X0):
    names = ("N", "CA", "C", "O")
    lines, k = ["REMARK   1 SYNTHETIC backbone-only chain for ABL-MYR-1 runner tests"], 0
    for (resid, resname), xyz in zip(res, X0):
        for nm, p in zip(names, xyz):
            k += 1
            lines.append(f"ATOM  {k:5d} {nm:<4s} {resname:>3s} A{resid:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}"
                         f"  1.00 20.00           {nm[0]}")
    return "\n".join(lines) + "\nEND\n"


def gro_text(res, X0_A, with_myr, shift=0.0):
    atoms = []
    names = ("N", "CA", "C", "O")
    for (resid, resname), xyz in zip(res, X0_A):
        for nm, p in zip(names, xyz):
            atoms.append((resid, resname, nm, p / 10.0 + shift))
    if with_myr:
        atoms += [(500, "MYR", "C1", np.array([1.0, 1.0, 1.0])), (500, "MYR", "C2", np.array([1.1, 1.0, 1.0]))]
    atoms += [(600, "SOL", "OW", np.array([2.0, 2.0, 2.0]))]
    out = ["SYNTHETIC system", f"{len(atoms)}"]
    for i, (resid, rn, nm, p) in enumerate(atoms, start=1):
        out.append(f"{resid:5d}{rn:<5s}{nm:>5s}{i:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
    return "\n".join(out) + "\n   5.0 5.0 5.0\n"


MDP = "; SYNTHETIC production parameters\nintegrator = md\ntcoupl = v-rescale\nref_t = 300\ngen_vel = no\ncontinuation = yes\ndt = 0.002\n"


def top_text(with_myr, sol=5000):
    inc = '#include "MYR.itp"\n' if with_myr else ""
    mol = "MYR 1\n" if with_myr else ""
    return (f'; SYNTHETIC topology\n#include "protein.itp"\n{inc}\n[ system ]\nSYN\n\n[ molecules ]\n'
            f"Protein_chain_A 1\n{mol}SOL {sol}\n")


def write_model_files(d: Path, cond: str, *, shift=0.0, mdp=MDP, sol=5000):
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    res, X0 = residues(), base_backbone()
    (d / "system.gro").write_text(gro_text(res, X0, cond == "Myr", shift))
    (d / "topol.top").write_text(top_text(cond == "Myr", sol))
    (d / "production_md.mdp").write_text(mdp)
    return {"gro": d / "system.gro", "top": d / "topol.top", "mdp_prod": d / "production_md.mdp"}


def write_xtc(path: Path, pdb: Path, X, t0_ps=1.0e6, dt_ps=2000.0, times=None):
    import MDAnalysis as mda
    u = mda.Universe(str(pdb))
    n = u.atoms.n_atoms
    with mda.Writer(str(path), n_atoms=n) as w:
        for f in range(len(X)):
            u.atoms.positions = X[f].reshape(-1, 3)
            u.dimensions = [200, 200, 200, 90, 90, 90]
            u.trajectory.ts.time = (t0_ps + dt_ps * f) if times is None else times[f]
            w.write(u.atoms)


def build_dataset(root: Path, cfg=SMALL, *, break_runs=None, n_frames=N_FRAMES):
    """SYNTHETIC archives with the member layout of the frozen spec. ``break_runs``: {run key: residue indices} displaced
    by 15 Å in every 4th frame (SYNTHETIC chain breaks)."""
    import MDAnalysis  # noqa: F401  (XTC writing)
    root = Path(root)
    stage = root / "stage"
    res = residues()
    pdb = stage / "PCA" / "KD_noIa_bb_ref.pdb"
    pdb.parent.mkdir(parents=True, exist_ok=True)
    pdb.write_text(pdb_text(res, base_backbone()))
    (stage / "PCA" / "average.pdb").write_text(pdb.read_text())
    (stage / "PCA" / "steps").write_text("SYNTHETIC\n")
    # per-coordinate σ in Å, small enough that Cα–Cα distances stay below the adapter's 4.2 Å break threshold
    sig = {"Myr": np.where(np.isin(np.arange(len(res)), list(REGION)), 0.06, 0.02), "noMyr": np.full(len(res), 0.02)}
    for s, reps in cfg["sets"].items():
        for c in cfg["conditions"]:
            for r in reps:
                rng = np.random.default_rng(1000 * int(s[-1]) + 10 * (c == "Myr") + r)
                X = frames_for(rng, sig[c], n_frames)
                for i in (break_runs or {}).get(f"{s}:{c}:{r}", ()):
                    X[::4, i] += 15.0
                write_xtc(stage / "PCA" / f"traj_kinked_{c}_{r}.xtc", pdb, X)
    for c, model in cfg["conditions"].items():
        for s, reps in cfg["sets"].items():
            base = stage / model / model / s
            write_model_files(base, c)
            (base / "system_em.gro").write_text("em\n")
            (base / "nvt.mdp").write_text(MDP)
            (base / "npt.mdp").write_text(MDP)
            for r in reps:
                (base / str(r)).mkdir(exist_ok=True)
                (base / str(r) / "md.tpr").write_bytes(f"SYNTHETIC tpr {c} {s} {r}".encode())
    data = root / "data"
    data.mkdir()
    spec = {"inputs": {}, "budget": {"cpu_hours_max": 2.0, "storage_bytes_max": 4_000_000_000}}
    members = {"PCA.zip": [stage / "PCA"]}
    for c, model in cfg["conditions"].items():
        members[f"{model}.zip"] = [stage / model / model]
    for arch, tops in members.items():
        with zipfile.ZipFile(data / arch, "w") as zf:
            for top in tops:
                for p in sorted(top.rglob("*")):
                    if p.is_file():
                        zf.write(p, str(p.relative_to(top.parent)))
        b = (data / arch).read_bytes()
        spec["inputs"][arch] = {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
    return spec, data


# --- pure tests: frozen records ----------------------------------------------------------------------------------------

def test_frozen_spec_is_unchanged_and_cfg_matches_its_text():
    frz = json.loads((ROOT / "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_freeze.json").read_text())
    assert R.sha256_file(R.SPEC) == frz["spec_sha256"] and frz["spec_sha256"].startswith("f8825840")
    text = R.SPEC.read_text()
    for needle in ("450", "fewer than 8 runs", "2 ns", "0.01 ns", "5 %", "1000 ns", "0.001 nm", "bins = 36", "max_lag = 50"):
        assert needle in text, needle
    assert (R.CFG["min_frames"], R.CFG["min_runs"], R.CFG["bins"], R.CFG["max_lag"]) == (450, 8, 36, 50)
    assert R.CFG["anchors"] == {"beta3_Lys": 290, "alphaC_Glu": 305, "HRD_H": 380, "DFG_Asp": 400}


def test_code_freeze_records_present_files_without_reading_them_and_refuses_after_extraction(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "HERE", tmp_path)
    monkeypatch.setattr(R, "CODE_FREEZE", tmp_path / "cf.json")
    data = tmp_path / "data_dryad"
    data.mkdir()
    f = data / "PCA.zip"
    f.write_bytes(b"x" * 10)
    f.chmod(0)                                    # contents cannot be read: the freeze must only list metadata
    try:
        assert R.cmd_freeze_code(None) == 0
    finally:
        f.chmod(0o600)
    rec = json.loads((tmp_path / "cf.json").read_text())
    assert [(q["name"], q["bytes"]) for q in rec["input_files_present_at_code_freeze"]] == [("PCA.zip", 10)]
    assert "sha256" not in rec["input_files_present_at_code_freeze"][0] and "not hashed" in rec["input_files_note"]
    assert rec["runner"]["sha256"] == R.sha256_file(R.SELF) and set(rec["tests_sha256"]) == set(R.TESTS)
    assert R.cmd_freeze_code(None) == 1           # an existing record is never overwritten
    (tmp_path / "cf.json").unlink()
    (data / "extracted").mkdir()
    with pytest.raises(SystemExit, match="already opened"):
        R.cmd_freeze_code(None)


def test_run_refuses_without_code_freeze(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "CODE_FREEZE", tmp_path / "missing.json")
    with pytest.raises(SystemExit, match="no code freeze"):
        R.verify_code_freeze()


# --- pure tests: G0 / G1 -----------------------------------------------------------------------------------------------

def _spec_for(tmp_path, content=b"abc"):
    spec = {"inputs": {}}
    for name in ("PCA.zip", "kinase_kinked_Myr.zip", "kinase_kinked_noMyr.zip"):
        (tmp_path / name).write_bytes(content)
        spec["inputs"][name] = {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
    return spec


def test_g0_accepts_matching_archives_and_names_each_failure(tmp_path):
    spec = _spec_for(tmp_path)
    assert set(R.check_g0(tmp_path, spec)) == set(spec["inputs"])
    (tmp_path / "PCA.zip").write_bytes(b"abd")                      # same size, different content
    with pytest.raises(R.GateError, match="SHA-256") as e:
        R.check_g0(tmp_path, spec)
    assert e.value.gate == "G0" and e.value.detail == {"archive": "PCA.zip"}
    (tmp_path / "PCA.zip").write_bytes(b"abcd")
    with pytest.raises(R.GateError, match="4 bytes, expected 3"):
        R.check_g0(tmp_path, spec)
    (tmp_path / "PCA.zip").unlink()
    with pytest.raises(R.GateError, match="not present"):
        R.check_g0(tmp_path, spec)


def test_member_list_matches_the_frozen_spec_counts():
    w = R.wanted_members(R.CFG)
    assert len(w["PCA.zip"]) == 3 + 40
    for m in ("kinase_kinked_Myr.zip", "kinase_kinked_noMyr.zip"):
        assert len(w[m]) == 2 * 6 + 20
    assert w["PCA.zip"]["xtc:Set2:noMyr:30"] == "PCA/traj_kinked_noMyr_30.xtc"
    assert w["kinase_kinked_Myr.zip"]["tpr:Set1:11"] == "Set1/11/md.tpr"


def test_member_selection_finds_nested_or_flat_layouts_and_refuses_missing_or_ambiguous():
    cfg = {**R.CFG, "sets": {"Set1": [11]}}
    nested = {"PCA.zip": ["PCA/KD_noIa_bb_ref.pdb", "PCA/average.pdb", "PCA/steps", "PCA/traj_kinked_Myr_11.xtc",
                          "PCA/traj_kinked_noMyr_11.xtc", "PCA/proj_kinked_Myr_11.xtc.xvg", "PCA/"]}
    for m in ("kinase_kinked_Myr", "kinase_kinked_noMyr"):
        nested[f"{m}.zip"] = [f"{m}/{m}/Set1/{f}" for f in ("system.gro", "system_em.gro", "topol.top", "production_md.mdp",
                                                              "nvt.mdp", "npt.mdp", "11/md.tpr", "11/nvt.tpr")]
    sel = R.select_members(nested, cfg)
    assert sel["PCA.zip"]["xtc:Set1:Myr:11"] == "PCA/traj_kinked_Myr_11.xtc"
    assert sel["kinase_kinked_Myr.zip"]["gro:Set1"] == "kinase_kinked_Myr/kinase_kinked_Myr/Set1/system.gro"
    assert sel["kinase_kinked_Myr.zip"]["tpr:Set1:11"].endswith("Set1/11/md.tpr")
    missing = {**nested, "PCA.zip": [n for n in nested["PCA.zip"] if "noMyr_11" not in n]}
    with pytest.raises(R.GateError, match="0 matches") as e:
        R.select_members(missing, cfg)
    assert e.value.gate == "G1"
    dup = {**nested, "PCA.zip": nested["PCA.zip"] + ["copy/PCA/traj_kinked_Myr_11.xtc"]}
    with pytest.raises(R.GateError, match="2 matches"):
        R.select_members(dup, cfg)
    # a differently named, look-alike member must not match
    with pytest.raises(R.GateError):
        R.resolve_member(["PCA/xtraj_kinked_Myr_11.xtc"], "PCA/traj_kinked_Myr_11.xtc")


def test_safe_extract_rejects_traversal_and_verifies_content(tmp_path):
    z = tmp_path / "a.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("PCA/ok.txt", b"hello")
        zf.writestr(zipfile.ZipInfo("../evil.txt"), b"x")
        zf.writestr(zipfile.ZipInfo("/abs.txt"), b"x")
    with zipfile.ZipFile(z) as zf:
        rec = R.safe_extract(zf, "PCA/ok.txt", tmp_path / "out")
        assert Path(rec["path"]).read_bytes() == b"hello" and rec["sha256"] == hashlib.sha256(b"hello").hexdigest()
        for bad in ("../evil.txt", "/abs.txt"):
            with pytest.raises(R.GateError, match="unsafe"):
                R.safe_extract(zf, bad, tmp_path / "out")
    assert not (tmp_path / "evil.txt").exists()


# --- pure tests: G2 / G4 / G5 ------------------------------------------------------------------------------------------

def test_topology_choice_and_atom_identity_checks():
    assert R.choose_topology(500, 400, 500) == "ref"
    assert R.choose_topology(480, 500, 500) == "avg"
    with pytest.raises(R.GateError, match="neither"):
        R.choose_topology(480, 490, 500)
    topo = {"n_atoms": 500, "names": ["N", "CA", "C", "O"], "residues": [[1, "", "ALA", None, True, True, True]]}
    facts = {"a": {"n_atoms": 500}, "b": {"n_atoms": 500}}
    assert R.check_g2(topo, facts)["phi_psi_comparator_available"] is True
    with pytest.raises(R.GateError, match="different atom counts"):
        R.check_g2(topo, {"a": {"n_atoms": 500}, "b": {"n_atoms": 499}})
    with pytest.raises(R.GateError, match="topology has 500 atoms, trajectories 496"):
        R.check_g2(topo, {"a": {"n_atoms": 496}})
    with pytest.raises(R.GateError, match="not backbone") as e:
        R.check_g2({**topo, "names": ["N", "CA", "CB"]}, facts)
    assert "CB" in e.value.message
    no_c = {**topo, "residues": [[1, "", "ALA", None, True, True, False]]}
    assert R.check_g2(no_c, facts)["phi_psi_comparator_available"] is False


def test_g4_accepts_identical_starts_and_names_each_difference(tmp_path):
    cfg = R.CFG
    mk = lambda tag, **kw: {c: write_model_files(tmp_path / tag / c, c, **kw.get(c, {}))  # noqa: E731
                            for c in ("Myr", "noMyr")}
    ok = R.check_g4_set(mk("ok"), cfg)
    assert ok["gro"]["n_protein_atoms"] == 4 * (LAST - FIRST + 1) and ok["gro"]["max_abs_dx_nm"] == 0.0
    assert ok["initialisation_recorded"]["gen_vel"] == "no" and ok["mdp_production_identical"]
    with pytest.raises(R.GateError, match="system.gro") as e:                       # 0.002 nm shift of every atom
        R.check_g4_set(mk("shift", noMyr={"shift": 0.002}), cfg)
    assert e.value.gate == "G4" and "0.002" in e.value.message
    within = R.check_g4_set(mk("tol", noMyr={"shift": 0.0009}), cfg)                # inside the 0.001 nm tolerance
    assert within["gro"]["ok"]
    with pytest.raises(R.GateError, match="production_md.mdp") as e:
        R.check_g4_set(mk("mdp", noMyr={"mdp": MDP.replace("gen_vel = no", "gen_vel = yes")}), cfg)
    assert any("gen_vel" in x for x in e.value.detail["differing_lines"])
    with pytest.raises(R.GateError, match="topol.top") as e:                        # water count differs: not only Myr
        R.check_g4_set(mk("sol", noMyr={"sol": 4990}), cfg)
    assert any("SOL" in x for x in e.value.detail["differing_lines"])
    # comments and whitespace in the mdp are ignored
    noisy = MDP.replace("integrator = md", "integrator   =   md   ; a comment") + "\n; trailing\n"
    assert R.check_g4_set(mk("noisy", noMyr={"mdp": noisy}), cfg)["mdp_production_identical"]


def test_g4_protein_atom_mismatch_is_reported(tmp_path):
    a = write_model_files(tmp_path / "a", "Myr")
    b = write_model_files(tmp_path / "b", "noMyr")
    lines = b["gro"].read_text().splitlines()
    kept = [lines[0], str(int(lines[1]) - 1)] + lines[2:6] + lines[7:]            # drop one protein atom
    b["gro"].write_text("\n".join(kept) + "\n")
    with pytest.raises(R.GateError, match="protein atom counts differ"):
        R.check_g4_set({"Myr": a, "noMyr": b}, R.CFG)


def test_frame_gate_window_interval_and_exclusion():
    cfg = SMALL
    t = 1.0e6 + 2000.0 * np.arange(60)
    g = R.frame_gate(t, cfg)
    assert (g["start_index"], g["n_analysis"], g["excluded"]) == (0, 60, False)
    t0 = 2000.0 * np.arange(540)                                                   # starts at 0 ns, 1.078 µs long
    g = R.frame_gate(t0, cfg)
    assert g["start_index"] == 500 and g["n_analysis"] == 40 and not g["excluded"]
    g = R.frame_gate(t0, {**cfg, "min_frames": 450})
    assert g["excluded"] and "40 analysis frames < 450" in g["cause"]
    g = R.frame_gate(2000.0 * np.arange(100), cfg)                                 # never reaches 1000 ns
    assert g["n_analysis"] == 0 and g["excluded"]
    with pytest.raises(R.GateError, match="interval") as e:
        R.frame_gate(np.r_[t[:20], t[20:] + 1000.0], cfg)
    assert e.value.gate == "G5" and e.value.detail["max_dt_ps"] == 3000.0
    R.frame_gate(np.r_[t[:20], t[20:] + 9.9], cfg)                                 # inside ± 0.01 ns
    with pytest.raises(R.GateError, match="interval"):
        R.frame_gate(np.r_[t[:20], t[20:] + 10.1], cfg)
    with pytest.raises(R.GateError, match="strictly increasing"):
        R.frame_gate(np.r_[t[:10], t[9], t[11:]], cfg)
    with pytest.raises(R.GateError, match="strictly increasing"):
        R.frame_gate(t[::-1], cfg)


# --- pure tests: statistics and outcome --------------------------------------------------------------------------------

def test_e_stat_matches_hand_computation_and_orders_separated_groups():
    d = lambda p, q: abs(p - q)  # noqa: E731
    a, b = [0.0, 1.0, 2.0], [10.0, 11.0, 12.0]
    # cross pair distances: 10,11,12,9,10,11,8,9,10 -> median 10; within a: 1,2,1 -> median 1; within b the same
    assert R.e_stat(a, b, d) == pytest.approx(10.0 - 1.0)
    assert R.e_stat(a, a, d) == pytest.approx(1.0 - 1.0)
    assert R.e_stat([0.0], b, d) is None
    assert R.e_stat(b, a, d) == R.e_stat(a, b, d)


def test_rows_to_e_and_scalar_and_hist_e():
    rows = [{"channel": "kappa", "window": 3, "between_member_js_median": 0.5, "within_member_js_median_a": 0.2,
             "within_member_js_median_b": 0.3},
            {"channel": "tau", "window": 3, "between_member_js_median": 0.5, "within_member_js_median_a": None,
             "within_member_js_median_b": 0.3}]
    assert R.rows_to_e(rows) == {("kappa", 3): pytest.approx(0.2)}
    vals = {"a1": np.array([1.0, 5.0]), "a2": np.array([1.2, 5.0]), "b1": np.array([2.0, np.nan]),
            "b2": np.array([2.2, 5.0])}
    e = R.scalar_e(vals, ["a1", "a2"], ["b1", "b2"])
    assert set(e) == {0} and e[0] == pytest.approx(1.0 - 0.2)                        # residue 1 has a NaN: skipped
    h = {"a1": {0: np.array([10, 0])}, "a2": {0: np.array([9, 1])}, "b1": {0: np.array([0, 10])}, "b2": {0: np.array([1, 9])}}
    assert R.hist_e(h, ["a1", "a2"], ["b1", "b2"])[0] > 0.5


def test_classify_outcome_covers_every_ordered_step():
    n = 100
    assert R.classify_outcome({}, n)["label"] == "inconclusive"
    assert R.classify_outcome({"kappa": {5}}, n, sets_analysed=False)["note"] == "a set was not analysed"
    o = R.classify_outcome({"kappa": {10, 30, 50}, "d_i_i2": set()}, n)                # three isolated TopoFold residues
    assert o["label"] == "TopoFold adds reproduced residues" and len(o["R_T_outside_N2_R_B"]) >= 3
    o = R.classify_outcome({"kappa": {10, 30, 50}, "rmsf": {10, 30, 50}}, n)
    assert o["label"] == "baselines cover TopoFold" and o["R_T_outside_N2_R_B"] == []
    o = R.classify_outcome({"kappa": {10}, "rmsf": {40}}, n)                          # κ window 10 = residues 9–11, all outside
    assert o["label"] == "TopoFold adds reproduced residues" and o["R_T_outside_N2_R_B"] == [9, 10, 11]
    o = R.classify_outcome({"kappa": {10}, "rmsf": {13}}, n)                          # N2 covers 11–15: residues 9, 10 outside
    assert o["label"] == "partial" and o["R_T_outside_N2_R_B"] == [9, 10]
    o = R.classify_outcome({"kappa": {10}, "d_i_i3": {10}}, n)                        # τ-span of the baseline covers it
    assert o["label"] == "baselines cover TopoFold"
    o = R.classify_outcome({"rmsf": {40}}, n)                                         # baselines only: R_T empty
    assert o["label"] == "baselines cover TopoFold" and o["R_B_outside_N2_R_T"] == [40]
    o = R.classify_outcome({"phi": {40}}, n)
    assert o["label"] == "baselines cover TopoFold" and o["R_phi_psi"] == [40]       # comparator only, reported separately
    assert R.window_residues("tau", 10) == range(9, 13) and R.window_residues("kappa", 10) == range(9, 12)


def test_g6_fractions_and_temporal_summary():
    rows = {"distribution": [{"channel": "kappa", "window": w} for w in range(8)] * 2,
            "exclusions": [{"channel": "kappa", "window": 8, "cause": "x"}, {"channel": "kappa", "window": 8, "cause": "y"},
                           {"channel": "tau", "window": 1, "cause": "x"}]}
    f = R.g6_fractions(rows, {"phi": 98}, {"phi": 2})
    assert f["kappa"] == {"windows": 9, "excluded": 1, "fraction": pytest.approx(1 / 9)}
    assert f["tau"]["excluded"] == 1 and f["tau"]["windows"] == 1 and f["d_i_i2"]["windows"] == 0
    assert f["phi"]["fraction"] == pytest.approx(0.02)
    t = R.temporal_summary([{"channel": "kappa", "lag1": 0.5, "decorrelation_lag": 3},
                            {"channel": "kappa", "lag1": 0.9, "decorrelation_lag": None}])
    assert t["kappa"]["lag1_median"] == pytest.approx(0.7) and t["kappa"]["fraction_decorrelated_within_max_lag"] == 0.5


def test_render_report_and_path_scrubbing():
    res = {"status": "complete", "outcome": {"label": "partial", "R_T": [1], "R_B": [2], "R_T_outside_N2_R_B": [1],
                                              "R_B_outside_N2_R_T": [2], "R_phi_psi": []}, "sets": {}}
    assert "partial" in R.render_report(res) and "not biological validation" in R.render_report(res)
    s = R.scrub({"a": f"{R.ROOT}/x", "b": [str(Path.home()) + "/y"]})
    assert s == {"a": "<repo>/x", "b": ["~/y"]}


# --- end-to-end on SYNTHETIC archives (core + MDAnalysis interpreter) ------------------------------------------------

def build_dataset_any(root: Path, **kw):
    """build_dataset in this process if MDAnalysis is importable, otherwise through the adapter interpreter."""
    if importlib.util.find_spec("MDAnalysis"):
        return build_dataset(root, **kw)
    import subprocess
    code = ("import sys, json; sys.path.insert(0, sys.argv[1]); import test_abl_myr_runner as T; "
            "spec, data = T.build_dataset(sys.argv[2], break_runs=json.loads(sys.argv[3])); print(json.dumps(spec))")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT / "python"), os.environ.get("PYTHONPATH", "")]))
    p = subprocess.run([ADAPTER, "-c", code, str(Path(__file__).parent), str(root), json.dumps(kw.get("break_runs") or {})],
                       env=env, capture_output=True, text=True, check=True)
    return json.loads(p.stdout.strip().splitlines()[-1]), Path(root) / "data"


def _run(tmp_path, **kw):
    spec, data = build_dataset_any(tmp_path, **kw)
    return R.run_pipeline(spec, SMALL, data, tmp_path / "out", ADAPTER, work=tmp_path / "work", log=lambda *_: None)


@e2e
def test_end_to_end_on_synthetic_archives(tmp_path):
    res = _run(tmp_path)
    assert res["status"] == "complete" and set(res["gates"]) == {"G0", "G1", "G2", "G3", "G4", "G5", "G6"}
    assert res["gates"]["G3"]["alignment"]["author_numbering_diagnostic"]["reference_index_minus_author_resid"] == {
        str(AUTHOR_SHIFT): LAST - FIRST + 1}
    assert res["gates"]["G6"]["pbc_used"] == "none" and res["gates"]["G5"]["excluded_runs"] == []
    assert all(s["exclusion_rows"] == 0 for s in res["sets"].values())
    n = LAST - FIRST + 1
    for s, d in res["sets"].items():
        assert set(d["runs"]) == {"Myr", "noMyr"} and all(len(v) == 4 for v in d["runs"].values())
        for c in R.CHANNELS_ALL:
            assert d["channels"][c]["windows_evaluated"] >= n - 4
        assert d["channels"]["rmsf"]["windows_exceeding"] > 0                      # SYNTHETIC effect is detectable
        assert d["temporal"]["kappa"]["series"] > 0
    rep = res["outcome"]["reproduced_windows"]
    inside = set(range(min(REGION) - 2, max(REGION) + 3))
    assert len([w for w in rep["rmsf"] if w in set(REGION)]) >= 8                 # the injected SYNTHETIC difference
    assert len([w for w in rep["rmsf"] if w not in inside]) <= 5                  # SYNTHETIC noise only
    assert res["outcome"]["label"] in {"TopoFold adds reproduced residues", "baselines cover TopoFold", "partial"}
    again = _run(tmp_path / "again")
    assert json.dumps(again["outcome"], sort_keys=True) == json.dumps(res["outcome"], sort_keys=True)


@e2e
def test_end_to_end_reports_exclusions_from_synthetic_breaks(tmp_path):
    res = _run(tmp_path, break_runs={"Set1:Myr:11": [90]})
    d = res["sets"]["Set1"]
    assert d["exclusion_rows"] > 0 and "ca_ca_distance" in d["exclusion_causes"]
    frac = res["gates"]["G6"]["passes"][0]["g6_fractions"]["Set1"]
    assert 0 < frac["kappa"]["fraction"] <= 0.05 and frac["phi"]["excluded"] >= 2   # φ/ψ invalid around the displaced residue
    assert res["sets"]["Set2"]["exclusion_rows"] == 0


@e2e
def test_end_to_end_g6_stops_when_more_than_five_percent_is_excluded(tmp_path):
    brk = {"Set1:noMyr:12": list(range(5, 120, 8))}                               # 15 displaced residues
    with pytest.raises(R.GateError, match="more than 5%") as e:
        _run(tmp_path, break_runs=brk)
    assert e.value.gate == "G6"
    assert [q["pbc"] for q in e.value.detail["passes"]] == ["none", "make_whole_sequential"]


@e2e
def test_end_to_end_stops_at_g4_when_the_start_differs(tmp_path):
    spec, data = build_dataset_any(tmp_path)
    # rebuild the noMyr archive with a changed starting structure; update the recorded hash like a re-download would
    z = data / "kinase_kinked_noMyr.zip"
    with zipfile.ZipFile(z) as zf:
        items = {n: zf.read(n) for n in zf.namelist()}
    key = next(n for n in items if n.endswith("Set2/system.gro"))
    items[key] = gro_text(residues(), base_backbone() + 0.5, False).encode()
    with zipfile.ZipFile(z, "w") as zf:
        for n, b in items.items():
            zf.writestr(n, b)
    spec["inputs"]["kinase_kinked_noMyr.zip"] = {"bytes": z.stat().st_size, "sha256": R.sha256_file(z)}
    with pytest.raises(R.GateError, match="system.gro") as e:
        R.run_pipeline(spec, SMALL, data, tmp_path / "out", ADAPTER, work=tmp_path / "work", log=lambda *_: None)
    assert e.value.gate == "G4"

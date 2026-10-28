"""MET wild-type template feature exporter: synthetic tests (CI) and real-input checks (skipped where the git-ignored PDB files are absent).

SYNTHETIC structures are generated here; no scientific result is computed.
"""
import builtins
import hashlib
import importlib.util
import json
import math
import random
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "benchmarks/audit/met_feature_export_2026-10-27"
spec = importlib.util.spec_from_file_location("met_export", D / "met_export.py")
me = importlib.util.module_from_spec(spec)
sys.modules["met_export"] = me
spec.loader.exec_module(me)

CONTRACT = me.load_contract()
REF = me.read_fasta(D / "data/P08581.fasta")
DRUG = {"code": "SYN", "pdb": "SYN1", "inhibitor_component": "LIG", "article_class": "synthetic"}
LO, HI = 1055, 1349


def _line(rec, serial, name, alt, resname, chain, resseq, icode, xyz, occ=1.0, element="C"):
    return f"{rec:<6}{serial:>5} {name:<4}{alt}{resname:>3} {chain}{resseq:>4}{icode}   {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{occ:6.2f}{20.0:6.2f}          {element:>2}"


def make_structure(remove=(), rename=None, drop_atoms=None, icodes=(), altlocs=(), shift_ca=None, ligand=True, chains=("A",), lig_copies=1, ca_override=None, hetero_extra=True):
    """SYNTHETIC chain(s) over LO..HI built from the pinned reference: helix-like Cα trace with jitter, standard heavy atoms placed near each Cα."""
    rename, drop_atoms, shift_ca = rename or {}, drop_atoms or {}, shift_ca or {}
    rng = random.Random(7)
    lines, serial = [], 1
    cas = {}
    for q in range(LO, HI + 1):
        a = math.radians(100.0 * (q - LO))
        cas[q] = np.array([2.3 * math.cos(a), 2.3 * math.sin(a), 1.5 * (q - LO)]) + np.array([rng.uniform(-.15, .15) for _ in range(3)])
    for q, v in shift_ca.items():
        cas[q] = cas[q] + np.array(v)
    for q, v in (ca_override or {}).items():
        cas[q] = np.array(v)
    for ch in chains:
        for q in range(LO, HI + 1):
            if q in remove:
                continue
            wt = REF[q - 1]
            rn = rename.get(q, me.ONE_TO_THREE[wt])
            names = list(me.BACKBONE) + list(me.SIDE[wt])
            ic = "A" if q in icodes else " "
            r2 = random.Random(q)
            for n in names:
                if n in drop_atoms.get(q, ()):
                    continue
                xyz = cas[q] if n == "CA" else cas[q] + np.array([r2.uniform(-2.5, 2.5) for _ in range(3)])
                if q in altlocs:
                    lines.append(_line("ATOM", serial, n, "A", rn, ch, q, ic, xyz, 0.6, n[0])); serial += 1
                    lines.append(_line("ATOM", serial, n, "B", rn, ch, q, ic, xyz + np.array([1.0, 0, 0]), 0.4, n[0])); serial += 1
                else:
                    lines.append(_line("ATOM", serial, n, " ", rn, ch, q, ic, xyz, 1.0, n[0])); serial += 1
    if ligand:
        for c in range(lig_copies):
            for i in range(6):
                lines.append(_line("HETATM", serial, f"C{i + 1}", " ", "LIG", chains[0], 2000 + c, " ", cas[1200] + np.array([3.0 + 0.7 * i, 0.5 * c, -1.0 + 0.4 * i]), 1.0, "C")); serial += 1
    if hetero_extra:
        lines.append(_line("HETATM", serial, "C1", " ", "GBL", chains[0], 2100, " ", cas[1201] + np.array([2.0, 2.0, 2.0]), 1.0, "C")); serial += 1
        lines.append(_line("HETATM", serial, "O", " ", "HOH", chains[0], 2200, " ", cas[1202] + np.array([0.5, 0.5, 0.5]), 1.0, "O")); serial += 1
    return "\n".join(lines) + "\nEND\n"


def entry(**kw):
    return me.build_entry(DRUG, make_structure(**kw), REF, CONTRACT)


def rows_by_pos(e):
    return {r["uniprot_position"]: r for r in e["rows"]}


def feat(r, cols):
    return np.array([r[c] for c in cols], float)


T_ALL = me.KAPPA_COLS + me.TAU_COLS
# ---------------------------------------------------------------------------- universe, mutation properties, contract


def test_row_universe_is_generated_from_the_reference_alone():
    u = me.row_universe(REF, CONTRACT)
    assert len(u) == 287 * 19 * 6 == 32718 and len({r["row_id"] for r in u}) == len(u)
    assert all(REF[r["uniprot_position"] - 1] == r["wt"] != r["mut"] for r in u) and {r["uniprot_position"] for r in u} == set(range(1059, 1346))
    assert hashlib.sha256((D / "data/P08581.fasta").read_bytes()).hexdigest() == CONTRACT["reference"]["sha256"]
    assert all(r["dms_position"] == r["uniprot_position"] - 1058 for r in u)


def test_mutation_properties_are_pinned_complete_and_antisymmetric():
    bl = me.read_blosum62(D / "data/BLOSUM62.txt")
    assert bl[("A", "R")] == -1 and bl[("W", "W")] == 11 and bl[("A", "A")] == 4
    m = me.mutation_table(bl)
    assert len(m) == 380 and len({(r["wt"], r["mut"]) for r in m}) == 380
    by = {(r["wt"], r["mut"]): r for r in m}
    for (a, b), r in by.items():
        s = by[(b, a)]
        assert r["dvol"] == pytest.approx(-s["dvol"]) and r["dkd"] == pytest.approx(-s["dkd"]) and r["dcharge"] == -s["dcharge"] and r["blosum62"] == s["blosum62"]
    assert by[("D", "K")]["dcharge"] == 2 and by[("G", "P")]["to_pro"] == 1 and by[("P", "G")]["from_pro"] == 1 and by[("A", "G")]["dvol"] == pytest.approx(60.1 - 88.6)


# ---------------------------------------------------------------------------- mapping, ligand selection, chain rule


def test_a_shifted_numbering_is_caught_as_engineered_residues():
    """A structure numbered one residue off the reference fails identity everywhere: no feature survives."""
    txt = make_structure()
    lines = []
    for ln in txt.splitlines():
        if ln.startswith("ATOM"):
            ln = ln[:22] + f"{int(ln[22:26]) + 1:>4}" + ln[26:]
        lines.append(ln)
    e = me.build_entry(DRUG, "\n".join(lines), REF, CONTRACT)
    bad = sum(r["T_reason"] != "ok" for r in e["rows"])
    assert bad > 0.9 * len(e["rows"]) and len(e["flagged_residues"]["engineered"]) > 200


def test_ligand_is_the_named_component_in_the_chosen_chain_only():
    e = entry()
    assert e["ligand"]["component"] == "LIG" and e["ligand"]["heavy_atoms"] == 6 and e["ligand"]["chain"] == "A"
    with pytest.raises(ValueError):
        entry(ligand=False)
    with pytest.raises(ValueError):
        entry(lig_copies=2)


def test_chain_rule_prefers_the_chain_with_more_resolved_residues_and_ties_go_to_the_first_letter():
    full = make_structure(chains=("B", "A"))
    res, order = me.parse_pdb(full)
    assert me.choose_chain(res, order, (LO, HI)) == "A"
    two = make_structure(chains=("A", "B"))
    res, order = me.parse_pdb(two)
    assert me.choose_chain(res, order, (LO, HI)) == "A"
    # chain A loses residues: B (complete) wins
    lines = [l for l in make_structure(chains=("A", "B")).splitlines() if not (l.startswith("ATOM") and l[21] == "A" and 1100 <= int(l[22:26]) <= 1110)]
    res, order = me.parse_pdb("\n".join(lines))
    assert me.choose_chain(res, order, (LO, HI)) == "B"


# ---------------------------------------------------------------------------- missingness propagation


def test_missing_residue_masks_every_window_that_contains_it_not_only_the_central_one():
    r = rows_by_pos(entry(remove={1100}))
    for p in range(1090, 1111):
        covered = 1096 <= p <= 1104
        for col in ("T_reason", "C_reason"):
            assert (r[p][col] != "ok") == covered, (p, col)
        assert np.isnan(feat(r[p], T_ALL + me.DIST_COLS)).all() if covered else np.isfinite(feat(r[p], T_ALL + me.DIST_COLS)).all()
    assert r[1100]["T_reason"].startswith("missing_residue") and r[1100]["L_reason"] == "missing_residue:1100" and r[1100]["residue_status"] == "missing"
    assert r[1099]["L_reason"] in ("ok", "no_side_chain")  # ligand context of a neighbour is unaffected


def test_engineered_residue_masks_its_windows_and_its_own_ligand_context():
    wt = REF[1119]
    other = "ALA" if wt != "A" else "GLY"
    r = rows_by_pos(entry(rename={1120: other}))
    for p in range(1110, 1131):
        assert (r[p]["T_reason"] != "ok") == (1116 <= p <= 1124), p
        assert (r[p]["C_reason"] != "ok") == (1116 <= p <= 1124), p
    assert r[1120]["L_reason"] == "engineered_residue:1120" and r[1119]["L_reason"] != "engineered_residue:1120"
    assert "engineered_residue:1120" in r[1120]["T_reason"]


def test_incomplete_side_chain_masks_ligand_context_of_that_residue_only():
    q = next(p for p in range(1130, 1200) if REF[p - 1] not in "GA")
    r = rows_by_pos(entry(drop_atoms={q: (me.SIDE[REF[q - 1]][-1],)}))
    assert r[q]["L_reason"] == f"incomplete_side_chain:{q}" and np.isnan(feat(r[q], me.L_COLS)).all()
    assert r[q]["T_reason"] == "ok" and r[q - 1]["L_reason"] in ("ok", "no_side_chain")


def test_chain_break_and_insertion_code_mask_all_windows_containing_them():
    r = rows_by_pos(entry(shift_ca={1140: (0.0, 0.0, 6.0)}))
    masked = [p for p in range(1130, 1151) if r[p]["T_reason"] != "ok"]
    assert masked == list(range(1136, 1145)) and "chain_break" in r[1140]["T_reason"]
    r = rows_by_pos(entry(icodes={1150}))
    assert [p for p in range(1140, 1161) if r[p]["T_reason"] != "ok"] == list(range(1146, 1155)) and r[1150]["L_reason"] == "insertion_code:1150"


def test_altloc_uses_one_conformer_blank_then_A_and_is_flagged_not_masked():
    base = rows_by_pos(entry())
    alt = rows_by_pos(entry(altlocs={1160}))
    for p in range(1156, 1165):
        assert alt[p]["T_reason"] == "ok" and alt[p]["altloc_in_window"] == 1
        np.testing.assert_allclose(feat(alt[p], me.DIST_COLS), feat(base[p], me.DIST_COLS), atol=1e-2)  # 'A' coordinates equal the unaltered ones (0.001 A text rounding)
    assert alt[1170]["altloc_in_window"] == 0


def test_glycine_has_an_explicit_no_side_chain_reason_and_other_features():
    g = next(p for p in range(1059, 1346) if REF[p - 1] == "G")
    r = rows_by_pos(entry())[g]
    assert r["L_reason"] == "no_side_chain" and math.isnan(r["L_min_sidechain_dist_A"]) and math.isfinite(r["L_min_heavy_dist_A"]) and r["T_reason"] == "ok"


# ---------------------------------------------------------------------------- exact coordinate support


def _window_table(ca):
    return {q: {"status": "ok", "ca": ca[i], "break_after": False, "altloc": False} for i, q in enumerate(range(100, 100 + len(ca)))}


@pytest.mark.parametrize("j", range(9))
def test_each_angle_and_distance_depends_on_exactly_its_atoms(j):
    rng = np.random.default_rng(j)
    ca = np.cumsum(rng.normal(size=(9, 3)) * 1.0 + np.array([3.0, 0.5, 0.2]), axis=0)
    t0 = _window_table(ca)
    k0, tau0, d0, reason, _, _ = me.window_block(t0, 104)
    assert reason == "ok"
    ca2 = ca.copy()
    ca2[j] += np.array([0.7, -0.4, 0.9])
    k1, tau1, d1, _, _, _ = me.window_block(_window_table(ca2), 104)
    for v in range(1, 8):  # kappa at vertex v uses atoms v-1, v, v+1
        assert (abs(k1[v - 1] - k0[v - 1]) > 1e-9) == (j in (v - 1, v, v + 1)), (j, "kappa", v)
    for h in range(1, 7):  # tau at hinge h uses atoms h-1 .. h+2
        changed = abs(np.sin(tau1[h - 1]) - np.sin(tau0[h - 1])) > 1e-9 or abs(np.cos(tau1[h - 1]) - np.cos(tau0[h - 1])) > 1e-9
        assert changed == (j in (h - 1, h, h + 1, h + 2)), (j, "tau", h)
    for n, (a, b) in enumerate(zip(*me.IU)):
        assert (abs(d1[n] - d0[n]) > 1e-9) == (j in (a, b)), (j, "dist", a, b)


def test_a_coordinate_outside_the_window_does_not_change_the_features():
    base, moved = rows_by_pos(entry()), rows_by_pos(entry(shift_ca={1200: (0.3, 0.2, 0.1)}))
    for p in range(1185, 1216):
        same = np.allclose(feat(base[p], T_ALL + me.DIST_COLS), feat(moved[p], T_ALL + me.DIST_COLS), atol=0)
        assert same == (not (1196 <= p <= 1204)), p


# ---------------------------------------------------------------------------- independent calculation, invariance


def test_distances_equal_independent_calculations():
    spd = pytest.importorskip("scipy.spatial.distance")
    e = entry()
    r = rows_by_pos(e)[1170]
    ca = np.array([[float(x) for x in ln[30:54].replace("-", " -").split()] for ln in make_structure().splitlines() if ln.startswith("ATOM") and ln[12:16].strip() == "CA" and 1166 <= int(ln[22:26]) <= 1174])
    np.testing.assert_allclose(feat(r, me.DIST_COLS), spd.pdist(ca), atol=2e-3)
    # ligand context by a plain-python loop
    res, order = me.parse_pdb(make_structure())
    lig = [res[k] for k in order if res[k]["resname"] == "LIG"][0]
    lpts = [me.atom(lig, n) for n in lig["atoms"]]
    q = 1200
    rq = res[("A", q, " ")]
    pts = {n: me.atom(rq, n) for n in list(me.BACKBONE) + list(me.SIDE[REF[q - 1]])}
    best = min(math.dist(p, l) for p in pts.values() for l in lpts)
    cnt = sum(1 for p in pts.values() if min(math.dist(p, l) for l in lpts) <= 4.5)
    rr = rows_by_pos(e)[q]
    assert rr["L_min_heavy_dist_A"] == pytest.approx(best) and rr["L_n_contact_atoms_4p5"] == cnt and rr["L_ca_min_dist_A"] == pytest.approx(min(math.dist(pts["CA"], l) for l in lpts))


def test_rigid_motion_leaves_features_unchanged_and_a_mirror_flips_torsion_sign_only():
    base = rows_by_pos(entry())
    rng = np.random.default_rng(3)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    t = np.array([5.0, -7.0, 3.0])

    def transform(txt, M):
        out = []
        for ln in txt.splitlines():
            if ln[:6] in ("ATOM  ", "HETATM"):
                v = M(np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
                ln = ln[:30] + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + ln[54:]
            out.append(ln)
        return "\n".join(out)

    txt = make_structure()
    rig = {r["uniprot_position"]: r for r in me.build_entry(DRUG, transform(txt, lambda v: Q @ v + t), REF, CONTRACT)["rows"]}
    mir = {r["uniprot_position"]: r for r in me.build_entry(DRUG, transform(txt, lambda v: np.array([-v[0], v[1], v[2]])), REF, CONTRACT)["rows"]}
    for p in (1100, 1170, 1250):
        cols = me.KAPPA_COLS + me.TAU_COLS + me.DIST_COLS + me.L_COLS
        np.testing.assert_allclose(feat(rig[p], cols), feat(base[p], cols), atol=2e-3)
        np.testing.assert_allclose(feat(mir[p], me.KAPPA_COLS + me.DIST_COLS + me.L_COLS), feat(base[p], me.KAPPA_COLS + me.DIST_COLS + me.L_COLS), atol=2e-3)
        sin_c = [c for c in me.TAU_COLS if "sin" in c]
        cos_c = [c for c in me.TAU_COLS if "cos" in c]
        np.testing.assert_allclose(feat(mir[p], sin_c), -feat(base[p], sin_c), atol=2e-3)
        np.testing.assert_allclose(feat(mir[p], cos_c), feat(base[p], cos_c), atol=2e-3)


# ---------------------------------------------------------------------------- identical descriptors across substitutions; outcome independence


def test_structural_descriptors_are_identical_across_the_19_substitutions_and_mutation_features_across_positions_and_drugs():
    e = entry()
    structure = [dict(r, drug="Crizo") for r in e["rows"]]
    universe = [u for u in me.row_universe(REF, CONTRACT) if u["drug"] == "Crizo"]
    mutation = me.mutation_table(me.read_blosum62(D / "data/BLOSUM62.txt"))
    flat = list(me.join_flat(universe, mutation, structure))
    assert len(flat) == 287 * 19
    s_cols = [c for c in structure[0] if c not in ("drug", "uniprot_position", "dms_position", "wt")]
    groups = {}
    for r in flat:
        groups.setdefault(r["uniprot_position"], []).append(r)
    for p, rs in groups.items():
        assert len(rs) == 19
        for c in s_cols:
            vals = [r[c] for r in rs]
            assert all((v == vals[0]) or (isinstance(v, float) and math.isnan(v) and math.isnan(vals[0])) for v in vals), (p, c)
    mcols = [c for c in mutation[0] if c not in ("wt", "mut")]
    by = {}
    for r in flat:
        by.setdefault((r["wt"], r["mut"]), set()).add(tuple(r[c] for c in mcols))
    assert all(len(v) == 1 for v in by.values())


def test_exporter_cannot_see_outcome_files(tmp_path, monkeypatch):
    src = (D / "met_export.py").read_text().lower()
    for token in ("rosace", "dmso", ".tsv", "met_scores", "wt_rosace"):
        assert token not in src
    outcomes = {"decoy_a.tsv": "variants\teffect\np.(A1059C)\t0.5\n", "decoy_b.tsv": "variants\teffect\np.(A1059C)\t-9\n"}
    opened = []
    real_open = builtins.open

    def spy(file, *a, **k):
        opened.append(str(file))
        return real_open(file, *a, **k)

    digests = []
    for scenario in ("absent", "decoy_a", "decoy_b"):
        wd = tmp_path / scenario
        wd.mkdir()
        if scenario != "absent":
            (wd / f"{scenario}.tsv").write_text(outcomes[f"{scenario}.tsv"])
            (wd / "outcomes.csv").write_text("row_id,y\nSYN:A1059C,1\n")
        monkeypatch.chdir(wd)
        monkeypatch.setattr(builtins, "open", spy)
        e = entry()
        monkeypatch.setattr(builtins, "open", real_open)
        digests.append(hashlib.sha256(json.dumps([{k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in r.items()} for r in e["rows"]], sort_keys=True).encode()).hexdigest())
    assert len(set(digests)) == 1 and not any("decoy" in f or "outcome" in f for f in opened)


# ---------------------------------------------------------------------------- real inputs (git-ignored files) and the freeze


REAL = all((D / d["input_file"]).is_file() for d in CONTRACT["drugs"])


@pytest.mark.skipif(not REAL, reason="git-ignored PDB files not available")
def test_real_entries_chain_ligand_and_engineered_residues_are_as_recorded():
    expected_engineered = {"3QTI": {1194: "PHE", 1234: "PHE"}}
    for d in CONTRACT["drugs"]:
        f = D / d["input_file"]
        assert hashlib.sha256(f.read_bytes()).hexdigest() == d["input_sha256"]
        e = me.build_entry(d, f.read_text(), REF, CONTRACT)
        assert len(e["rows"]) == 287 and e["ligand"]["component"] == d["inhibitor_component"]
        assert e["engineered_detail"] == expected_engineered.get(d["pdb"], {}), d["pdb"]
        assert not e["flagged_residues"]["insertion_code"], d["pdb"]
        assert e["flagged_residues"]["missing_CA"] == ({"2WGJ": [1345]}.get(d["pdb"], [])), d["pdb"]  # 2WGJ: residue 1345 present without a Cα


def test_the_freeze_record_matches_the_frozen_files_once_it_exists():
    fz = D / "freeze_exporter.json"
    if not fz.is_file():
        pytest.skip("not frozen yet")
    for f, h in json.loads(fz.read_text())["files"].items():
        if f.endswith("test_met_export.py"):
            continue  # this file's own hash is checked by `met_export.py run`, which would otherwise make the test self-referential
        assert hashlib.sha256((D / f).read_bytes()).hexdigest() == h, f

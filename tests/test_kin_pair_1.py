"""KIN-PAIR-1 helpers (benchmarks/audit/kin_pair_1): SYNTHETIC fixtures with known truth, analytic references, the selection
procedure on a mock RCSB. No network. The structpair-dependent test needs an interpreter with MDAnalysis and is skipped without."""

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KP_DIR = ROOT / "benchmarks/audit/kin_pair_1"
sys.path.insert(0, str(KP_DIR))
import kinpair as kp  # noqa: E402
import select_panel as sel  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
THREE = {v: k for k, v in kp.THREE2ONE.items() if len(k) == 3 and k not in ("MSE", "HYP", "CME", "HID", "HIE", "HIP", "HSD", "HSE", "HSP")}


def ref_sequence(n=300, seed=1):
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list(AA), n))


def helix(n, shift=0.0):
    t = np.arange(n) * math.radians(100)
    return np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)] + shift


def pdb_text(seq, first_resseq, start_idx, end_idx, chain="A", acc="P99999", db_begin=None, skip=(), seqadv=(), coords=None,
             tag="", altloc_at=None):
    """SYNTHETIC chain: residues start_idx..end_idx (1-based UniProt indices) numbered from first_resseq, with DBREF."""
    n = end_idx - start_idx + 1
    X = helix(n) if coords is None else coords
    lines = [f"DBREF  XXXX {chain} {first_resseq:4d} {first_resseq + n - 1:5d}  UNP    {acc:<8s}   SYN_HUMAN {start_idx if db_begin is None else db_begin:6d} {end_idx:6d}"]
    for (name, num, conflict) in seqadv:
        lines.append(f"SEQADV XXXX {name:>3s} {chain} {num:4d}  UNP  {acc:<9s} ALA {num:5d} {conflict}")
    k = 0
    for j in range(n):
        if j in skip:
            continue
        resseq = first_resseq + j
        name = THREE[seq[start_idx - 1 + j]]
        u = (X[min(j + 1, n - 1)] - X[max(j - 1, 0)])
        u = u / np.linalg.norm(u)
        for an, p in (("N", X[j] - 1.25 * u), ("CA", X[j]), ("C", X[j] + 1.25 * u)):
            alt = " "
            if altloc_at == j and an == "CA":
                alt = "A"
            k += 1
            lines.append(f"ATOM  {k:5d} {an:^4s}{alt}{name} {chain}{resseq:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00 20.00           {an[0]}")
        if altloc_at == j:
            k += 1
            p = X[j] + 0.3
            lines.append(f"ATOM  {k:5d}  CA BALA {chain}{resseq:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  0.40 20.00           C")
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


# --- parsing and mapping ----------------------------------------------------------------------------------------------

def test_mapping_recovers_uniprot_indices_with_a_numbering_offset_a_gap_and_a_tag():
    ref = ref_sequence()
    P = kp.parse_pdb(pdb_text(ref, 1001, 21, 290, skip=(50, 51, 52)))                 # author numbers offset by 980, residues missing
    m = kp.map_chain(P, "A", ref, "P99999")
    assert m["problems"] == [] and m["identity"] == 1.0
    assert all(idx == key[0] - 980 for key, idx in m["index"].items())
    assert len(m["index"]) == 270 - 3 and m["n_dbref_disagreements"] == 0


def test_mapping_refuses_a_wrong_accession_wrong_dbref_offset_and_undocumented_mutation():
    ref = ref_sequence()
    assert any("DBREF" in p for p in kp.map_chain(kp.parse_pdb(pdb_text(ref, 1, 21, 290)), "A", ref, "P00000")["problems"])
    bad_off = kp.map_chain(kp.parse_pdb(pdb_text(ref, 1, 21, 290, db_begin=25)), "A", ref, "P99999")
    assert any("disagree with the DBREF" in p for p in bad_off["problems"])
    other = ref[:100] + ("A" if ref[100] != "A" else "C") + ref[101:]
    txt = pdb_text(other, 1, 21, 290)
    m = kp.map_chain(kp.parse_pdb(txt), "A", ref, "P99999")
    assert any("without a SEQADV" in p for p in m["problems"])
    num = 101 - 21 + 1
    documented = pdb_text(other, 1, 21, 290, seqadv=[(kp.THREE2ONE and THREE[other[100]], num, "ENGINEERED MUTATION")])
    assert kp.map_chain(kp.parse_pdb(documented), "A", ref, "P99999")["problems"] == []


def test_shifted_mapping_control_collapses_identity():
    ref = ref_sequence()
    P = kp.parse_pdb(pdb_text(ref, 1, 21, 290))
    m = kp.map_chain(P, "A", ref, "P99999")
    assert kp.shifted_identity(P, "A", ref, m, 0) == 1.0
    assert kp.shifted_identity(P, "A", ref, m, 1) < 0.2


def test_altloc_policies_choose_the_documented_location():
    ref = ref_sequence()
    P = kp.parse_pdb(pdb_text(ref, 1, 21, 290, altloc_at=5), "blank_then_A")
    Q = kp.parse_pdb(pdb_text(ref, 1, 21, 290, altloc_at=5), "highest_occupancy")
    r, q = P.residues["A"][5], Q.residues["A"][5]
    assert r.altloc_ca and np.allclose(r.ca, helix(270)[5], atol=1e-3)               # 'A' (occupancy 1.0)
    assert np.allclose(q.ca, helix(270)[5], atol=1e-3)                                # 'A' also has the highest occupancy
    assert not P.residues["A"][4].altloc_ca


def test_unread_modified_residues_are_not_read_but_keep_their_place_in_the_sequence():
    ref = ref_sequence()
    txt = pdb_text(ref, 1, 21, 290)
    lines = txt.splitlines()
    for i, ln in enumerate(lines):                                                   # residue 10 becomes phosphothreonine
        if ln.startswith("ATOM") and ln[22:26].strip() == "10":
            lines[i] = ln[:17] + "TPO" + ln[20:]
    lines.insert(1, f"MODRES XXXX TPO A   10  {THREE[ref[29]]}  MODIFIED RESIDUE")
    P = kp.parse_pdb("\n".join(lines) + "\n")
    r = [x for x in P.residues["A"] if x.resseq == 10][0]
    assert not r.read and r.ca is not None
    assert kp.map_chain(P, "A", ref, "P99999")["problems"] == []


def test_correspondence_uses_residues_read_in_both_chains_by_uniprot_index():
    ref = ref_sequence()
    PA = kp.parse_pdb(pdb_text(ref, 1, 21, 290, skip=(10,)))
    PB = kp.parse_pdb(pdb_text(ref, 500, 31, 295, skip=(20,)))
    ma, mb = kp.map_chain(PA, "A", ref, "P99999"), kp.map_chain(PB, "A", ref, "P99999")
    c = kp.correspondence(PA, ma, "A", PB, mb, "A")
    idx = [i for _, _, i in c["pairs"]]
    # shared UniProt span 31..290; index 31 is the first residue of B; A lacks index 21+10=31 hmm: A skip=(10,) is index 31, B skip=(20,) is index 51
    assert idx == [i for i in range(32, 291) if i != 51]
    assert c["n_shared"] == 258 and 31 not in idx and 51 not in idx
    assert all(a[0] == i - 21 + 1 for a, _, i in c["pairs"]) and all(b[0] == i - 31 + 500 for _, b, i in c["pairs"])


# --- independent numerical references ----------------------------------------------------------------------------------

def test_kappa_matches_the_analytic_helix_value():
    X = helix(10)
    s = 2 * 2.3 ** 2 * (1 - math.cos(math.radians(100)))
    want = math.acos((s * math.cos(math.radians(100)) + 1.5 ** 2) / (s + 1.5 ** 2))
    assert kp.kappa(*X[2:5]) == pytest.approx(want, abs=1e-12)


def test_tau_textbook_dihedrals_and_chirality():
    zig = [np.array(p, float) for p in ((0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 0))]       # planar trans zig-zag: |tau| = pi
    assert abs(kp.tau(*zig)) == pytest.approx(math.pi, abs=1e-12)
    p = [np.array(v, float) for v in ((1, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, 1))]         # dihedral of +-90 degrees
    t = kp.tau(*p)
    assert abs(t) == pytest.approx(math.pi / 2, abs=1e-12)
    assert kp.tau(*[q * np.array([-1, 1, 1]) for q in p]) == pytest.approx(-t, abs=1e-12)  # mirror image flips the sign
    assert kp.kappa(np.zeros(3), np.array([1.0, 0, 0]), np.array([2.0, 0, 0])) is None     # collinear: undefined


def test_windows_are_rigid_motion_invariant_and_circular_near_pi():
    ref = ref_sequence()
    rng = np.random.default_rng(5)
    X = helix(270) + rng.normal(0, 0.4, (270, 3))
    PA = kp.parse_pdb(pdb_text(ref, 1, 21, 290, coords=X))
    Y = X + rng.normal(0, 0.3, X.shape)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    Q *= np.sign(np.linalg.det(Q))
    Z = Y @ Q.T + np.array([50.0, -20.0, 7.0])
    PB1, PB2 = kp.parse_pdb(pdb_text(ref, 1, 21, 290, coords=Y)), kp.parse_pdb(pdb_text(ref, 1, 21, 290, coords=Z))
    w1 = kp.independent_windows(PA.residues["A"], PB1.residues["A"])
    w2 = kp.independent_windows(PA.residues["A"], PB2.residues["A"])
    for k in ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4"):
        a = np.array([np.nan if v is None else v for v in w1[k]])
        b = np.array([np.nan if v is None else v for v in w2[k]])
        assert np.array_equal(np.isnan(a), np.isnan(b)) and np.nanmax(np.abs(a - b)) < 2e-3      # file coordinates carry 3 decimals
    assert w1["global_rmsd"] == pytest.approx(w2["global_rmsd"], abs=2e-3)
    assert all(v is None or 0 <= v <= math.pi + 1e-12 for v in w1["d_tau"])


def test_a_numbering_gap_or_long_link_invalidates_exactly_the_windows_that_contain_it():
    ref = ref_sequence()
    A = pdb_text(ref, 1, 21, 290)
    B = pdb_text(ref, 1, 21, 290, skip=(100,))                                         # residue 21+100 absent in B
    PA, PB = kp.parse_pdb(A), kp.parse_pdb(B)
    ma, mb = kp.map_chain(PA, "A", ref, "P99999"), kp.map_chain(PB, "A", ref, "P99999")
    c = kp.correspondence(PA, ma, "A", PB, mb, "A")
    ra = {r.key: r for r in PA.residues["A"]}
    rb = {r.key: r for r in PB.residues["A"]}
    w = kp.independent_windows([ra[a] for a, _, _ in c["pairs"]], [rb[b] for _, b, _ in c["pairs"]])
    j = 100                                                                            # first position after the gap
    assert w["links_broken"] == [j - 1]
    assert [i for i, v in enumerate(w["valid3"]) if not v and 1 <= i <= len(c["pairs"]) - 2] == [j - 1, j]
    assert [i for i, v in enumerate(w["valid4"]) if not v and 1 <= i <= len(c["pairs"]) - 3] == [j - 2, j - 1, j]


def test_trimmed_displacement_and_ties():
    rng = np.random.default_rng(2)
    XA = rng.normal(0, 10, (50, 3))
    XB = XA.copy()
    XB[:5] += 30.0                                                                     # 5 displaced positions, 10 % of 50
    d, d0 = kp.trimmed_displacement(XA, XB, 0.2)
    assert d[5:].max() < 1e-9 and d[:5].min() > 25 and d0.max() > d0.min()             # the trimmed fit ignores the displaced group
    d_same, _ = kp.trimmed_displacement(XA, XA.copy(), 0.2)
    assert d_same.max() < 1e-9


def test_classification_rule_redundant_complementary_insufficient():
    n = 120
    rng = np.random.default_rng(0)
    base = rng.random(n)
    win = {"d_kappa": list(base), "d_i_i2": list(base * 1.0), "rmsd3": list(base * 2), "d_tau": list(rng.random(n)),
           "d_i_i3": list(rng.random(n)), "rmsd4": list(rng.random(n))}
    c = kp.classify_pair(win)
    assert c["3"]["class"] == "redundant" and c["3"]["overlap_O"] == 10
    assert c["4"]["class"] == "complementary" and c["4"]["overlap_O"] <= 7
    short = {k: v[:30] for k, v in win.items()}
    assert kp.classify_pair(short)["3"]["class"] == "insufficient_evidence"
    assert kp.top_k([1.0, 5.0, 5.0, None, 2.0], 2) == [1, 2]                           # ties: lower position first


# --- selection on a mock RCSB -----------------------------------------------------------------------------------------

class MockNet(sel.Net):
    """Entries per target: ids with (resolution, ligand, pdb text). Records what was requested."""

    def __init__(self, targets):
        super().__init__()
        self.t = targets
        self.requested = []

    def search(self, accession, protocol):
        ids = [e for e, _ in sorted(self.t[accession]["entries"].items(), key=lambda kv: (kv[1]["res"], kv[0]))]
        return ids, len(ids)

    def json(self, url):
        parts = url.split("/")
        entry = parts[-1] if "entry/" in url else parts[-2]
        for acc, t in self.t.items():
            e = t["entries"].get(entry)
            if not e:
                continue
            if "/core/entry/" in url:
                return {"rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1"], "non_polymer_entity_ids": ["2"]},
                        "rcsb_entry_info": {"resolution_combined": [e["res"]]}, "rcsb_accession_info": {"initial_release_date": "2020-01-01T00:00:00"},
                        "symmetry": {"space_group_name_H_M": "P 1"}}
            if "nonpolymer_entity" in url:
                return {"pdbx_entity_nonpoly": {"comp_id": e["lig"]}, "rcsb_nonpolymer_entity": {"formula_weight": 0.4}}
            if "/polymer_entity/" in url:
                return {"rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ["A"], "reference_sequence_identifiers": [
                    {"database_name": "UniProt", "database_accession": acc}]}, "entity_poly": {"rcsb_sample_sequence_length": 270, "rcsb_mutation_count": 0}}
        return None

    def pdb_text(self, entry):
        self.requested.append(entry)
        self.files += 1
        for t in self.t.values():
            if entry in t["entries"]:
                return t["entries"][entry]["text"]            # None: no PDB-format file
        return None

    def uniprot_fasta(self, acc):
        return self.t[acc]["ref"], "sha"


def _protocol():
    return json.loads((KP_DIR / "protocol.json").read_text())


def test_selection_is_deterministic_excludes_exposed_entries_and_targets_and_records_every_reason():
    ref = ref_sequence()
    def entries(acc, specs):
        return {e: {"res": r, "lig": lig, "text": pdb_text(ref, 1, 21, 290, acc=acc) if ok else None} for e, r, lig, ok in specs}
    targets = {
        "P00001": {"ref": ref, "entries": entries("P00001", [("1AAA", 1.5, "L1", True), ("1BBB", 1.6, "L1", True), ("1CCC", 1.7, "L2", True),
                                                             ("1DDD", 1.8, "L3", True), ("1EEE", 1.9, "L4", True), ("1FFF", 2.0, "L5", True)])},
        "P00002": {"ref": ref, "entries": entries("P00002", [("2AAA", 1.5, "L1", True), ("2BBB", 1.6, "L2", True)])},          # < 4 entries
        "P00003": {"ref": ref, "entries": entries("P00003", [("3AAA", 1.4, "L1", False), ("3BBB", 1.5, "L2", False),
                                                             ("3CCC", 1.6, "L3", False), ("3DDD", 1.7, "L4", False)])},       # no usable file
    }
    protocol = _protocol()
    runs = []
    for _ in range(2):
        net = MockNet({a: {"ref": v["ref"], "entries": dict(v["entries"])} for a, v in targets.items()})
        rec = sel.select(net, protocol, [("T1", "P00001"), ("T2", "P00002"), ("T3", "P00003"), ("T4", "P00004")], {"1AAA"}, {"P00004"}, lambda s: None)
        runs.append(rec)
    assert json.dumps(runs[0], default=str) == json.dumps(runs[1], default=str)
    rec = runs[0]
    assert rec["targets"]["T4"]["status"] == "excluded" and "exposed" in rec["targets"]["T4"]["reason"]
    assert rec["targets"]["T1"]["excluded_exposed_entries"] == ["1AAA"]
    assert rec["targets"]["T2"]["status"] == "excluded" and "fewer than 4" in rec["targets"]["T2"]["reason"]
    assert rec["targets"]["T3"]["status"] == "excluded"
    assert all(x["status"] == "excluded" and "no PDB-format file" in x["reason"] for x in rec["targets"]["T3"]["examined_entries"])
    t1 = rec["targets"]["T1"]
    assert t1["status"] == "selected" and [(p["A"], p["B"]) for p in rec["panel"]] == [("1BBB", "1CCC"), ("1DDD", "1EEE")]
    assert "1AAA" not in {e["entry"] for e in t1["examined_entries"]}                   # exposed entries are never examined or downloaded
    assert [r["rank"] for r in rec["ranking"]] == [1, 2, 3] and rec["ranking"][0]["gene"] == "T1"
    assert all(p["n_shared_positions"] == 270 for p in rec["panel"])


def test_selection_stops_after_three_targets_not_three_pairs():
    ref = ref_sequence()
    def ents(acc, pre):
        return {f"{pre}{k}": {"res": 1.5 + 0.1 * k, "lig": f"L{k}", "text": pdb_text(ref, 1, 21, 290, acc=acc)} for k in range(1, 5)}
    accs = ["P00001", "P00002", "P00003", "P00004"]
    net = MockNet({a: {"ref": ref, "entries": ents(a, str(i + 1))} for i, a in enumerate(accs)})
    rec = sel.select(net, _protocol(), [(f"T{i + 1}", a) for i, a in enumerate(accs)], set(), set(), lambda s: None)
    assert len(rec["panel"]) == 6 and len({p["gene"] for p in rec["panel"]}) == 3          # 3 targets, 2 pairs each
    assert [t["status"] for t in rec["targets"].values()].count("selected") == 3
    assert rec["targets"]["T4"]["status"] == "not examined (panel complete)"


def test_pair_partner_needs_a_different_primary_ligand():
    ref = ref_sequence()
    es = {e: {"res": r, "lig": lig, "text": pdb_text(ref, 1, 21, 290, acc="P00001")} for e, r, lig in
          [("1AAA", 1.5, "L1"), ("1BBB", 1.6, "L1"), ("1CCC", 1.7, "L1"), ("1DDD", 1.8, "L1")]}
    net = MockNet({"P00001": {"ref": ref, "entries": es}})
    rec = sel.select(net, _protocol(), [("T1", "P00001")], set(), set(), lambda s: None)
    assert rec["panel"] == [] and rec["targets"]["T1"]["status"] == "excluded"


# --- structpair end to end on a SYNTHETIC pair ------------------------------------------------------------------------

def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


def test_independent_windows_agree_with_structpair_on_a_synthetic_pair_with_a_gap(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    import csv
    ref = ref_sequence()
    rng = np.random.default_rng(3)
    X = helix(270) + rng.normal(0, 0.3, (270, 3))
    Y = X + rng.normal(0, 0.5, X.shape)
    (tmp_path / "SYNTHETIC_A.pdb").write_text(pdb_text(ref, 1, 21, 290, coords=X, skip=(60, 61)))
    (tmp_path / "SYNTHETIC_B.pdb").write_text(pdb_text(ref, 1, 21, 290, coords=Y))
    PA, PB = (kp.parse_pdb((tmp_path / f"SYNTHETIC_{k}.pdb").read_text()) for k in "AB")
    ma, mb = kp.map_chain(PA, "A", ref, "P99999"), kp.map_chain(PB, "A", ref, "P99999")
    c = kp.correspondence(PA, ma, "A", PB, mb, "A")
    spec = {"name": "SYNTHETIC kin-pair check", "structures": {"A": {"label": "SYNTHETIC A", "file": "SYNTHETIC_A.pdb", "chain": "A"},
                                                         "B": {"label": "SYNTHETIC B", "file": "SYNTHETIC_B.pdb", "chain": "A"}},
            "correspondence": {"pairs": [[a[0], b[0]] for a, b, _ in c["pairs"]]}, "on_missing": "exclude", "resname_mismatch": "flag"}
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=py)
    r = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out")],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader(open(tmp_path / "out/windows.csv")))
    ra = {r_.key: r_ for r_ in PA.residues["A"]}
    rb = {r_.key: r_ for r_ in PB.residues["A"]}
    w = kp.independent_windows([ra[a] for a, _, _ in c["pairs"]], [rb[b] for _, b, _ in c["pairs"]])
    cols = {"d_kappa": "d_kappa_rad", "d_i_i2": "d_i_i2_angstrom", "rmsd3": "rmsd3_angstrom", "d_tau": "d_tau_rad",
            "d_i_i3": "d_i_i3_angstrom", "rmsd4": "rmsd4_angstrom"}
    worst = 0.0
    for k, col in cols.items():
        for j, row in enumerate(rows):
            got = None if row[col] == "" else float(row[col])
            assert (got is None) == (w[k][j] is None), (k, j)
            if got is not None:
                worst = max(worst, abs(got - w[k][j]))
    assert worst < 1e-4

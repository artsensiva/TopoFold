"""KRAS K1: provenance-check helpers (SYNTHETIC inputs), and consistency of the frozen protocol, parameter file and schematic."""

import hashlib
import io
import json
import math
import re
import struct
import sys
import zipfile
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
K1 = ROOT / "benchmarks/audit/kras_k1"
sys.path.insert(0, str(K1))
import check_k1_provenance as C  # noqa: E402
import plot_design_schematic as D  # noqa: E402

P = json.loads((K1 / "protocol_params.json").read_text())
TEXT = (K1 / "PROTOCOL.md").read_text()


def _frame(natoms=483, step=0, t=0.0, box=(7.0, 6.9, 6.8), nbytes=100):
    b = struct.pack(">iiif", 1995, natoms, step, t) + struct.pack(">9f", box[0], 0, 0, 0, box[1], 0, 0, 0, box[2])
    b += struct.pack(">i", natoms) + struct.pack(">f", 1000.0) + struct.pack(">6i", 0, 0, 0, 1, 1, 1) + struct.pack(">ii", 5, nbytes)
    return b


def test_frame_fields_reads_atoms_time_box_and_padded_size():
    f = C.frame_fields(_frame(step=500, t=1000.0, nbytes=101))
    assert (f["natoms"], f["step"], f["time_ps"]) == (483, 500, 1000.0)
    assert f["box_nm"] == pytest.approx([7.0, 6.9, 6.8])
    assert f["frame_bytes"] == 92 + 104


def test_last_frame_accepts_only_a_frame_that_ends_at_end_of_file():
    f1, f2 = _frame(t=0.0, nbytes=8) + bytes(8), _frame(step=1, t=1000.0, nbytes=12) + bytes(12)
    data = bytes(40) + f1 + f2                                   # SYNTHETIC tail window
    base = 1000
    size = base + len(data)
    got = C.last_frame(data, base, size, 483)
    assert got["time_ps"] == 1000.0 and got["offset"] == base + 40 + len(f1)
    assert C.last_frame(data + bytes(4), base, size + 4, 483) is None


S1_TEXT = """
                  G12D            GDP                5             32 741           36           28
                  G12D             GTP               10             32 743           37           28
                Wild-type          GTP               10             32 740           36           28
   Total          (170 µs)           85
"""


def test_parse_s1_rows_and_names():
    s1 = C.parse_s1(S1_TEXT)
    assert set(s1) == {"GDP_G12D", "GTP_G12D", "GTP_WT"}
    assert s1["GTP_G12D"] == {"replicas": 10, "atoms": 32743, "K+": 37, "Cl-": 28}


def test_pyinteraph_identity_uses_labels_only():
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("PyInteraph_data/GDP_G12D/salt-bridges.dat", "A-16LYS_sc.lys.NZp:A-12ASP_sc.asp.COOn\t99.9\n")
        z.writestr("PyInteraph_data/GDP_G12D/hydrophobic-clusters.dat", "A-5LYS_x:A-3GLU_y\t12\n")
        z.writestr("PyInteraph_data/GTP_WT/salt-bridges.dat", "A-5LYS_sc:A-3GLU_sc\t1.0\n")
    out = C.pyinteraph_identity(buf.getvalue())
    assert out["GDP_G12D"]["residue_12_names"] == {"ASP": 1}
    assert out["GTP_WT"]["residue_12_names"] == {}
    assert "99.9" not in json.dumps(out)


def test_committed_replica_counts_match_s1_and_zenodo():
    r = json.loads((K1 / "results/k1_provenance_checks.json").read_text())
    rc = r["replica_counts"]
    assert (rc["compared"], rc["equal"]) == (14, 14)
    assert r["s1_table"]["total_replicas"] == rc["zenodo_trajectory_files"] == 85
    assert r["s1_table"]["sha256"] == C.S1_SHA256
    g = r["g12d_summary"]
    assert g["files"] == 15 and g["natoms"] == [483] and g["interval_ps"] == [1000.0] and g["first_time_ps"] == [0.0]
    assert [k for k, v in g["last_time_ps"].items() if v == 2315000.0] == ["GTP_G12D_1.xtc"]
    ids = r["pyinteraph_label_identity"]["systems"]
    assert ids["GDP_G12D"]["residue_12_names"] == ids["GTP_G12D"]["residue_12_names"] == {"ASP": 3}
    assert not any(v["nucleotide_or_mg_labels"] for v in ids.values())


# --- the frozen decision parameters ------------------------------------------------------------------------------------

DEC = P["decision"]


def test_parameter_file_is_the_frozen_version():
    assert P["status"].startswith("FROZEN") and P["version"] == "1" and "DRAFT" not in P["status"]
    assert "not held out" in P["endpoint"]["evidence_class"] and "not functional ground truth" in P["endpoint"]["evidence_class"]
    assert P["channel_selection"]["fixed_in_sensitivity"].startswith("T* and B* stay fixed")
    assert "gate" in P["endpoint"] and "does not show that the endpoint can discriminate methods" in P["endpoint"]["gate"]["meaning"]


def test_gap_threshold_integer_form_equals_the_exact_fraction_for_every_k_and_match_count():
    d = Fraction(DEC["delta"].split()[0])
    assert d == Fraction(3, 20) == Fraction(str(DEC["delta_value"]))
    for k in range(8, 39):
        for diff in range(-k, k + 1):
            assert (20 * diff >= 3 * k) == (Fraction(diff, k) >= d)


def test_correspondence_and_jaccard_integer_forms():
    for k in range(8, 39):
        for m in range(0, k + 1):
            assert (2 * m >= k) == (Fraction(m, k) >= Fraction(1, 2))
    for inter in range(0, 40):
        for union in range(max(inter, 1), 80):
            assert (5 * inter >= 4 * union) == (Fraction(inter, union) >= Fraction(4, 5))


def test_balanced_subset_counts():
    b = DEC["balanced_subsets"]
    assert math.comb(10, 5) == b["count_when_all_retained"] == 252
    assert math.ceil(b["min_fraction"] * 252 - 1e-9) == b["min_count_when_N_252"] == 227


def test_gate_bounds_from_the_parameter_file():
    g = P["endpoint"]["gate"]
    assert (g["min_positive"], g["min_segments"], g["min_segment_length"]) == (8, 2, 3)
    assert math.floor(g["max_positive_fraction"] * 153 + 1e-9) == 38


def test_prose_states_the_corrected_rules():
    must = ["not held out", "not independent validation", "not functional ground truth",
            "It does not demonstrate that the endpoint can discriminate methods",
            "chosen once", "No winning channel is reselected inside any of them",
            "at least 0.15 in every one of", "A merely positive gap is insufficient",
            "at least 227 of 252", "*Jaccard condition* (evaluated in the primary analysis only)",
            "*Segment condition* (evaluated in the primary analysis only)", "5 |∩| < 4 |∪|",
            "lower residue number first", "separately frozen, genuinely prospective evaluation on another system",
            "not evidence of functional or clinical validity"]
    missing = [m for m in must if m not in TEXT]
    assert not missing, missing


def test_same_sign_wording_survives_only_in_the_change_log():
    sec = TEXT.split("## 8. Changes from the 2b75a0d draft")
    assert "same sign" not in sec[0] and "same sign" in sec[1]


def test_prose_parameters_and_schematic_agree():
    prov = json.loads((K1 / "results/k1_provenance_checks.json").read_text())
    feas = json.loads(D.FEAS.read_text())
    f = D.facts(P, prov, feas)
    assert (f["gdp"], f["gtp"], f["files"], f["natoms"], f["subsets"], f["mb"]) == (5, 10, 18, 483, 252, 60)
    assert (f["delta"], f["corr"], f["jac"], f["bal_pct"], f["pos_min"], f["pos_max"], f["segs"], f["seglen"]) == \
        (0.15, 0.5, 0.8, 90, 8, 38, 2, 3)
    for needle in ("0.15", "0.5", "0.8", "90 %", "252", "227", "38", "300–2000"):
        assert needle in TEXT
    src = (K1 / "plot_design_schematic.py").read_text()
    assert "Gap ≥ {f['delta']} in all of" in src and "every leave-one-out run" in src and "same sign" not in src


def test_schematic_record_matches_files_on_disk():
    rec = json.loads((D.OUT_DIR / f"{D.STEM}_provenance.json").read_text())
    for rel, h in rec["inputs_sha256"].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == h
    assert hashlib.sha256((D.OUT_DIR / f"{D.STEM}.png").read_bytes()).hexdigest() == rec["output"]["sha256"]
    assert hashlib.sha256(Path(D.__file__).read_bytes()).hexdigest() == rec["generator"]["sha256"]
    assert rec["output"]["bytes"] < 5 * 1024 * 1024


def test_no_banned_wording_in_the_protocol_files():
    low = (TEXT + json.dumps(P)).lower()
    assert not [w for w in ("authentic", "pristine", "revolutionary", "obliterate", "definitive", "zero hyperparameters") if re.search(w, low)]

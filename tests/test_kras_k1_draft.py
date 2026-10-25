"""Tests for the KRAS K1 draft helpers (SYNTHETIC inputs) and for the consistency of the committed K1 records."""

import hashlib
import io
import json
import math
import struct
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
K1 = ROOT / "benchmarks/audit/kras_k1"
sys.path.insert(0, str(K1))
import check_k1_provenance as C  # noqa: E402
import plot_design_schematic as D  # noqa: E402


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
    base = 1000                                                   # file offset of data[0] (4-byte aligned)
    size = base + len(data)
    got = C.last_frame(data, base, size, 483)
    assert got["time_ps"] == 1000.0 and got["offset"] == base + 40 + len(f1)
    assert C.last_frame(data + bytes(4), base, size + 4, 483) is None   # nothing ends exactly at EOF


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
    assert sorted(set(g["last_time_ps"].values())) == [2000000.0, 2315000.0]
    assert [k for k, v in g["last_time_ps"].items() if v == 2315000.0] == ["GTP_G12D_1.xtc"]
    ids = r["pyinteraph_label_identity"]["systems"]
    assert ids["GDP_G12D"]["residue_12_names"] == ids["GTP_G12D"]["residue_12_names"] == {"ASP": 3}
    assert not any(v["nucleotide_or_mg_labels"] for v in ids.values())


def test_params_are_a_draft_and_the_figure_numbers_follow_them():
    p = json.loads((K1 / "protocol_params_DRAFT.json").read_text())
    assert p["status"].startswith("DRAFT") and p["metric"] == "precision_at_k"
    prov = json.loads((K1 / "results/k1_provenance_checks.json").read_text())
    feas = json.loads(D.FEAS.read_text())
    f = D.facts(p, prov, feas)
    assert (f["gdp"], f["gtp"], f["files"], f["natoms"]) == (5, 10, 18, 483)
    assert f["subsets"] == math.comb(10, 5) == 252
    assert f["mb"] == round(59986475 / 1e6) == 60
    assert math.floor(p["endpoint"]["max_positive_fraction"] * 153) == 38


def test_schematic_record_matches_files_on_disk():
    rec = json.loads((D.OUT_DIR / f"{D.STEM}_provenance.json").read_text())
    for rel, h in rec["inputs_sha256"].items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == h
    assert hashlib.sha256((D.OUT_DIR / f"{D.STEM}.png").read_bytes()).hexdigest() == rec["output"]["sha256"]
    assert hashlib.sha256(Path(D.__file__).read_bytes()).hexdigest() == rec["generator"]["sha256"]
    assert rec["output"]["bytes"] < 5 * 1024 * 1024

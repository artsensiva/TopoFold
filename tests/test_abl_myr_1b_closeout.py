"""Tests of the ABL-MYR-1b closeout inventory helpers on SYNTHETIC files (no Dryad data, no trajectory)."""

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/audit/abl_myr_1b"))
import closeout_inventory as C  # noqa: E402


def _h(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def test_sha256_file_matches_hashlib(tmp_path):
    p = tmp_path / "SYNTHETIC.bin"
    data = bytes(range(256)) * 9000          # larger than one read chunk
    p.write_bytes(data)
    assert C.sha256_file(p, chunk=1000) == _h(data)


def test_hash_pairs_finds_both_record_shapes(tmp_path):
    (tmp_path / "a").mkdir()
    rec = {"code": {"path": "a/run.py", "sha256": "0" * 64},
           "deps": {"a/lib.py": "1" * 64},
           "spec_sha256": "2" * 64,                          # no path: not a (path, hash) pair
           "lineage": {"original specification": "3" * 64},  # a name, not a path
           "nested": [{"a/x.json": "4" * 64}]}
    assert sorted(C.hash_pairs(rec, tmp_path)) == [("a/lib.py", "1" * 64), ("a/run.py", "0" * 64), ("a/x.json", "4" * 64)]


def test_check_freezes_reports_mismatch_and_missing(tmp_path):
    (tmp_path / "d").mkdir()
    (tmp_path / "d/ok.py").write_bytes(b"SYNTHETIC ok")
    (tmp_path / "d/changed.py").write_bytes(b"SYNTHETIC changed")
    fz = {"x": {"d/ok.py": _h(b"SYNTHETIC ok"), "d/changed.py": _h(b"SYNTHETIC original"), "d/gone.py": "5" * 64}}
    (tmp_path / "d/freeze.json").write_text(json.dumps(fz))
    out = C.check_freezes(["d/freeze.json"], tmp_path)
    assert (out["pairs_checked"], out["pairs_matching"]) == (3, 1)
    assert sorted(m["path"] for m in out["mismatches"]) == ["d/changed.py", "d/gone.py"]
    assert [m["now"] for m in out["mismatches"] if m["path"] == "d/gone.py"] == [None]


def test_members_against_g1(tmp_path):
    ext = tmp_path / "extracted"
    (ext / "A/A").mkdir(parents=True)
    (ext / "A/A/good.pdb").write_bytes(b"SYNTHETIC good")
    (ext / "A/A/bad.pdb").write_bytes(b"SYNTHETIC tampered")
    res = {"gates": {"G1": {"files": {
        "A.zip:good": {"member": "A/good.pdb", "bytes": 14, "sha256": _h(b"SYNTHETIC good")},
        "A.zip:bad": {"member": "A/bad.pdb", "bytes": 14, "sha256": _h(b"SYNTHETIC bad")},
        "A.zip:missing": {"member": "A/missing.pdb", "bytes": 1, "sha256": "6" * 64}}}}}
    old = C.ROOT
    C.ROOT = tmp_path
    try:
        m = C.members(res, ext)
    finally:
        C.ROOT = old
    assert (m["recorded_in_g1"], m["present"], m["matching"]) == (3, 2, 1)
    assert m["missing"] == ["A.zip:missing"] and m["mismatching"] == ["A.zip:bad"]


def test_tree_manifest_is_order_stable_and_empty_safe(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "ROOT", tmp_path)
    d = tmp_path / "w"
    (d / "sub").mkdir(parents=True)
    (d / "b.npz").write_bytes(b"SYNTHETIC b")
    (d / "sub/a.json").write_bytes(b"SYNTHETIC a")
    t1 = C.tree(d)
    assert (t1["files"], t1["bytes"]) == (2, 22)
    assert C.tree(d)["manifest_sha256"] == t1["manifest_sha256"]
    (d / "b.npz").write_bytes(b"SYNTHETIC B")
    assert C.tree(d)["manifest_sha256"] != t1["manifest_sha256"]
    empty = C.tree(tmp_path / "absent")
    assert empty == {"present": False, "files": 0, "bytes": 0, "manifest_sha256": None, "entries": []}


def test_sums_text_lists_committed_files_and_used_archives_only():
    inv = {"committed": {"entries": [{"path": "x/a.json", "sha256": "7" * 64}]},
           "archives": [{"name": "PCA.zip", "sha256_source": "8" * 64, "used_by_abl_myr_1b_comparison": True},
                        {"name": "FDA.zip", "sha256_source": "9" * 64, "used_by_abl_myr_1b_comparison": False}]}
    assert C.sums_text(inv) == (f"{'7' * 64}  x/a.json\n"
                                f"{'8' * 64}  benchmarks/audit/abl_myr_1/data_dryad/PCA.zip\n")


@pytest.mark.parametrize("name", C.USED_ARCHIVES)
def test_used_archives_are_in_the_intake_inventory(name):
    intake = json.loads(C.INTAKE.read_text())
    assert name in {f["path"] for f in intake["raw_inventory"]["files"]}


def test_reconcile_extracted_counts_duplicates_and_reader_caches():
    mem = {"files": [{"key": "A.zip:gro:Set2", "local": "x/system_em.gro"},
                     {"key": "A.zip:gro_em:Set2", "local": "x/system_em.gro"},
                     {"key": "P.zip:xtc:1", "local": "x/traj_1.xtc"}]}
    tree = {"entries": [{"path": "x/system_em.gro", "bytes": 5}, {"path": "x/traj_1.xtc", "bytes": 7},
                        {"path": "x/.traj_1.xtc_offsets.npz", "bytes": 3}, {"path": "x/.traj_1.xtc_offsets.lock", "bytes": 0},
                        {"path": "x/SYNTHETIC_stray.txt", "bytes": 2}]}
    r = C.reconcile_extracted(mem, tree)
    assert (r["g1_entries"], r["unique_member_files"]) == (3, 2)
    assert r["files_listed_under_several_g1_keys"] == {"x/system_em.gro": ["A.zip:gro:Set2", "A.zip:gro_em:Set2"]}
    assert r["files_not_in_g1"] == {"MDAnalysis XTC offset cache": {"files": 2, "bytes": 3}, "other": {"files": 1, "bytes": 2}}

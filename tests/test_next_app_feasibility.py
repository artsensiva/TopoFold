"""Tests of the feasibility metadata helpers on SYNTHETIC inputs (no network)."""

import struct
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/audit/next_app_feasibility"))
import collect_metadata as M  # noqa: E402


def _xtc_header(natoms=483, step=500, time_ps=1000.0, nbytes=1866):
    box = struct.pack(">9f", *([5.0] * 9))
    head = struct.pack(">iiif", 1995, natoms, step, time_ps) + box + struct.pack(">i", natoms)
    rest = struct.pack(">f", 1000.0) + struct.pack(">6i", 0, 0, 0, 1, 1, 1) + struct.pack(">ii", 5, nbytes)
    return head + rest


def test_xtc_frame_header_reads_atoms_time_and_padded_frame_size():
    h = M.xtc_frame_header(_xtc_header(nbytes=1866))
    assert (h["natoms"], h["step"], h["time_ps"]) == (483, 500, 1000.0)
    assert h["frame_bytes"] == 92 + 1868          # compressed bytes padded to a multiple of 4


def test_xtc_frame_header_rejects_other_data():
    bad = bytearray(_xtc_header())
    bad[0:4] = struct.pack(">i", 1234)
    with pytest.raises(ValueError):
        M.xtc_frame_header(bytes(bad))


def _pdb_line(serial, name, resname, resid, chain="A"):
    return f"ATOM  {serial:5d} {name:<4} {resname:>3} {chain}{resid:4d}    {0:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00"


def test_pdb_sequence_one_entry_per_residue_in_file_order():
    lines = [_pdb_line(1, "N", "MET", 1), _pdb_line(2, "CA", "MET", 1), _pdb_line(3, "N", "HIE", 2),
             _pdb_line(4, "CA", "HIE", 2), _pdb_line(5, "N", "XYZ", 3), "TER"]
    assert M.pdb_sequence("\n".join(lines)) == [(1, "M"), (2, "H"), (3, "X")]


def test_map_to_reference_lists_every_mismatch_at_offset_zero():
    ref = "MTEYKLVVVG"                                     # SYNTHETIC reference
    res = [(4, "Y"), (5, "K"), (6, "A"), (11, "G")]       # residue 6 differs; 11 is beyond the reference
    out = M.map_to_reference(res, ref)
    assert (out["n"], out["first"], out["last"]) == (4, 4, 11)
    assert out["mismatches"] == [{"resid": 6, "topology": "A", "reference": "L"},
                                 {"resid": 11, "topology": "G", "reference": None}]


def test_fasta_seq_joins_sequence_lines():
    assert M.fasta_seq(">sp|SYNTHETIC|X\nMTEY\nKLVV\n") == "MTEYKLVV"

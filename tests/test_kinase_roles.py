"""Tests for topofold_inputs.kinase_roles (no MDAnalysis needed).

SYNTHETIC sequences for missing residues, insertions, mutations, insertion codes, ambiguous motifs and ambiguous
alignments; the committed KIN-ROB-1 panel fixture (public PDB/UniProt data, tests/fixtures/
kin_rob_1_panel_observed_chains.json) for the frozen panel, including the original PKA 1ATP specification error.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.atom_adapter import KINASE_DESCRIPTOR_ROLES  # noqa: E402
from topofold_inputs.ca_adapter import AdapterError  # noqa: E402
from topofold_inputs.kinase_roles import (DERIVED_ROLES, STANDARD, MotifWindow, Numbering, ObservedResidue,  # noqa: E402
                                          Reference, RoleMappingError, map_kinase_roles)

THREE = {v: k for k, v in STANDARD.items()}
POS = {"beta3_Lys": 20, "alphaC_Glu": 38, "HRD_H": 80, "DFG_Asp": 100}   # SYNTHETIC reference indices
AUTHOR = 500                                                              # SYNTHETIC author = reference index + 500


def synthetic_reference(extra=None, length=130, seed=1):
    """SYNTHETIC kinase-like reference: K20, E38, H80-R81-D82, D100-F101-G102; filler residues contain none of
    D, E, F, G, H, K, R, Y, so every motif occurs once unless ``extra`` adds one."""
    seq = list(np.random.default_rng(seed).choice(list("ACILMNPQSTVW"), length))
    for p, aa in {20: "K", 38: "E", 80: "H", 81: "R", 82: "D", 100: "D", 101: "F", 102: "G", **(extra or {})}.items():
        seq[p - 1] = aa
    return Reference("SYNTHETIC-KIN", "SYNTHETIC-KIN-1", "".join(seq), "SYNTHETIC")


def chain_of(ref, *, drop=(), bonded_at_gaps=False, mutate=None, insert_after=None, icode_for=None):
    """SYNTHETIC observed chain: every reference residue except ``drop`` (reference indices), author = index + 500."""
    out, prev = [], None
    for i, aa in enumerate(ref.sequence, start=1):
        if i in drop:
            continue
        bonded = None if prev is None else (prev == i - 1 or bonded_at_gaps)
        name = THREE[(mutate or {}).get(i, aa)]
        icode = ""
        resid = i + AUTHOR
        if icode_for and i in icode_for:
            resid, icode = icode_for[i]
        out.append(ObservedResidue(resid, icode, name, bonded))
        if insert_after and i == insert_after:
            out += [ObservedResidue(resid, ic, "GLY", True) for ic in "ABC"]
        prev = i
    return out


def mapped(m):
    return {r: (None if s is None else (s.resid, s.icode)) for r, s in m.roles.items()}


def expected(shift=AUTHOR):
    return {r: (POS[a] + k + shift, "") for r, (a, k) in DERIVED_ROLES.items()}


def test_role_tables_match_descriptor_roles():
    assert set(DERIVED_ROLES) == set(KINASE_DESCRIPTOR_ROLES)


def test_correct_reference_numbering_is_accepted_without_using_author_numbers():
    ref = synthetic_reference()
    m = map_kinase_roles(chain_of(ref), "A", ref, Numbering.reference_index(), POS)
    assert mapped(m) == expected()
    aln = m.report["alignment"]
    assert aln["identity"] == 1.0 and aln["mutations"] == [] and aln["insertions"] == []
    assert aln["author_numbering_diagnostic"]["reference_index_minus_author_resid"] == {"-500": 130}
    assert all(v["status"] == "observed" for v in m.report["roles"].values())


def test_literature_numbers_labelled_as_reference_are_refused_and_offset_is_not_applied():
    ref = synthetic_reference()
    lit = {a: p - 1 for a, p in POS.items()}
    with pytest.raises(RoleMappingError) as e:
        map_kinase_roles(chain_of(ref), "A", ref, Numbering.reference_index(), lit)
    assert any("beta3_Lys" in p and "expected K" in p for p in e.value.problems)
    assert e.value.report["diagnostic_offsets_not_applied"]["d"] == [1]
    assert e.value.report["accepted"] is False
    # the same labels with the offset stated explicitly are accepted
    m = map_kinase_roles(chain_of(ref), "A", ref, Numbering("reference", 1, "SYNTHETIC literature numbering"), lit)
    assert mapped(m) == expected()


def test_numbering_must_be_explicit():
    with pytest.raises(AdapterError, match="explicitly stated"):
        Numbering("reference")
    with pytest.raises(AdapterError, match="explicitly stated"):
        Numbering("reference", True)
    with pytest.raises(AdapterError, match="unknown numbering scheme"):
        Numbering("uniprot", 0)
    with pytest.raises(AdapterError, match="no reference_minus_label"):
        Numbering("structure_author", 0)
    with pytest.raises(AdapterError, match="must both be stated"):
        Reference("P00000", "", "ACDK")
    ref = synthetic_reference()
    with pytest.raises(AdapterError, match="anchors must be exactly"):
        map_kinase_roles(chain_of(ref), "A", ref, Numbering.reference_index(), {**POS, "DFG_Phe": 101})


def test_fasta_identifier_must_match_stated_accession_or_isoform():
    text = ">sp|P99999|SYNTH_HUMAN SYNTHETIC\nACDEFGHIK\n"
    assert Reference.from_fasta(text, "P99999", "P99999-1").sequence == "ACDEFGHIK"
    assert Reference.from_fasta(text.replace("P99999|", "P99999-2|"), "P99999", "P99999-2").isoform == "P99999-2"
    with pytest.raises(AdapterError, match="neither the stated accession"):
        Reference.from_fasta(text, "Q11111", "Q11111-1")


def test_missing_residues_at_a_chain_break_are_unobserved_roles():
    ref = synthetic_reference()
    m = map_kinase_roles(chain_of(ref, drop=range(102, 108)), "A", ref, Numbering.reference_index(), POS)
    got = mapped(m)
    assert got["DFG_Gly"] is None and got["DFGxxX"] is None
    assert got["DFG_Phe"] == (601, "")
    assert m.report["roles"]["DFG_Gly"]["status"] == "unobserved"
    assert "chain break" in m.report["roles"]["DFG_Gly"]["reason"]
    assert m.report["alignment"]["unobserved_reference_ranges"] == [[102, 107]]
    assert set(m.residue_specs()) == set(DERIVED_ROLES) - {"DFG_Gly", "DFGxxX"}
    assert {s.resid for s, _ in m.atom_request()} == {v[0] for k, v in got.items() if v}


def test_residues_absent_between_bonded_neighbours_are_construct_deletions():
    ref = synthetic_reference()
    with pytest.raises(RoleMappingError, match="deleted in the construct"):   # the deletion removes DFG-Gly
        map_kinase_roles(chain_of(ref, drop=range(102, 108), bonded_at_gaps=True), "A", ref,
                         Numbering.reference_index(), POS)
    m = map_kinase_roles(chain_of(ref, drop=range(55, 61), bonded_at_gaps=True), "A", ref,
                         Numbering.reference_index(), POS)       # a deletion away from the roles is recorded
    assert mapped(m) == expected()
    aln = m.report["alignment"]
    assert aln["absent_between_bonded_residues"] == [[55, 60]] and aln["unobserved_reference_ranges"] == []


def test_insertions_are_recorded_and_do_not_shift_roles():
    ref = synthetic_reference()
    obs = chain_of(ref, insert_after=60)
    m = map_kinase_roles(obs, "A", ref, Numbering.reference_index(), POS)
    assert mapped(m) == expected()
    assert [(r["resid"], r["icode"]) for r in m.report["alignment"]["insertions"]] == [(560, c) for c in "ABC"]


def test_role_mutation_is_refused_unless_allowed_and_other_mutations_are_recorded():
    ref = synthetic_reference()
    with pytest.raises(RoleMappingError, match="role or motif mutation"):
        map_kinase_roles(chain_of(ref, mutate={20: "M"}), "A", ref, Numbering.reference_index(), POS)
    with pytest.raises(RoleMappingError, match="role or motif mutation"):   # HRD Arg is a motif residue, not a role
        map_kinase_roles(chain_of(ref, mutate={81: "A"}), "A", ref, Numbering.reference_index(), POS)
    m = map_kinase_roles(chain_of(ref, mutate={20: "M"}), "A", ref, Numbering.reference_index(), POS,
                         allow_role_mutations=True)
    assert m.report["roles"]["beta3_Lys"]["status"] == "observed_mutated" and mapped(m) == expected()
    m = map_kinase_roles(chain_of(ref, mutate={60: "A"}), "A", ref, Numbering.reference_index(), POS)
    assert [(x["reference_index"], x["observed_aa"]) for x in m.report["alignment"]["mutations"]] == [(60, "A")]


def test_insertion_codes_on_role_residues_are_carried_through():
    ref = synthetic_reference()
    obs = chain_of(ref, icode_for={100: (599, "A"), 101: (599, "B")})
    m = map_kinase_roles(obs, "A", ref, Numbering.reference_index(), POS)
    assert mapped(m)["DFG_Asp"] == (599, "A") and mapped(m)["DFG_Phe"] == (599, "B")
    assert len(m.report["alignment"]["insertion_code_residues"]) == 2


def test_ambiguous_motif_window_is_refused_and_unique_window_accepted():
    ref = synthetic_reference(extra={110: "D", 111: "F", 112: "G"})
    win = {**POS, "DFG_Asp": MotifWindow(95, 115)}
    with pytest.raises(RoleMappingError, match="found 2 times"):
        map_kinase_roles(chain_of(ref), "A", ref, Numbering.reference_index(), win)
    m = map_kinase_roles(chain_of(ref), "A", ref, Numbering.reference_index(), {**POS, "DFG_Asp": MotifWindow(95, 105)})
    assert mapped(m)["DFG_Asp"] == (600, "")
    with pytest.raises(RoleMappingError, match="found 0 times"):
        map_kinase_roles(chain_of(ref), "A", ref, Numbering.reference_index(), {**POS, "DFG_Asp": MotifWindow(1, 30)})


def test_ambiguous_alignment_of_a_role_is_refused_and_bond_information_resolves_it():
    ref = synthetic_reference(extra={103: "G"})           # ...D F G G..., one G unobserved
    obs = chain_of(ref, drop={103})
    unknown = [ObservedResidue(o.resid, o.icode, o.resname, None) for o in obs]   # no bonding information
    with pytest.raises(RoleMappingError, match="ambiguous"):
        map_kinase_roles(unknown, "A", ref, Numbering.reference_index(), POS)
    m = map_kinase_roles(obs, "A", ref, Numbering.reference_index(), POS)   # G102 bonded to F101; break before 104
    assert mapped(m)["DFG_Gly"] == (602, "")
    assert m.report["ambiguity_checks"]["102"]["unique"] is True


def test_anchor_order_contradiction_is_refused():
    ref = synthetic_reference(extra={50: "D", 51: "F", 52: "G"})
    with pytest.raises(RoleMappingError, match="anchor order"):
        map_kinase_roles(chain_of(ref), "A", ref, Numbering.reference_index(), {**POS, "DFG_Asp": 50})


def test_wrong_reference_is_refused_on_identity():
    ref = synthetic_reference()
    other = synthetic_reference(seed=99)
    with pytest.raises(RoleMappingError, match="reference or isoform may be wrong"):
        map_kinase_roles(chain_of(other), "A", ref, Numbering.reference_index(), POS)


def test_structure_author_numbering():
    ref = synthetic_reference()
    obs = chain_of(ref)
    m = map_kinase_roles(obs, "A", ref, Numbering("structure_author"), {a: p + AUTHOR for a, p in POS.items()})
    assert mapped(m) == expected()
    with pytest.raises(RoleMappingError, match="expected K"):
        map_kinase_roles(obs, "A", ref, Numbering("structure_author"), {**{a: p + AUTHOR for a, p in POS.items()},
                                                                        "beta3_Lys": 519})
    with pytest.raises(RoleMappingError, match="not in chain"):
        map_kinase_roles(obs, "A", ref, Numbering("structure_author"), {**{a: p + AUTHOR for a, p in POS.items()},
                                                                        "beta3_Lys": 9999})


# --- frozen KIN-ROB-1 panel (committed public-data fixture) -------------------------------------------------------------

FX = json.loads((ROOT / "tests/fixtures/kin_rob_1_panel_observed_chains.json").read_text())
ROB = json.loads((ROOT / "benchmarks/audit/kin_rob_1/run_spec.json").read_text())["role_positions_uniprot"]
ROB_A = json.loads((ROOT / "benchmarks/audit/kin_rob_1a/run_spec.json").read_text())["role_positions_uniprot"]
ROB_RES = {(c["pdb"], c["chain"]): c for c in
           json.loads((ROOT / "benchmarks/audit/kin_rob_1/results/KIN-ROB-1_results.json").read_text())["chains"]}
ROB_A_RES = json.loads((ROOT / "benchmarks/audit/kin_rob_1a/results/KIN-ROB-1_results.json").read_text())["chains"][0]


def panel(pdb, chain):
    c = next(x for x in FX["chains"] if (x["pdb"], x["chain"]) == (pdb, chain))
    u = FX["uniprot"][c["uniprot"]]
    ref = Reference.from_fasta(u["header"] + "\n" + u["sequence"], c["uniprot"], c["uniprot"] + "-1", "UniProt")
    return [ObservedResidue(*r) for r in c["residues"]], ref


def anchors_of(positions):
    return {"beta3_Lys": positions["beta3_Lys"], "alphaC_Glu": positions["alphaC_Glu"], "HRD_H": positions["XHRD"] + 1,
            "DFG_Asp": positions["DFG_Asp"]}


def recorded_roles(rec):
    """Role -> author residue used by the frozen run, None where it reported the residue as not found."""
    c3 = rec["C3_role_residue_identity_in_structure"]["roles"]
    gone = {int(s.split(":")[0]) for s in rec["C4_required_atom_coverage"]["missing"] if "residue_not_found" in s}
    return {r: (None if v["pdb_resid"] in gone else (v["pdb_resid"], "")) for r, v in c3.items()}


@pytest.mark.parametrize("pdb,chain", [(c["pdb"], c["chain"]) for c in FX["chains"] if c["pdb"] != "1ATP"])
def test_panel_frozen_spec_accepted_and_matches_kin_rob_1(pdb, chain):
    obs, ref = panel(pdb, chain)
    m = map_kinase_roles(obs, chain, ref, Numbering.reference_index("UniProt numbering of the canonical isoform"),
                         anchors_of(ROB[ref.accession]))
    assert {r: v["reference_index"] for r, v in m.report["roles"].items()} == ROB[ref.accession]
    assert mapped(m) == recorded_roles(ROB_RES[(pdb, chain)])


def test_panel_original_pka_spec_is_refused():
    obs, ref = panel("1ATP", "E")
    with pytest.raises(RoleMappingError) as e:
        map_kinase_roles(obs, "E", ref, Numbering.reference_index("as labelled in the frozen KIN-ROB-1 spec"),
                         anchors_of(ROB["P05132"]))
    assert len([p for p in e.value.problems if "expected" in p]) == 4
    assert e.value.report["diagnostic_offsets_not_applied"]["d"] == [1]


@pytest.mark.parametrize("positions,numbering", [
    (anchors_of(ROB_A["P05132"]), Numbering.reference_index("UniProt P05132 numbering (KIN-ROB-1a)")),
    (anchors_of(ROB["P05132"]), Numbering("reference", 1, "PKA literature numbering, without the initiator Met")),
    ({"beta3_Lys": 72, "alphaC_Glu": 91, "HRD_H": 164, "DFG_Asp": 184}, Numbering("structure_author")),
])
def test_panel_correctly_declared_pka_mapping_is_accepted(positions, numbering):
    obs, ref = panel("1ATP", "E")
    m = map_kinase_roles(obs, "E", ref, numbering, positions)
    assert {r: v["reference_index"] for r, v in m.report["roles"].items()} == ROB_A["P05132"]
    assert mapped(m) == recorded_roles(ROB_A_RES)
    assert m.report["alignment"]["author_numbering_diagnostic"]["reference_index_minus_author_resid"] == {"1": 336}

"""Kinase residue-role mapping from a declared reference sequence to one structure chain (engineering capability).

Resolves the residues the kinase descriptors need (:data:`topofold_inputs.atom_adapter.KINASE_DESCRIPTOR_ROLES`) and
verifies them. No state label, threshold or classification is involved. Nothing about numbering is assumed:

* **Declared reference.** The caller states the accession, the isoform and the sequence (:class:`Reference`). The
  sequence defines reference indices (1-based); for a UniProt isoform these are UniProt numbers of that isoform.
* **Declared numbering.** Anchor labels are interpreted only as :class:`Numbering` says: ``"reference"`` (reference index =
  label + ``reference_minus_label``, an integer the caller must state, 0 included) or ``"structure_author"`` (author
  residue number and insertion code in the mapped chain). No offset is inferred or applied silently. When the motif
  checks fail, offsets under which they would pass are *reported* as a diagnostic and never used.
* **Anchors and derived roles.** Four anchors (:data:`ANCHOR_MOTIFS`): β3 Lys (``K``), αC Glu (``E``), the first residue
  of the HRD/YRD motif (``H`` or ``Y``, then ``R``, ``D``) and the DFG Asp (``D``, ``F``, ``G``). An anchor is a position
  or a :class:`MotifWindow` searched for exactly one motif occurrence. The ten descriptor roles are fixed sequence
  neighbours of the anchors in the reference (:data:`DERIVED_ROLES`).
* **Alignment, not numbering.** The observed chain is aligned to the reference (local alignment, integer scores). Author
  residue numbers are never used to place residues; their offsets to the reference are reported for inspection only.
  Physics: residues that are present in the construct but not modelled leave a chain break (peptide C–N is about
  1.33 Å, while one missing residue leaves its flanking C and N at least about 2.4 Å apart; the bond test is
  C(i−1)–N(i) ≤ :data:`PEPTIDE_BOND_MAX`). Skipping reference residues at a break, or where bonding is unknown, costs
  little and does not depend on length. Skipping them between peptide-bonded residues means the construct itself lacks
  them (an engineered deletion); it is scored like a sequence indel and reported separately.
* **Ambiguity.** For each role and motif position the alignment is recomputed with its mapping forbidden (or, for an
  unobserved position, with a match required). An equal score means another optimal alignment places it differently,
  and the mapping is refused.
* **Refusals** raise :class:`RoleMappingError` with every problem stated and the full report attached: motif or order
  contradictions, window searches with zero or several hits, low identity or too few aligned residues (wrong reference
  or isoform), ambiguous placement, role or motif residues that differ from the reference in the structure (unless
  ``allow_role_mutations``), role or motif residues absent between peptide-bonded residues, anchors outside the
  reference, or author residues that are absent or unaligned.

Recorded, not refused: mutations elsewhere, internal insertions, unaligned terminal residues (tags), unobserved
reference residues, chain breaks, insertion codes and modified or force-field residue names.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence, Tuple, Union

import numpy as np

from .atom_adapter import KINASE_DESCRIPTOR_ROLES, role_request
from .ca_adapter import AdapterError, ResidueSpec, _chain_of, _chain_source

KINASE_ROLES_VERSION = "0.1.0"
PEPTIDE_BOND_MAX = 2.0  # Å, C(i-1)–N(i)

# anchor -> ((offset from anchor, allowed one-letter residues), ...) checked in the reference
ANCHOR_MOTIFS: dict[str, tuple[tuple[int, str], ...]] = {
    "beta3_Lys": ((0, "K"),),
    "alphaC_Glu": ((0, "E"),),
    "HRD_H": ((0, "HY"), (1, "R"), (2, "D")),
    "DFG_Asp": ((0, "D"), (1, "F"), (2, "G")),
}
ANCHOR_ORDER = ("beta3_Lys", "alphaC_Glu", "HRD_H", "DFG_Asp")  # N- to C-terminal order in the kinase fold

# descriptor role -> (anchor, offset in the reference)
DERIVED_ROLES: dict[str, tuple[str, int]] = {
    "beta3_Lys": ("beta3_Lys", 0),
    "alphaC_Glu": ("alphaC_Glu", 0),
    "alphaC_Glu_plus4": ("alphaC_Glu", 4),
    "XHRD": ("HRD_H", -1),
    "XDFG_minus1": ("DFG_Asp", -2),
    "XDFG": ("DFG_Asp", -1),
    "DFG_Asp": ("DFG_Asp", 0),
    "DFG_Phe": ("DFG_Asp", 1),
    "DFG_Gly": ("DFG_Asp", 2),
    "DFGxxX": ("DFG_Asp", 5),   # sixth residue counting DFG-Asp as 1
}
assert set(DERIVED_ROLES) == set(KINASE_DESCRIPTOR_ROLES)

STANDARD = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
            "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
            "TYR": "Y", "VAL": "V"}
# modified residues and force-field protonation/variant names -> parent residue; each use is recorded
VARIANT_NAMES = {"MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C", "CME": "C", "CSD": "C", "KCX": "K",
                 "MLY": "K", "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "CYX": "C",
                 "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K"}
ONE_LETTER = {**STANDARD, **VARIANT_NAMES}
REFERENCE_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYU")

# integer alignment scores; see the module docstring for the deletion rule
MATCH, MISMATCH = 5, -3
INS_OPEN, INS_EXT = -8, -2      # observed residue absent from the reference (counted inside the alignment only)
DEL_OPEN, DEL_EXT = -1, 0       # reference residues unobserved at a chain break: any length is equally plausible
BONDED_DEL_OPEN, BONDED_DEL_EXT = -8, -2   # reference residues absent between peptide-bonded residues (construct)
MIN_IDENTITY, MIN_ALIGNED = 0.90, 30
_NEG = -(10 ** 12)


class RoleMappingError(AdapterError):
    """A refused mapping. ``problems`` lists every reason; ``report`` is the full mapping report."""

    def __init__(self, problems: Sequence[str], report: Mapping[str, Any]):
        self.problems = list(problems)
        self.report = dict(report)
        super().__init__("kinase role mapping refused: " + "; ".join(self.problems))


@dataclass(frozen=True)
class Reference:
    """A stated reference sequence. Reference index i (1-based) is ``sequence[i - 1]``."""

    accession: str
    isoform: str
    sequence: str
    source: str = ""

    def __post_init__(self):
        if not str(self.accession).strip() or not str(self.isoform).strip():
            raise AdapterError("reference accession and isoform must both be stated")
        bad = sorted(set(self.sequence) - REFERENCE_ALPHABET)
        if not self.sequence or bad:
            raise AdapterError(f"reference sequence must be non-empty one-letter amino acids; invalid {bad}")

    @classmethod
    def from_fasta(cls, text: str, accession: str, isoform: str, source: str = "") -> "Reference":
        """One-record FASTA. A UniProt header's identifier must equal ``accession`` or ``isoform``."""
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        if not lines or not lines[0].startswith(">") or sum(ln.startswith(">") for ln in lines) != 1:
            raise AdapterError("expected exactly one FASTA record")
        parts = lines[0][1:].split("|")
        if len(parts) >= 2 and parts[1] not in (accession, isoform):
            raise AdapterError(f"FASTA identifier {parts[1]!r} is neither the stated accession {accession!r} nor "
                               f"isoform {isoform!r}")
        return cls(accession, isoform, "".join(lines[1:]), source)

    def describe(self) -> dict[str, Any]:
        return {"accession": self.accession, "isoform": self.isoform, "length": len(self.sequence),
                "sequence_sha256": hashlib.sha256(self.sequence.encode()).hexdigest(), "source": self.source}


@dataclass(frozen=True)
class Numbering:
    """Numbering convention of the anchor labels (see module docstring). ``reference_minus_label`` must be stated for
    ``scheme="reference"`` and absent for ``"structure_author"``."""

    scheme: str
    reference_minus_label: int | None = None
    description: str = ""

    def __post_init__(self):
        if self.scheme == "reference":
            if isinstance(self.reference_minus_label, bool) or not isinstance(self.reference_minus_label, int):
                raise AdapterError("scheme 'reference' needs an explicitly stated integer reference_minus_label "
                                   "(0 when labels are reference indices); it is never inferred")
        elif self.scheme == "structure_author":
            if self.reference_minus_label is not None:
                raise AdapterError("scheme 'structure_author' takes no reference_minus_label")
        else:
            raise AdapterError(f"unknown numbering scheme {self.scheme!r}; expected 'reference' or 'structure_author'")

    @classmethod
    def reference_index(cls, description: str = "labels are 1-based indices of the stated reference sequence"):
        return cls("reference", 0, description)


@dataclass(frozen=True)
class MotifWindow:
    """Search an anchor's motif between two labels (inclusive, in the declared reference numbering)."""

    first: int
    last: int


@dataclass(frozen=True)
class ObservedResidue:
    """One residue of the structure chain, in chain order. ``bonded_to_previous`` is None when unknown."""

    resid: int
    icode: str
    resname: str
    bonded_to_previous: bool | None = None


AnchorLabel = Union[int, Tuple[int, str], MotifWindow]


@dataclass
class RoleMapping:
    """An accepted mapping: role -> structure residue (None when the reference residue is unobserved)."""

    chain: str
    roles: dict[str, ResidueSpec | None]
    report: dict[str, Any] = field(default_factory=dict)

    def residue_specs(self) -> dict[str, ResidueSpec]:
        """Observed roles only (for :func:`topofold_inputs.kinase_geometry.kinase_descriptors`; unmapped roles give
        NaN with reason ``role_not_mapped``)."""
        return {r: s for r, s in self.roles.items() if s is not None}

    def atom_request(self) -> list[tuple[ResidueSpec, list[str]]]:
        """Atom request for :func:`topofold_inputs.atom_adapter.extract_atoms` covering the observed roles."""
        specs = self.residue_specs()
        return role_request(specs, {r: KINASE_DESCRIPTOR_ROLES[r] for r in specs})


# --- structure side (MDAnalysis universe) ------------------------------------------------------------------------------

def observed_chain(universe, chain: str, *, selection: str | None = None, chain_attribute: str = "auto",
                   peptide_bond_max: float = PEPTIDE_BOND_MAX) -> list[ObservedResidue]:
    """Residues of one chain in topology order with peptide-bond flags from the first frame's coordinates.

    The default selection is ``protein`` plus the residue names in :data:`VARIANT_NAMES` (MDAnalysis' ``protein``
    omits e.g. SEP and TPO). Bonding uses the shortest C(i−1)–N(i) distance over alternate locations; a missing C or N
    gives ``None``."""
    sel_text = selection or "protein or resname " + " ".join(sorted(VARIANT_NAMES))
    sel = universe.select_atoms(sel_text)
    if sel.n_atoms == 0:
        raise AdapterError(f"selection {sel_text!r} matched no atoms")
    src = _chain_source(universe, sel, chain_attribute)
    universe.trajectory[0]
    out, seen, prev_c = [], set(), None
    for res in sel.residues:
        atoms = res.atoms & sel
        if _chain_of(atoms, src) != chain:
            continue
        key = (int(res.resid), str(getattr(res, "icode", "") or ""))
        if key in seen:
            raise AdapterError(f"chain {chain}: residue {key} occurs twice; chain identity or insertion codes may have "
                               f"been lost")
        seen.add(key)
        n = atoms.select_atoms("name N").positions
        bonded = None
        if out and prev_c is not None and len(prev_c) and len(n):
            d = np.linalg.norm(prev_c[:, None, :] - n[None, :, :], axis=-1).min()
            bonded = bool(d <= peptide_bond_max)
        out.append(ObservedResidue(key[0], key[1], str(res.resname), bonded))
        prev_c = atoms.select_atoms("name C").positions
    if not out:
        raise AdapterError(f"chain {chain!r} has no residues in selection {sel_text!r}")
    return out


def map_kinase_roles_in_universe(universe, chain: str, reference: Reference, numbering: Numbering,
                                 anchors: Mapping[str, AnchorLabel], **kwargs: Any) -> RoleMapping:
    """:func:`observed_chain` followed by :func:`map_kinase_roles`."""
    return map_kinase_roles(observed_chain(universe, chain), chain, reference, numbering, anchors, **kwargs)


# --- mapping -------------------------------------------------------------------------------------------------------------

def map_kinase_roles(observed: Sequence[ObservedResidue], chain: str, reference: Reference, numbering: Numbering,
                     anchors: Mapping[str, AnchorLabel], *, min_identity: float = MIN_IDENTITY,
                     min_aligned: int = MIN_ALIGNED, allow_role_mutations: bool = False) -> RoleMapping:
    """Map the descriptor roles of :data:`DERIVED_ROLES` to residues of ``observed`` (see module docstring).

    Returns a :class:`RoleMapping` or raises :class:`RoleMappingError`."""
    ref = reference.sequence
    obs = list(observed)
    if not obs:
        raise AdapterError("no observed residues")
    if set(anchors) != set(ANCHOR_MOTIFS):
        raise AdapterError(f"anchors must be exactly {list(ANCHOR_ORDER)}; got {sorted(anchors)} (the other roles are "
                           f"derived from the anchors)")
    keys = [(o.resid, o.icode) for o in obs]
    if len(set(keys)) != len(keys):
        raise AdapterError("observed residues contain duplicate (resid, icode) keys")
    letters, variants = [], []
    for o in obs:
        letters.append(ONE_LETTER.get(o.resname.upper(), "X"))
        if o.resname.upper() in VARIANT_NAMES:
            variants.append(_res(o) | {"maps_to": VARIANT_NAMES[o.resname.upper()]})
    obs_seq = "".join(letters)
    at_break = np.array([False] + [o.bonded_to_previous is not True for o in obs[1:]], dtype=bool)  # before obs[j]

    problems: list[str] = []
    report: dict[str, Any] = {
        "schema": "topofold_inputs.kinase_roles/1",
        "engineering_note": "Engineering output; no biological or scientific validation is implied.",
        "version": KINASE_ROLES_VERSION, "chain": chain, "reference": reference.describe(),
        "numbering": {"scheme": numbering.scheme, "reference_minus_label": numbering.reference_minus_label,
                      "description": numbering.description},
        "anchors_declared": {a: _label_json(v) for a, v in anchors.items()},
        "config": {"min_identity": min_identity, "min_aligned": min_aligned,
                   "allow_role_mutations": allow_role_mutations,
                   "scores": {"match": MATCH, "mismatch": MISMATCH, "insertion_open": INS_OPEN,
                              "insertion_extend": INS_EXT, "deletion_open_at_break": DEL_OPEN,
                              "deletion_extend_at_break": DEL_EXT, "deletion_open_between_bonded": BONDED_DEL_OPEN,
                              "deletion_extend_between_bonded": BONDED_DEL_EXT}},
    }

    score, pairs = _align(ref, obs_seq, at_break)
    o2r = {j: i for i, j in pairs}
    r2o = {i: j for i, j in pairs}
    report["alignment"] = aln = _alignment_report(ref, obs, obs_seq, pairs, score, variants)
    if aln["aligned_pairs"] < min_aligned or aln["identity"] < min_identity:
        problems.append(f"alignment to {reference.accession} ({reference.isoform}) has {aln['aligned_pairs']} aligned "
                        f"residues at identity {aln['identity']:.3f} (required ≥ {min_aligned} and ≥ {min_identity}); "
                        f"the reference or isoform may be wrong for this chain")

    # anchors -> reference indices (0-based internally, 1-based in the report)
    anchor_idx: dict[str, int] = {}
    checks: dict[str, Any] = {}
    for a in ANCHOR_ORDER:
        lab = anchors[a]
        if isinstance(lab, MotifWindow):
            if numbering.scheme != "reference":
                problems.append(f"{a}: a motif window needs scheme 'reference'")
                continue
            lo, hi = lab.first + numbering.reference_minus_label, lab.last + numbering.reference_minus_label
            hits = [p for p in range(max(lo, 1), min(hi, len(ref)) + 1) if _motif_ok(ref, a, p - 1)]
            checks[a] = {"window_reference_index": [lo, hi], "hits_reference_index": hits}
            if len(hits) != 1:
                problems.append(f"{a}: motif {_motif_text(a)} found {len(hits)} times in reference window {lo}–{hi} "
                                f"({hits}); {'ambiguous' if hits else 'not found'}")
                continue
            anchor_idx[a] = hits[0] - 1
        elif numbering.scheme == "reference":
            if not isinstance(lab, int) or isinstance(lab, bool):
                problems.append(f"{a}: scheme 'reference' takes an integer label, got {lab!r}")
                continue
            anchor_idx[a] = lab + numbering.reference_minus_label - 1
        else:
            resid, icode = (lab, "") if isinstance(lab, int) and not isinstance(lab, bool) else tuple(lab)
            if (resid, icode) not in keys:
                problems.append(f"{a}: author residue {resid}{icode} is not in chain {chain}")
                continue
            j = keys.index((resid, icode))
            if j not in o2r:
                problems.append(f"{a}: author residue {resid}{icode} is not aligned to the reference "
                                f"(insertion or unaligned terminus)")
                continue
            anchor_idx[a] = o2r[j]
            alt = _align(ref, obs_seq, at_break, forbid_pair=(o2r[j], j))[0]
            checks[a] = {"author_residue_unique_placement": alt < score}
            if alt >= score:
                problems.append(f"{a}: placement of author residue {resid}{icode} in the reference is ambiguous "
                                f"(an equally scoring alignment places it elsewhere)")
    for a, p in anchor_idx.items():
        motif = "".join(ref[p + k] if 0 <= p + k < len(ref) else "-" for k, _ in ANCHOR_MOTIFS[a])
        ok = _motif_ok(ref, a, p)
        checks.setdefault(a, {}).update(reference_index=p + 1, expected=_motif_text(a), found=motif, passed=ok)
        if not ok:
            problems.append(f"{a}: reference {reference.accession} has {motif} at index {p + 1}, expected "
                            f"{_motif_text(a)}")
    report["anchor_checks"] = checks
    if len(anchor_idx) == 4:
        order_ok = all(anchor_idx[x] < anchor_idx[y] for x, y in zip(ANCHOR_ORDER, ANCHOR_ORDER[1:]))
        report["anchor_order"] = {"expected": " < ".join(ANCHOR_ORDER), "passed": order_ok}
        report["anchor_spacing_reported_not_checked"] = {
            "alphaC_Glu - beta3_Lys": anchor_idx["alphaC_Glu"] - anchor_idx["beta3_Lys"],
            "DFG_Asp - HRD_H": anchor_idx["DFG_Asp"] - anchor_idx["HRD_H"]}
        if not order_ok:
            problems.append(f"anchor order contradicts the kinase fold: reference indices "
                            f"{ {a: anchor_idx[a] + 1 for a in ANCHOR_ORDER} } are not increasing in the order "
                            f"{' < '.join(ANCHOR_ORDER)}")
    if numbering.scheme == "reference" and any(not c.get("passed", True) for c in checks.values()):
        explicit = {a: v for a, v in anchors.items() if isinstance(v, int) and not isinstance(v, bool)}
        alts = [d for d in range(-10, 11) if d and len(explicit) == 4 and all(
            _motif_ok(ref, a, v + numbering.reference_minus_label + d - 1) for a, v in explicit.items())]
        report["diagnostic_offsets_not_applied"] = {
            "note": "extra offsets d (|d| ≤ 10) under which all four anchor motifs would hold; reported only, never "
                    "applied. If the labels use another numbering, declare it explicitly.",
            "d": alts}
        if alts:
            problems.append(f"all anchor motifs would hold with reference_minus_label = "
                            f"{[numbering.reference_minus_label + d for d in alts]}; not applied (declare the "
                            f"numbering explicitly if that is intended)")

    # roles and motif positions in the structure
    roles: dict[str, ResidueSpec | None] = {}
    role_rep: dict[str, Any] = {}
    positions: dict[int, list[str]] = {}
    if len(anchor_idx) == 4:
        for r, (a, k) in DERIVED_ROLES.items():
            positions.setdefault(anchor_idx[a] + k, []).append(r)
        for a, motif in ANCHOR_MOTIFS.items():
            for k, _ in motif:
                positions.setdefault(anchor_idx[a] + k, []).append(f"{a}_motif+{k}")
    ambiguity = {}
    for p, names in sorted(positions.items()):
        if not 0 <= p < len(ref):
            problems.append(f"{'/'.join(names)}: reference index {p + 1} is outside the reference (1–{len(ref)})")
            continue
        j = r2o.get(p)
        if j is not None:
            alt = _align(ref, obs_seq, at_break, forbid_pair=(p, j))[0]
        else:
            alt = _align(ref, obs_seq, at_break, require_match=p)[0]
        unique = alt < score
        ambiguity[str(p + 1)] = {"for": names, "unique": unique, "best_alternative_score": int(alt) if alt > _NEG // 2
                                 else None, "score": int(score)}
        if not unique:
            problems.append(f"{'/'.join(names)} (reference {p + 1}{ref[p]}): placement is ambiguous; an equally "
                            f"scoring alignment {'maps it elsewhere' if j is not None else 'observes it'}")
        if j is None and _unobserved_reason(p, pairs, obs).startswith("absent between"):
            problems.append(f"{'/'.join(names)} (reference {p + 1}{ref[p]}): absent between peptide-bonded structure "
                            f"residues, i.e. deleted in the construct; the roles cannot be resolved")
        if j is not None and obs_seq[j] != ref[p]:
            msg = (f"{'/'.join(names)}: structure residue {obs[j].resname} {obs[j].resid}{obs[j].icode} differs from "
                   f"reference {ref[p]}{p + 1}")
            if not allow_role_mutations:
                problems.append(msg + " (role or motif mutation; pass allow_role_mutations=True to accept it)")
        for r in (n for n in names if n in DERIVED_ROLES):
            if j is None:
                roles[r] = None
                role_rep[r] = {"reference_index": p + 1, "reference_aa": ref[p], "status": "unobserved",
                               "reason": _unobserved_reason(p, pairs, obs)}
            else:
                o = obs[j]
                roles[r] = ResidueSpec(chain, o.resid, o.icode, label=r)
                role_rep[r] = {"reference_index": p + 1, "reference_aa": ref[p], "resid": o.resid, "icode": o.icode,
                               "resname": o.resname,
                               "status": "observed" if obs_seq[j] == ref[p] else "observed_mutated"}
    report["roles"] = {r: role_rep[r] for r in DERIVED_ROLES if r in role_rep}
    report["ambiguity_checks"] = ambiguity
    report["problems"] = problems
    report["accepted"] = not problems
    if problems:
        raise RoleMappingError(problems, report)
    return RoleMapping(chain=chain, roles={r: roles[r] for r in DERIVED_ROLES}, report=report)


def _unobserved_reason(p: int, pairs, obs) -> str:
    """Why reference index p has no structure residue."""
    after = next((j for i, j in pairs if i > p), None)
    if not pairs or p < pairs[0][0] or p > pairs[-1][0]:
        return "outside the aligned span (not in the construct or not modelled)"
    if obs[after].bonded_to_previous is True:
        return "absent between peptide-bonded residues (deleted in the construct)"
    return ("inside the aligned span, at a chain break" if obs[after].bonded_to_previous is False
            else "inside the aligned span, bonding unknown")


def _motif_ok(ref: str, anchor: str, p: int) -> bool:
    return all(0 <= p + k < len(ref) and ref[p + k] in allowed for k, allowed in ANCHOR_MOTIFS[anchor])


def _motif_text(anchor: str) -> str:
    return "-".join("/".join(allowed) for _, allowed in ANCHOR_MOTIFS[anchor])


def _label_json(v: AnchorLabel):
    if isinstance(v, MotifWindow):
        return {"window": [v.first, v.last]}
    return v if isinstance(v, int) else list(v)


def _res(o: ObservedResidue) -> dict[str, Any]:
    return {"resid": o.resid, "icode": o.icode, "resname": o.resname}


def _alignment_report(ref, obs, obs_seq, pairs, score, variants) -> dict[str, Any]:
    o2r = {j: i for i, j in pairs}
    matched = sum(1 for i, j in pairs if ref[i] == obs_seq[j])
    first_j, last_j = (pairs[0][1], pairs[-1][1]) if pairs else (0, -1)
    unobserved, bonded_gaps, prev = [], [], None
    for i, j in pairs:
        if prev is not None and i > prev + 1:  # 1-based inclusive range of reference indices
            (bonded_gaps if obs[j].bonded_to_previous is True else unobserved).append([prev + 2, i])
        prev = i
    offsets: dict[int, int] = {}
    for i, j in pairs:
        offsets[i + 1 - obs[j].resid] = offsets.get(i + 1 - obs[j].resid, 0) + 1
    return {
        "score": int(score), "aligned_pairs": len(pairs), "identical_pairs": matched,
        "identity": matched / len(pairs) if pairs else 0.0, "observed_residues": len(obs),
        "reference_span": [pairs[0][0] + 1, pairs[-1][0] + 1] if pairs else None,
        "mutations": [_res(obs[j]) | {"reference_index": i + 1, "reference_aa": ref[i], "observed_aa": obs_seq[j]}
                      for i, j in pairs if ref[i] != obs_seq[j]],
        "insertions": [_res(obs[j]) for j in range(first_j, last_j + 1) if j not in o2r],
        "unaligned_terminal": [_res(obs[j]) for j in range(len(obs)) if j < first_j or j > last_j],
        "unobserved_reference_ranges": unobserved,
        "absent_between_bonded_residues": bonded_gaps,
        "chain_breaks_before": [_res(o) for o in obs[1:] if o.bonded_to_previous is False],
        "bonding_unknown": sum(1 for o in obs[1:] if o.bonded_to_previous is None),
        "insertion_code_residues": [_res(o) for o in obs if o.icode],
        "variant_residue_names": variants,
        "author_numbering_diagnostic": {
            "reference_index_minus_author_resid": {str(k): v for k, v in sorted(offsets.items())},
            "note": "reported for inspection only; author numbers are never used to place residues"},
    }


def _align(ref: str, obs: str, at_break: np.ndarray, *, forbid_pair: tuple[int, int] | None = None,
           require_match: int | None = None) -> tuple[int, list[tuple[int, int]]]:
    """Local alignment with integer scores. ``at_break[j]``: obs[j] is not known to be bonded to obs[j-1].

    Returns (score, [(ref index, obs index), ...] in order). ``forbid_pair`` removes one match; ``require_match=p``
    keeps only alignments in which reference index p is matched."""
    n, m = len(ref), len(obs)
    r = np.frombuffer(ref.encode(), dtype=np.uint8)
    o = np.frombuffer(obs.encode(), dtype=np.uint8)
    S = np.where(r[:, None] == o[None, :], MATCH, MISMATCH).astype(np.int64)
    S[:, o == ord("X")] = MISMATCH
    if forbid_pair is not None:
        S[forbid_pair] = _NEG
    Mx = np.full((n + 1, m + 1), _NEG, dtype=np.int64)
    Dx = np.full((n + 1, m + 1), _NEG, dtype=np.int64)
    Ix = np.full((n + 1, m + 1), _NEG, dtype=np.int64)
    # D[i, j]: ref i-1 skipped between obs j-1 and obs j (internal only); cost depends on whether they are bonded
    d_open = np.full(m + 1, _NEG, dtype=np.int64)
    d_ext = np.full(m + 1, _NEG, dtype=np.int64)
    d_open[1:m] = np.where(at_break[1:m], DEL_OPEN, BONDED_DEL_OPEN)
    d_ext[1:m] = np.where(at_break[1:m], DEL_EXT, BONDED_DEL_EXT)
    for i in range(1, n + 1):
        start = 0 if require_match is None or i - 1 <= require_match else _NEG
        prev = np.maximum(np.maximum(Mx[i - 1, :-1], Dx[i - 1, :-1]), np.maximum(Ix[i - 1, :-1], start))
        Mx[i, 1:] = prev + S[i - 1]
        Dx[i] = np.maximum(np.maximum(Mx[i - 1] + d_open, Dx[i - 1] + d_ext), Ix[i - 1] + d_open)
        if require_match is not None and i - 1 == require_match:
            Dx[i] = _NEG
        a = np.maximum(Mx[i], Dx[i])
        b = a[:-1] + INS_OPEN - INS_EXT * np.arange(1, m + 1)
        Ix[i, 1:] = INS_EXT * np.arange(1, m + 1) + np.maximum.accumulate(b)
        Ix[i, 0] = _NEG
        np.maximum(Mx[i], _NEG, out=Mx[i])
        np.maximum(Dx[i], _NEG, out=Dx[i])
        np.maximum(Ix[i], _NEG, out=Ix[i])
    ends = Mx.copy()
    if require_match is not None:
        ends[: require_match + 1] = _NEG   # the last matched reference index must be ≥ p
    best = int(ends.max())
    if best <= _NEG // 2:
        return _NEG, []
    i, j = (int(x) for x in np.unravel_index(int(np.argmax(ends)), ends.shape))
    pairs, state = [], "M"
    while True:
        if state == "M":
            pairs.append((i - 1, j - 1))
            v = Mx[i, j] - S[i - 1, j - 1]
            start = 0 if require_match is None or i - 1 <= require_match else _NEG
            if start == 0 and v == 0:   # the alignment starts at this match
                break
            i, j = i - 1, j - 1
            state = "M" if Mx[i, j] == v else "D" if Dx[i, j] == v else "I"
        elif state == "D":
            v = Dx[i, j]
            i -= 1
            state = "M" if Mx[i, j] + d_open[j] == v else "D" if Dx[i, j] + d_ext[j] == v else "I"
        else:
            v = Ix[i, j]
            j -= 1
            state = "M" if Mx[i, j] + INS_OPEN == v else "I" if Ix[i, j] + INS_EXT == v else "D"
    return best, pairs[::-1]


__all__ = ["ANCHOR_MOTIFS", "ANCHOR_ORDER", "DERIVED_ROLES", "KINASE_ROLES_VERSION", "MotifWindow", "Numbering",
           "ObservedResidue", "PEPTIDE_BOND_MAX", "Reference", "RoleMapping", "RoleMappingError", "map_kinase_roles",
           "map_kinase_roles_in_universe", "observed_chain"]

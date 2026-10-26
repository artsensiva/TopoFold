"""Runs the existing kinase_roles helper on PIM1 chain A of each input (needs the adapter interpreter with MDAnalysis).

    python run_roles.py CONFIG.json OUT.json

CONFIG: accession, isoform, sequence, source, files {id: path}, anchors {beta3_Lys, alphaC_Glu, HRD_H: int; DFG_Asp: [first, last] window}, allow_role_mutations.
OUT: per file, the helper's full report (accepted or refused with every problem). Nothing is changed to make a mapping pass.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.ca_adapter import load_universe  # noqa: E402
from topofold_inputs.kinase_roles import MotifWindow, Numbering, Reference, RoleMappingError, map_kinase_roles_in_universe  # noqa: E402


def main():
    cfg = json.loads(Path(sys.argv[1]).read_text())
    ref = Reference(cfg["accession"], cfg["isoform"], cfg["sequence"], cfg["source"])
    a = cfg["anchors"]
    anchors = {"beta3_Lys": a["beta3_Lys"], "alphaC_Glu": a["alphaC_Glu"], "HRD_H": a["HRD_H"], "DFG_Asp": MotifWindow(*a["DFG_Asp"])}
    out = {}
    for fid, path in cfg["files"].items():
        u = load_universe(path)
        try:
            m = map_kinase_roles_in_universe(u, "A", ref, Numbering.reference_index(), anchors, allow_role_mutations=cfg.get("allow_role_mutations", False))
            out[fid] = {"accepted": True, "report": m.report}
        except RoleMappingError as e:
            out[fid] = {"accepted": False, "problems": e.problems, "report": e.report}
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1, default=str) + "\n")
    print({k: v["accepted"] for k, v in out.items()})


if __name__ == "__main__":
    main()

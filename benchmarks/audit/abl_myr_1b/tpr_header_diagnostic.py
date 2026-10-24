"""Post-run diagnostic for the ABL-MYR-1b K1 failure: what each md.tpr header declares, and what the MDAnalysis Universe
exposes. Not part of any frozen protocol; reads the a1 extraction only, computes no comparison.

    PYTHONPATH=python <python with MDAnalysis> benchmarks/audit/abl_myr_1b/tpr_header_diagnostic.py
"""

import json
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EXT = ROOT / "benchmarks/audit/abl_myr_1/data_dryad/extracted_1b_a1"
SETS = {"Set1": range(11, 21), "Set2": range(21, 31)}
CONDITIONS = {"Myr": "kinase_kinked_Myr", "noMyr": "kinase_kinked_noMyr"}


def main() -> int:
    import MDAnalysis as mda
    from MDAnalysis.topology.tpr import utils as tu
    out = {"purpose": "diagnostic after ABL-MYR-1b-CHECK-a1 failed K1; not part of a frozen protocol",
           "mdanalysis_version": mda.__version__, "runs": {}}
    for s, reps in SETS.items():
        for c, model in CONDITIONS.items():
            for r in reps:
                p = EXT / model / model / s / str(r) / "md.tpr"
                th = tu.read_tpxheader(tu.TPXUnpacker(p.read_bytes()))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    u = mda.Universe(str(p))
                out["runs"][f"{s}:{c}:{r}"] = {
                    "gromacs_version": th.ver_str.decode(), "tpx_version": int(th.fver), "natoms": int(th.natoms),
                    "ngtc": int(th.ngtc), "header_bBox": int(th.bBox), "header_bX": int(th.bX), "header_bV": int(th.bV),
                    "universe_dimensions": None if u.dimensions is None else [float(x) for x in u.dimensions],
                    "universe_n_bonds": int(len(u.atoms.bonds))}
    runs = out["runs"].values()
    out["summary"] = {"runs": len(out["runs"]), "header_box_declared": sum(v["header_bBox"] == 1 for v in runs),
                      "universe_box_exposed": sum(v["universe_dimensions"] is not None for v in runs),
                      "gromacs_versions": sorted({v["gromacs_version"] for v in runs}),
                      "natoms": sorted({v["natoms"] for v in runs}), "n_bonds": sorted({v["universe_n_bonds"] for v in runs})}
    (HERE / "results/ABL-MYR-1b_tpr_header_diagnostic.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["summary"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

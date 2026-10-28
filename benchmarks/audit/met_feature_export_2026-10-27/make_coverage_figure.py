"""Coverage figure from results/features_structure.csv (the T and C blocks share one window support, so one row per drug). Run with the repository .venv."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(HERE / "results/features_structure.csv")))
drugs = list(dict.fromkeys(r["drug"] for r in rows))
pdb = {r["drug"]: r["pdb"] for r in rows}
pos = sorted({int(r["uniprot_position"]) for r in rows})
CAT = ["ok", "no_side_chain", "missing_residue", "engineered_residue", "incomplete_side_chain", "missing_CA", "other"]
COL = ["#2a78d6", "#9cc3f0", "#eb6834", "#b3261e", "#8a5cc2", "#2b2a27", "#6b6a66"]


def cat(reason):
    k = reason.split(";")[0].split(":")[0]
    return CAT.index(k) if k in CAT else CAT.index("other")


grid, labels = [], []
for d in drugs:
    rs = {int(r["uniprot_position"]): r for r in rows if r["drug"] == d}
    for blk, col in (("ligand context L", "L_reason"), ("window T and C", "T_reason")):
        grid.append([cat(rs[p][col]) for p in pos])
        labels.append(f"{d} ({pdb[d]}) {blk}")
fig, ax = plt.subplots(figsize=(10.5, 4.6))
ax.imshow(np.array(grid), aspect="auto", cmap=ListedColormap(COL), vmin=0, vmax=len(CAT) - 1, interpolation="nearest", extent=[pos[0] - 0.5, pos[-1] + 0.5, len(grid) - 0.5, -0.5])
ax.set_yticks(range(len(grid)))
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel("MET residue (UniProt P08581 numbering; DMS span 1059-1345)")
ax.legend(handles=[Patch(color=c, label=n) for c, n in zip(COL, CAT)], fontsize=7, ncol=7, loc="upper center", bbox_to_anchor=(0.5, 1.13), frameon=False)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(HERE / "figures/met_feature_coverage.png", dpi=120)
print("figure written")

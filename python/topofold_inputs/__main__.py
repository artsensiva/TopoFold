"""Command line: ``python -m topofold_inputs TOPOLOGY [TRAJECTORY ...] --run-id ID --out DIR``."""

import argparse
import json
import sys

from .ca_adapter import ALTLOC_POLICIES, MODIFIED_POLICIES, PBC_POLICIES, AdapterError, AltlocTieError, extract_ca, load_universe


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m topofold_inputs", description=__doc__)
    ap.add_argument("topology")
    ap.add_argument("trajectories", nargs="*")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--out", required=True, help="output directory for <run-id>.npz and <run-id>.json")
    ap.add_argument("--selection", default="protein")
    ap.add_argument("--correspondence", help="JSON file: list of {chain, resid, icode, label}")
    ap.add_argument("--on-missing", choices=("mark", "reject"), default="mark")
    ap.add_argument("--chain-attribute", choices=("auto", "chainID", "segid"), default="auto")
    ap.add_argument("--altloc", choices=ALTLOC_POLICIES, default="reject")
    ap.add_argument("--pbc", choices=PBC_POLICIES, default="none")
    ap.add_argument("--max-ca-ca", type=float, default=4.2, help="geometric break threshold, Å")
    ap.add_argument("--start", type=int, default=None, help="first source frame (slice start)")
    ap.add_argument("--stop", type=int, default=None, help="slice stop (exclusive)")
    ap.add_argument("--step", type=int, default=None, help="slice step")
    ap.add_argument("--modified-residues", choices=MODIFIED_POLICIES, default="none",
                    help="phospho_pdb: admit SEP, TPO and PTR documented by MODRES records of a PDB-format file (default none: historical behaviour)")
    ap.add_argument("--error-json", default=None,
                    help="if the altloc policy meets alternate locations tied on occupancy, also write them here (JSON, kind altloc_tie)")
    args = ap.parse_args(argv)

    corr = None
    if args.correspondence:
        with open(args.correspondence) as fh:
            corr = json.load(fh)
    try:
        u = load_universe(args.topology, *args.trajectories)
        res = extract_ca(u, args.run_id, selection=args.selection, correspondence=corr, on_missing=args.on_missing,
                         chain_attribute=args.chain_attribute, altloc=args.altloc, pbc=args.pbc, modified_residues=args.modified_residues,
                         max_ca_ca_angstrom=args.max_ca_ca, start=args.start, stop=args.stop, step=args.step)
    except AltlocTieError as exc:
        print(f"rejected: {exc}", file=sys.stderr)
        if args.error_json:
            with open(args.error_json, "w") as fh:
                json.dump({"kind": "altloc_tie", "policy": args.altloc, "message": str(exc), "ties": exc.ties}, fh, indent=1)
        return 2
    except AdapterError as exc:
        print(f"rejected: {exc}", file=sys.stderr)
        return 2
    npz, js = res.export(args.out)
    s = res.report["summary"]
    print(f"{s['frames']} frames x {s['residues']} residues ({s['present']} present, {s['missing']} missing), "
          f"{s['breaks']} breaks, {s['segments']} segments")
    print(f"wrote {npz}\nwrote {js}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

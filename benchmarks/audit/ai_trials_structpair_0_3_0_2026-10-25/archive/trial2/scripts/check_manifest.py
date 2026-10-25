"""Self-devised: re-hash every file listed in manifest outputs_sha256 and inputs (trial helper, not part of TopoFold)."""
import hashlib, json, sys, pathlib
d = pathlib.Path(sys.argv[1]); m = json.loads((d/"manifest.json").read_text())
ok = bad = 0
for name, h in m["outputs_sha256"].items():
    p = d/name
    got = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
    s = "OK" if got == h else f"MISMATCH (got {got})"
    ok += got == h; bad += got != h
    print(f"output {name}: {s}")
print(json.dumps(m["inputs"], indent=1)[:1500])
print(f"outputs: {ok} match, {bad} mismatch")

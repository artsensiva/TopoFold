"""Self-devised pre-run inspection: Calpha residues, gaps, altlocs, resname differences (plain text parsing)."""
import sys, collections
def ca(path, chain):
    res = collections.OrderedDict(); alts = collections.defaultdict(list)
    for l in open(path):
        if l.startswith(("ATOM","HETATM")) and l[12:16].strip()=="CA" and l[21]==chain and l[17:20]!=" CA":
            key=(int(l[22:26]), l[26].strip())
            res.setdefault(key, l[17:20]); 
            if l[16]!=" ": alts[key].append((l[16], l[17:20], float(l[54:60])))
        if l.startswith("ENDMDL"): break
    return res, alts
def ranges(nums):
    out=[]; 
    for n in nums:
        if out and n==out[-1][1]+1: out[-1][1]=n
        else: out.append([n,n])
    return ",".join(f"{a}-{b}" for a,b in out)
sets={}
for path,ch in [("inputs/1OPJ.pdb","A"),("inputs/1OPJ.pdb","B"),("inputs/1OPK.pdb","A")]:
    r,a=ca(path,ch); sets[(path,ch)]=r
    print(f"{path} {ch}: {len(r)} CA, ranges {ranges([k[0] for k in r])}, inscodes {[k for k in r if k[1]]}")
    print("   altloc CA:", {f'{k[0]}{k[1]}':v for k,v in a.items()})
common=set.intersection(*[set(r) for r in sets.values()])
print("common CA positions:", len(common), ranges(sorted(k[0] for k in common)))
ref=sets[("inputs/1OPK.pdb","A")]
for key in sorted(common):
    names={sets[s][key] for s in sets}
    if len(names)>1: print("resname differs at", key, {f"{s[0][-8:-4]}{s[1]}":sets[s][key] for s in sets})

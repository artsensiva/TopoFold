"""Self-devised independent check (numpy only): Ca-Ca bond lengths and displacements after a fit that excludes 395-420."""
import numpy as np, sys
def ca(path, ch):
    d={}
    for l in open(path):
        if l.startswith("ATOM") and l[12:16]==" CA " and l[21]==ch: d[int(l[22:26])]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return d
S={"1OPJ_A":ca("inputs/1OPJ.pdb","A"),"1OPJ_B":ca("inputs/1OPJ.pdb","B"),"1OPK_A":ca("inputs/1OPK.pdb","A")}
R=list(range(248,530))
for k,s in S.items():
    b=[np.linalg.norm(s[i+1]-s[i]) for i in range(398,418)]
    print(k,"CA(i)-CA(i+1) 398..417 min/max %.2f/%.2f"%(min(b),max(b)))
def fit(P,Q):
    pc,qc=P.mean(0),Q.mean(0); H=(P-pc).T@(Q-qc); U,_,Vt=np.linalg.svd(H); d=np.sign(np.linalg.det(Vt.T@U.T))
    Rm=Vt.T@np.diag([1,1,d])@U.T; return lambda X:(X-pc)@Rm.T+qc
core=[i for i in R if not 395<=i<=420 and i>255]
for a,b in [("1OPJ_A","1OPK_A"),("1OPJ_B","1OPK_A"),("1OPJ_A","1OPJ_B")]:
    f=fit(np.array([S[a][i] for i in core]),np.array([S[b][i] for i in core]))
    dev=np.array([np.linalg.norm(f(S[a][i][None])[0]-S[b][i]) for i in R])
    cr=np.sqrt(np.mean([dev[R.index(i)]**2 for i in core]))
    print(f"{a} vs {b}: core RMSD {cr:.2f} A; disp at 400..415:", " ".join(f"{i}:{dev[R.index(i)]:.1f}" for i in range(400,416)))

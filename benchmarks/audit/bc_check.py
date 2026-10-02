"""Audit M2/M3: bimodality coefficient (port of bimodality.rs) on unimodal, rare-event and angle-wrapping cases."""
import numpy as np
def bc_topofold(x):
    # port of StreamingMoments::bimodality_coefficient (population g1,g2, clamp [0,1])
    x=np.asarray(x,float); n=len(x); m=x.mean(); d=x-m
    m2=(d**2).mean(); m3=(d**3).mean(); m4=(d**4).mean()
    g=m3/m2**1.5; k=m4/m2**2-3
    c=3*(n-1)**2/((n-2)*(n-3))
    return float(np.clip((g*g+1)/(k+c),0,1))
wrap=lambda a:(a+np.pi)%(2*np.pi)-np.pi
rng=np.random.default_rng(1); N=2500; W=5  # 8-residue window -> 5 torsions averaged
print("== 1. Branch-cut artefact: unimodal beta-strand torsion near +/-180 deg ==")
for mu,sd in [(-170,15),(-175,15),(180,10),(-120,15),(50,10)]:
    tau=wrap(np.radians(rng.normal(mu,sd,size=(N,W))))   # per-residue tau, wrapped to (-pi,pi]
    wm=tau.mean(axis=1)                                   # TopoFold: linear mean over window
    circ=np.angle(np.exp(1j*tau).mean(axis=1))            # correct circular mean
    print(f"tau ~ N({mu:+d} deg, {sd} deg): BC(linear mean, as in TopoFold)={bc_topofold(wm):.3f} | BC(circular mean)={bc_topofold(wrap(circ-np.radians(mu))):.3f}")
print("\n== 2. BC = 1 for ANY two-point / rare-event distribution (no bimodality needed) ==")
for p in [0.5,0.1,0.02,0.005]:
    x=rng.normal(0,0.05,N); x[rng.random(N)<p]+=1.0
    print(f"Gaussian + {p*100:.1f}% rare excursions: BC={bc_topofold(x):.3f}")
print("\n== 3. Unimodal skewed distributions exceed the 0.555 / 0.6 thresholds ==")
for name,x in [("exponential",rng.exponential(1,N)),("lognormal s=0.8",rng.lognormal(0,0.8,N)),("gamma k=1",rng.gamma(1,1,N)),("chi2 df=1",rng.chisquare(1,N)),("Gaussian",rng.normal(size=N))]:
    print(f"{name:18s}: BC={bc_topofold(x):.3f}")
print("\n== 4. Branch-cut: one residue of an 8-residue window fluctuating around 180 deg, rest rigid ==")
for sd in [5,10,15,20]:
    tau=np.zeros((N,W)); tau[:,:]=np.radians(rng.normal(50,5,size=(N,W)))
    tau[:,2]=wrap(np.radians(rng.normal(180,sd,N)))
    wm=tau.mean(axis=1)
    fixed=tau.copy(); fixed[:,2]=wrap(fixed[:,2]-np.pi)+np.pi   # unwrapped around pi
    print(f"sd={sd:2d} deg: BC(TopoFold linear)={bc_topofold(wm):.3f}  vs  BC(properly unwrapped)={bc_topofold(fixed.mean(axis=1)):.3f}")

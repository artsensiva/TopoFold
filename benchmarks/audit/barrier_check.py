"""Audit C2/C3: the former compute_pmf_barrier estimator on unimodal data, and label circularity."""
import numpy as np
from sklearn.metrics import silhouette_score
# port of compute_pmf_barrier() from benchmarks/run_honest_real_md_suite.py
def pmf_barrier(c2d,labels):
    sil=float(silhouette_score(c2d,labels))
    c1=c2d[labels==0].mean(0); c2=c2d[labels==1].mean(0); u=(c2-c1)/np.linalg.norm(c2-c1)
    h,_=np.histogram(c2d@u,bins=50,density=True); h=h[h>0]
    return sil,float(-np.log(max(h.min()/h.max(),1e-12)))
rng=np.random.default_rng(7)
for N in [1200,2500]:
    res=[]
    for rep in range(20):
        X=rng.normal(size=(N,2))                       # strictly UNIMODAL Gaussian, no states at all
        lab=(X[:,0]>np.median(X[:,0])).astype(int)     # median split, as in the BPTI/Mpro scripts
        res.append(pmf_barrier(X,lab))
    s=np.array(res)
    print(f"N={N}: unimodal Gaussian -> 'barrier' = {s[:,1].mean():.2f} +/- {s[:,1].std():.2f} kT, silhouette = {s[:,0].mean():.2f}")
# circular labels: labels defined by the same features used for silhouette (run_real_bpti_validation.py)
fa=rng.normal(1,0.3,2500); fb=rng.normal(1,0.3,2500)   # independent noise, no structure
lab=(fb<fa).astype(int)
print(f"Labels = (frechet_b < frechet_a) on pure noise: silhouette in (frechet_a, frechet_b) space = {silhouette_score(np.c_[fa,fb],lab):.2f}")

"""Historical audit reproduction of the TopoFold v0.8.3 pre-fix writhe implementation.

This script does not validate the current production implementation. It preserves the
original normalization defect (4π instead of 2π) for v0.8.3 bug reproduction."""
import numpy as np
# Port of TopoFold writhe.rs
def solid(a,b,c):
    det = a@np.cross(b,c); den = 1+a@b+b@c+c@a
    return 2*np.arctan2(det,den)
def seg_w(r1,r2,r3,r4):
    u=lambda v: v/np.linalg.norm(v)
    u13,u14,u24,u23=u(r3-r1),u(r4-r1),u(r4-r2),u(r3-r2)
    return (solid(u13,u14,u24)+solid(u13,u24,u23))/(4*np.pi)
def wr_topofold(X):
    n=len(X); s=0
    for i in range(n-2):
        for j in range(i+2,n-1):
            s+=seg_w(X[i],X[i+1],X[j],X[j+1])
    return s
# Reference: Klenin & Langowski 2000 exact segment-pair formula (method 1a)
def kl_pair(r1,r2,r3,r4):
    r13,r14,r23,r24=r3-r1,r4-r1,r3-r2,r4-r2
    n1=np.cross(r13,r14); n2=np.cross(r14,r24); n3=np.cross(r24,r23); n4=np.cross(r23,r13)
    ns=[v/np.linalg.norm(v) for v in (n1,n2,n3,n4)]
    om=np.arcsin(np.clip(ns[0]@ns[1],-1,1))+np.arcsin(np.clip(ns[1]@ns[2],-1,1))+np.arcsin(np.clip(ns[2]@ns[3],-1,1))+np.arcsin(np.clip(ns[3]@ns[0],-1,1))
    sgn=np.sign(np.cross(r4-r3,r2-r1)@r13)
    return om*sgn/(4*np.pi)
def wr_kl(X):
    n=len(X); s=0
    for i in range(n-2):
        for j in range(i+2,n-1):
            s+=kl_pair(X[i],X[i+1],X[j],X[j+1])
    return 2*s   # Wr = sum over ordered pairs i!=j = 2*sum_{i<j}
# Brute-force numerical Gauss double integral as an independent check
def wr_numeric(X,m=40):
    n=len(X); s=0; t=(np.arange(m)+0.5)/m
    for i in range(n-1):
        a0,a1=X[i],X[i+1]; da=a1-a0
        P=a0+np.outer(t,da)
        for j in range(n-1):
            if abs(i-j)<=1: continue
            b0,b1=X[j],X[j+1]; db=b1-b0
            Q=b0+np.outer(t,db)
            D=P[:,None,:]-Q[None,:,:]
            num=np.einsum('k,ijk->ij',np.cross(da,db),D)
            s+=(num/np.linalg.norm(D,axis=2)**3).sum()/m/m
    return s/(4*np.pi)
th=np.radians(100); X=np.array([[2.3*np.cos(i*th),2.3*np.sin(i*th),1.5*i] for i in range(25)])
rng=np.random.default_rng(0); Y=np.cumsum(rng.normal(size=(30,3))*2,axis=0)
for name,C in [("ideal alpha-helix CA (25)",X),("random walk (30)",Y)]:
    a,b,c=wr_topofold(C),wr_kl(C),wr_numeric(C)
    print(f"{name}: TopoFold={a:+.4f}  Klenin-Langowski={b:+.4f}  numeric Gauss={c:+.4f}  ratio TF/numeric={a/c:.3f}")

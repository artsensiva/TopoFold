"""Audit C2: global vs local baselines on the project's own synthetic bistable benchmark.

Uses benchmarks/generate_bistable_trajectory.py unchanged (true labels: frames 0..999 = A, 1000..1999 = B).
Pure NumPy/scikit-learn; does not require the compiled topofold module.
"""
import os
import sys

import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import generate_bistable_trajectory as g  # noqa: E402

out = g.generate_bistable_ensemble()
traj = out[0]
labels = np.array([0] * (len(traj) // 2) + [1] * (len(traj) - len(traj) // 2))


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    return Pc @ (U @ np.diag([1, 1, d]) @ Vt)


def sil(X):
    return silhouette_score(PCA(2).fit_transform(X), labels)


def kappa_tau(X):
    T = np.diff(X, axis=0)
    T /= np.linalg.norm(T, axis=1)[:, None]
    k = np.arccos(np.clip((T[:-1] * T[1:]).sum(1), -1, 1))
    B = np.cross(T[:-1], T[1:])
    B /= np.linalg.norm(B, axis=1)[:, None]
    t = np.arctan2((np.cross(B[:-1], B[1:]) * T[1:-1]).sum(1), (B[:-1] * B[1:]).sum(1))
    return np.r_[k, np.sin(t), np.cos(t)]


ref = traj[0]
loop = list(range(24, 37))
iu = np.triu_indices(len(loop), 2)
glob = np.array([kabsch(f, ref) for f in traj]).reshape(len(traj), -1)
loc = np.array([kabsch(f[loop], ref[loop]) for f in traj]).reshape(len(traj), -1)
dist = np.array([np.linalg.norm(f[loop][:, None] - f[loop][None], axis=2)[iu] for f in traj])
inv = np.array([kappa_tau(f[loop]) for f in traj])

print(f"Global Cartesian PCA (all residues, global Kabsch): S = {sil(glob):.3f}")
print(f"Local Cartesian PCA (loop 25-37, local Kabsch):     S = {sil(loc):.3f}")
print(f"Loop Ca-Ca distances + PCA:                         S = {sil(dist):.3f}")
print(f"Loop (kappa, tau) + PCA:                            S = {sil(inv):.3f}")

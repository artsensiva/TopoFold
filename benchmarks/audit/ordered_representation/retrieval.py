"""Label-independent, tie-aware k-nearest-neighbour purity (protocol §6, amendment B).

The query is excluded by index. Candidates within tol = REL_TOL·(1 + d_k) of the k-th smallest
distance d_k form the boundary tie set; the score is the expected purity under uniform selection
of the remaining neighbours from that set. Neither labels nor frame order break ties.
"""
import numpy as np

REL_TOL = 1e-12


def knn_purity(D, labels, k, rel_tol=REL_TOL, per_query=False):
    D = np.array(D, dtype=np.float64, copy=True)
    labels = np.asarray(labels)
    n = D.shape[0]
    assert D.shape == (n, n) and 1 <= k <= n - 1
    np.fill_diagonal(D, np.inf)
    dk = np.partition(D, k - 1, axis=1)[:, k - 1]
    tol = rel_tol * (1.0 + dk)
    less = D < (dk - tol)[:, None]
    eq = np.abs(D - dk[:, None]) <= tol[:, None]
    same = labels[None, :] == labels[:, None]
    np.fill_diagonal(same, False)
    n_less, n_eq = less.sum(1), eq.sum(1)
    need = k - n_less
    assert np.all(need >= 1) and np.all(n_eq >= need)
    purity = ((same & less).sum(1) + need * (same & eq).sum(1) / n_eq) / k
    chance = same.sum(1) / (n - 1)
    out = {'purity': float(purity.mean()), 'tie_frequency': float(np.mean(n_eq > need)),
           'mean_boundary_tie_size': float(n_eq.mean()), 'chance': float(chance.mean()),
           'n_queries': int(n)}
    if per_query:
        out['per_query_purity'] = purity
    return out

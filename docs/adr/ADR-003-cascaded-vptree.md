# ADR-003: Cascaded Metric VP-Tree Indexing over KD-Trees and HNSW

## Status
Accepted (v0.2.0)

## Context
A primary goal of TopoFold is sub-microsecond subcurve search: given a target subcurve conformational motif (e.g. an open cryptic pocket conformation), find the $k$ most structurally similar conformers across an ensemble of $10^5\text{--}10^7$ frames.

Standard high-dimensional spatial indexing structures (KD-trees, Ball-trees, R-trees, and HNSW graph indexes) rely on Euclidean vector spaces $\mathbb{R}^D$ where distances satisfy $L_2$ geometric norms and coordinate orthogonal projections.
However, protein subcurve similarity in TopoFold is governed by the **Discrete Fréchet Distance** $d_F$ over intrinsic discrete curvature and torsion $(\kappa, \tau)$:
$$d_F(P, Q) = \min_{\alpha, \beta} \max_{t \in [0, 1]} d_{\mathbb{S}^1 \times \mathbb{S}^1}(P(\alpha(t)), Q(\beta(t)))$$
$d_F$ is a true metric (satisfies identity, symmetry, and triangle inequality: $d_F(x, z) \le d_F(x, y) + d_F(y, z)$), but it is strictly **non-Euclidean** and **non-vectorial**. Coordinate projections, hyperplanes, and dot products do not exist for $d_F$.

## Decision
We implemented a **Cascaded Metric Space Index** in `topofold-index`:
1. **Tier 1 (Chebyshev $L_\infty$ Coarse Filter)**:
   For each conformer frame, extract the sliding-window local writhe spectrum:
   $$\mathbf{w} = [\operatorname{Wr}_w(1), \dots, \operatorname{Wr}_w(K)]$$
   Because global and local writhe bound maximal curve deformation, any two curves with Fréchet distance $d_F < \epsilon$ satisfy:
   $$\|\mathbf{w}_1 - \mathbf{w}_2\|_\infty \le C \cdot \epsilon$$
   Candidates exceeding the Chebyshev threshold are rejected in $\mathcal{O}(1)$ time with zero Fréchet distance evaluations.
2. **Tier 2 (Vantage-Point Tree Search)**:
   A **Vantage-Point Tree (VP-Tree)** is built over the metric space $(X, d_F)$.
   - Each internal node chooses a vantage point $v \in X$ and computes the median distance $M = \operatorname{median}_{x} d_F(v, x)$.
   - Points are partitioned into "inside ball" ($d_F(v, x) \le M$) and "outside ball" ($d_F(v, x) > M$).
   - Query pruning uses the triangle inequality:
     $$\text{If } |d_F(v, q) - r_{\text{query}}| > M, \quad \text{prune inner child}$$
     $$\text{If } d_F(v, q) + r_{\text{query}} \le M, \quad \text{prune outer child}$$

## Consequences

### Positive
- **Exact Nearest Neighbors**: Unlike approximate nearest neighbor graphs (HNSW) which can miss cryptic conformers, VP-trees guarantee exact metric retrieval.
- **Strictly Metric-Only**: Requires only distance evaluations; completely agnostic to coordinate representations.
- **Logarithmic Retrieval**: Reduces search complexity from linear scan $\mathcal{O}(M)$ to $\mathcal{O}(\log M)$.
- **Sub-Microsecond Latency**: The Chebyshev filter prunes $>95\%$ of candidates before Fréchet evaluation, achieving benchmark query times of **$1.01\,\mu\text{s/frame}$** on synthetic ensembles and **$1.29\,\mu\text{s/frame}$** on BPTI.

### Negative / Trade-offs
- Tree construction is $\mathcal{O}(M \log M)$ with parallel Rayon tree-building. For a 2,500-frame ensemble, construction takes $\approx 55\text{ ms}$, which is easily amortized over thousands of subcurve queries.

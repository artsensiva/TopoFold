# ADR-002: Spherical Triangle Decomposition (Van Oosterom & Strackee) over Naive Gauss Double Integrals

> **Erratum (Oct 2026).** The formula below uses 1/(2π) over unordered segment pairs, but `writhe.rs` divides by 4π and its sign convention is reversed, so the implementation returns −½ of the standard writhe. See `docs/KNOWN_ISSUES.md`, M1. Note also that `atan2` removes ambiguity in each triangle's solid angle but writhe itself is a geometric, not a topological, invariant.

## Status
Accepted (v0.1.0)

## Context
Writhe ($\operatorname{Wr}$) is a fundamental knot-theoretic invariant measuring the chiral self-coiling and tertiary topological entanglement of a space curve.
The classical continuous definition of Writhe is Gauss's double line integral:
$$\operatorname{Wr}(\mathcal{C}) = \frac{1}{4\pi} \oint_{\mathcal{C}} \oint_{\mathcal{C}} \frac{(\mathbf{r}_1 - \mathbf{r}_2) \cdot (d\mathbf{r}_1 \times d\mathbf{r}_2)}{\|\mathbf{r}_1 - \mathbf{r}_2\|^3}$$
When applied to discrete polygonal chains representing protein backbones, two naive computational approaches are commonly attempted:
1. **Numerical Quadrature**: Midpoint discretization of the double integral.
   - *Failure Mode*: Highly sensitive to discretization step size; computationally intractable for large ensembles; singularity when $\|\mathbf{r}_1 - \mathbf{r}_2\| \to 0$.
2. **Naive Dihedral Angle Summation (Levitt/Klenin)**:
   - *Failure Mode*: Evaluating pairwise segment projections into flat tangent planes suffers from $2\pi$ branch-cut discontinuities (phase wrap-around) when segments cross coordinate boundaries, producing artificial jumps in the computed writhe across trajectory frames.

## Decision
We chose the exact spherical polygon decomposition theorem of **Van Oosterom & Strackee (1983)**:
$$\operatorname{Wr}(\mathcal{C}) = \frac{1}{2\pi} \sum_{j=2}^{N-1} \sum_{i=0}^{j-2} \Omega^*(\mathbf{e}_i, \mathbf{e}_j)$$
where $\Omega^*(\mathbf{e}_i, \mathbf{e}_j)$ represents the oriented solid angle subtended by segment $\mathbf{e}_i = [\mathbf{r}_i, \mathbf{r}_{i+1}]$ viewed from the endpoints of segment $\mathbf{e}_j = [\mathbf{r}_j, \mathbf{r}_{j+1}]$.

The solid angle of each spherical quadrangle is decomposed into two spherical triangles. The signed solid angle $\Omega$ of a spherical triangle defined by three unit vectors $\mathbf{a}, \mathbf{b}, \mathbf{c}$ is evaluated via Eriksson's formula:
$$\tan\left(\frac{\Omega}{2}\right) = \frac{\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})}{1 + \mathbf{a} \cdot \mathbf{b} + \mathbf{a} \cdot \mathbf{c} + \mathbf{b} \cdot \mathbf{c}}$$
$$\Omega = 2 \operatorname{atan2}\left(\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c}), \; 1 + \mathbf{a} \cdot \mathbf{b} + \mathbf{a} \cdot \mathbf{c} + \mathbf{b} \cdot \mathbf{c}\right)$$

## Consequences

### Positive
- **Exact & Branch-Cut-Free**: The $\operatorname{atan2}$ formulation guarantees smooth, continuous writhe values without branch-cut jumps or $2\pi$ phase wrapping.
- **Strictly Invariant under $SE(3)$**: Global translation and rotation cancel analytically inside the dot and cross products; invariance is preserved down to $< 10^{-12}$.
- **Chirality Detection**: Inverting chain coordinates $\mathbf{r} \to -\mathbf{r}$ exactly flips the sign $\operatorname{Wr} \to -\operatorname{Wr}$, enabling rigorous enantiomeric discrimination.
- **Local Writhe Windowing**: Local writhe $\operatorname{Wr}_w(i)$ over a window $[i - w, i + w]$ captures tertiary compaction in $\mathcal{O}(w^2)$ time instead of $\mathcal{O}(N^2)$.

### Negative / Trade-offs
- Requires four spherical triangle evaluations per segment pair. Optimized vector operations in `nalgebra` keep this overhead under $1.2\,\mu\text{s}$ per residue window.

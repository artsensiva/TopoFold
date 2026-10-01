# ADR-004: Pébay's $O(1)$ Single-Pass Streaming Central Moments over Full-Trajectory Sorting

## Status
Accepted (v0.6.5)

## Context
Autonomous detection of cryptic pockets and bistable functional loops requires analyzing the statistical distribution of local differential geometry invariants ($\kappa, \tau, \theta_\beta$) across millions of trajectory frames for sliding subcurve windows along the protein backbone.
To quantify whether a distribution is unimodal (rigid backbone, Gaussian thermal noise) or bimodal (two-state cryptic pocket opening, allosteric switch), TopoFold employs **Sarle's Bimodality Coefficient (BC)**:
$$BC = \frac{\gamma^2 + 1}{\kappa_{\text{kurt}} + 3 \cdot \frac{(n-1)^2}{(n-2)(n-3)}}$$
where $\gamma$ is sample skewness and $\kappa_{\text{kurt}}$ is sample excess kurtosis.

Evaluating skewness and kurtosis naively presents significant scalability bottlenecks:
1. **Memory Ingestion Overhead**: Storing full feature distributions for all sliding windows requires allocating $\mathcal{O}(W \cdot F)$ arrays in memory, causing severe cache thrashing and memory exhaustion on long trajectories ($F > 10^6$ frames).
2. **Multi-Pass / Sorting Algorithms**: Two-pass algorithms (pass 1: mean $\bar{x}$; pass 2: higher powers $\sum (x_i - \bar{x})^k$) or sorting-based non-parametric tests (Hartigan's Dip test) require storing the entire sequence and are non-streaming, preventing real-time online analysis during MD production runs.

## Decision
We implemented **Pébay's (2008) extension of Welford's algorithm** for streaming calculation of sample central moments ($M_1, M_2, M_3, M_4$) in `topofold-core::bimodality::StreamingMoments`.

For each new observation $x$, moments are updated online in $\mathcal{O}(1)$ time and $\mathcal{O}(1)$ space:
$$\delta = x - M_1, \quad \delta_n = \frac{\delta}{n}, \quad \delta_n^2 = \delta_n^2, \quad \text{term}_1 = \delta \cdot \delta_n \cdot (n - 1)$$
$$M_4 \leftarrow M_4 + \text{term}_1 \cdot \delta_n^2 \cdot (n^2 - 3n + 3) + 6 \delta_n^2 M_2 - 4 \delta_n M_3$$
$$M_3 \leftarrow M_3 + \text{term}_1 \cdot \delta_n \cdot (n - 2) - 3 \delta_n M_2$$
$$M_2 \leftarrow M_2 + \text{term}_1$$
$$M_1 \leftarrow M_1 + \delta_n$$

From these accumulators, central variance $m_2 = M_2 / n$, skewness $\gamma = (M_3 / n) / m_2^{1.5}$, and excess kurtosis $\kappa_{\text{kurt}} = (M_4 / n) / m_2^2 - 3.0$ are obtained in $\mathcal{O}(1)$ operations.

## Consequences

### Positive
- **$\mathcal{O}(1)$ Constant Memory**: Memory consumption is strictly 40 bytes per window accumulator (five `f64` numbers), independent of whether the trajectory contains $10^3$ or $10^9$ frames.
- **Single-Pass Streaming**: Trajectories can be streamed directly from disk (or directly from an active simulation socket) without storing past frames.
- **Numerical Stability**: Pébay's formulas avoid catastrophic floating-point cancellation inherent in textbook formulas $\sum x^2 - (\sum x)^2 / n$.
- **Ultra-High Throughput**: Scans an entire 2,500-frame BPTI trajectory in **$11.7\text{ ms}$** ($4.69\,\mu\text{s/frame}$).

### Negative / Trade-offs
- Small sample size sensitivity: For $n < 4$, skewness and kurtosis are mathematically undefined; the implementation explicitly guards against $n < 4$ by returning $BC = 0.0$.

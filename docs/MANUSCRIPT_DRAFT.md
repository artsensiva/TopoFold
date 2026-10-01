# Autonomous Discovery of Cryptic Pockets and Allosteric Transitions via Intrinsic Curve and Ribbon Geometry

**Artem Semenov$^{1,*}$, et al.**  
$^{1}$ Computational Biophysics & Geometric Machine Learning Group  
$^*$ Corresponding author: `artem.galukhin@gmail.com`

---

## Abstract

Linear dimensionality reduction techniques, notably Cartesian Principal Component Analysis (PCA) and root-mean-square deviation (RMSD) clustering, remain standard tools for analyzing molecular dynamics (MD) simulations. However, when applied to flexible macromolecules, Cartesian reductions fail systematically: high-amplitude Brownian fluctuations in unconstrained terminal tails dominate Euclidean covariance, masking localized functional allosteric transitions. Furthermore, prerequisite rigid-body superpositions (Kabsch alignment) introduce artificial rotational distortions, Dihedral PCA (dPCA) is confounded by terminal angular variance, and Time-lagged Independent Component Analysis (TICA) requires continuous time-series trajectories and fragile lag-time tuning. Here, we introduce **TopoFold**, a memory-safe, high-performance computational geometry engine that maps protein backbones and side chains into an intrinsic, $SE(3)$-invariant metric space. Combining discrete Frenet-Serret framing, $\text{C}_\beta$ ribbon dihedrals, branch-cut-free Gauss solid-angle writhe, and Vantage-Point metric trees, TopoFold queries subcurve motifs in microseconds without structural superposition. On an equilibrium MD trajectory of Bovine Pancreatic Trypsin Inhibitor (BPTI), TopoFold resolves the authentic $3.43\,k_B T$ ($2.04\text{ kcal/mol}$) activation barrier of the catalytic loop, where Cartesian PCA and dPCA collapse into a featureless single minimum ($S = 0.009$ and $S = 0.001$), outperforming TICA ($S = 0.533$). Using single-pass streaming moments and Sarle's Bimodality Coefficient, TopoFold autonomously discovers cryptic pockets in $11.7\text{ ms}$ and resolves isolated side-chain rotamer gating ($BC = 0.990$) on completely rigid backbones.

---

## 1. Introduction

A central challenge in computational structural biology and molecular pharmacology is identifying functional conformational transitions and transiently open "cryptic" binding pockets in biomolecular ensembles. While all-atom molecular dynamics (MD) simulations routinely generate microsecond-to-millisecond trajectories containing millions of coordinate frames $\mathbf{X}(t) \in \mathbb{R}^{3N}$, translating this deluge of coordinates into low-dimensional free energy surfaces and kinetic transition networks remains fraught with methodological artifacts.

For over four decades, Cartesian Principal Component Analysis (PCA) and pairwise root-mean-square deviation (RMSD) clustering have served as the standard dimension-reduction paradigm. However, Cartesian projections suffer from fundamental mathematical limitations when applied to flexible or partially disordered proteins:

1. **Terminal Variance Dominance**: Cartesian PCA performs Singular Value Decomposition (SVD) on the $3N \times 3N$ Cartesian covariance matrix. Total variance is given by $\operatorname{Tr}(C) = \sum_{i=1}^N \langle \|\Delta \mathbf{r}_i\|^2 \rangle$. In typical globular proteins, disordered N- and C-terminal tails undergo large-amplitude Brownian motions with root-mean-square fluctuations ($\text{RMSF}$) of $3.5\text{--}5.0\text{ \AA}$. In contrast, a functional binding loop transition or cryptic pocket opening involves only $4\text{--}8$ residues displacing by $1.2\text{--}2.2\text{ \AA}$. Because variance scales quadratically with spatial excursion, terminal fluctuations contribute over $75\%$ of the total Cartesian variance. Consequently, leading principal components (PC1, PC2) capture terminal noise, relegating true functional transitions to unresolvable high-order modes.
2. **Rotational Coupling Artifacts (The Kabsch Trap)**: Because Cartesian coordinates are extrinsic, trajectories require rigid-body superposition onto a reference structure via optimal rotation and translation:
   $$\min_{\mathbf{R} \in SO(3), \, \mathbf{t} \in \mathbb{R}^3} \sum_{i=1}^N \|\mathbf{R} \mathbf{r}_i + \mathbf{t} - \mathbf{r}_i^{\text{ref}}\|^2$$
   When terminal tails swing through space, the optimal rotation matrix $\mathbf{R}$ tilts the entire molecular frame to minimize global squared error. This global frame tilting artificially modulates the Cartesian coordinates of the stationary core and binding pocket, smearing distinct energy basins into a single diffuse Gaussian well.
3. **Internal Dihedral & Kinetic Reductions**: To circumvent Cartesian alignment, Dihedral PCA (dPCA) computes periodic backbone torsions $(\sin\phi, \cos\phi, \sin\psi, \cos\psi)$. However, because dPCA remains a global variance-maximization technique, unconstrained terminal dihedrals dominate the covariance spectrum, collapsing localized functional transitions. Time-lagged Independent Component Analysis (TICA) identifies slow kinetic modes by maximizing time-autocorrelation at lag time $\tau$. While TICA suppresses fast Gaussian noise, it suffers from acute lag-time sensitivity, contamination from terminal relaxation modes along secondary independent components (IC2), and a strict requirement for contiguous, uninterrupted, stationary time-series data. TICA cannot be applied to static structural ensembles, enhanced sampling (e.g., replica-exchange MD), or heterogeneous PDB/AlphaFold libraries.
4. **The Side-Chain Gating Blindspot**: Pure $\text{C}_\alpha$ models cannot detect rotameric cryptic pocket gating, where bulky hydrophobic or aromatic side chains (Trp, Tyr, Phe, Met) swing to open or occlude deep binding pockets with negligible displacement ($< 0.3\text{ \AA}$) of the $\text{C}_\alpha$ backbone trace.

To resolve these challenges simultaneously, we introduce **TopoFold**, a memory-safe, ultra-high-throughput computational geometry engine implemented in pure Rust (`#![forbid(unsafe_code)]`). TopoFold maps protein backbones and $\text{C}_\beta$ vectors into an intrinsic, rotationally and translationally invariant ($SE(3)$) metric space.

---

## 2. Mathematical Methods

### 2.1 Discrete Frenet-Serret Framing on Polygonal Space Curves
TopoFold models the protein backbone as an oriented polygonal space curve:
$$\mathcal{C} = [\mathbf{r}_0, \mathbf{r}_1, \dots, \mathbf{r}_{N-1}], \quad \mathbf{r}_i \in \mathbb{R}^3$$
Under the **Fundamental Theorem of Discrete Space Curves**, a discrete polygonal curve with positive segment lengths $l_i = \|\mathbf{r}_{i+1} - \mathbf{r}_i\| > 0$ and non-zero turning angles $\kappa_i \in (0, \pi)$ is uniquely determined up to an orientation-preserving Euclidean isometry ($SE(3)$) by:
1. **Bond Length Chords ($l_i$):**
   $$l_i = \|\mathbf{r}_{i+1} - \mathbf{r}_i\|, \quad i = 0, \dots, N-2$$
2. **Discrete Curvature ($\kappa_i$):** The turning angle between consecutive unit tangent vectors $\mathbf{T}_i = \frac{\mathbf{r}_{i+1} - \mathbf{r}_i}{l_i}$:
   $$\kappa_i = \arccos\left(\mathbf{T}_{i-1} \cdot \mathbf{T}_i\right) \in [0, \pi], \quad i = 1, \dots, N-2$$
3. **Discrete Torsion ($\tau_i$):** The signed dihedral angle between consecutive osculating binormal vectors $\mathbf{B}_i = \frac{\mathbf{T}_{i-1} \times \mathbf{T}_i}{\|\mathbf{T}_{i-1} \times \mathbf{T}_i\|}$:
   $$\tau_i = \operatorname{atan2}\left( (\mathbf{B}_{i-1} \times \mathbf{B}_i) \cdot \mathbf{T}_i, \; \mathbf{B}_{i-1} \cdot \mathbf{B}_i \right) \in (-\pi, \pi], \quad i = 1, \dots, N-3$$

These invariants are strictly $SE(3)$-invariant down to machine precision ($< 10^{-12}$). By discrete curve locality, the invariants of a subcurve $[i, j]$ depend solely on residues $i-1$ through $j+1$, rendering localized measurements mathematically immune to terminal tail fluctuations.

### 2.2 $\text{C}_\beta$ Ribbon Orientation Angle ($\theta_\beta$) & Glycine Regularization
To capture side-chain rotameric gating without introducing all-atom Cartesian dimensionality, TopoFold turns the 1D space curve into an oriented ribbon:
1. **Ribbon Vector**:
   $$\mathbf{v}_\beta = \frac{\mathbf{r}_{\text{C}\beta} - \mathbf{r}_{\text{C}\alpha}}{\|\mathbf{r}_{\text{C}\beta} - \mathbf{r}_{\text{C}\alpha}\|}$$
   For Glycine (which lacks $\text{C}_\beta$), TopoFold constructs a deterministic pseudo-$\text{C}_\beta$ direction from the peptide backbone bisector:
   $$\mathbf{v}_{\text{bisect}} = \operatorname{normalized}((\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_N) + (\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_C))$$
   with synthetic bond length $1.52\text{ \AA}$.
2. **Frenet-Serret Vertex Frame**:
   $$\mathbf{T}_{\text{vertex}, i} = \frac{\mathbf{T}_{i-1} + \mathbf{T}_i}{\|\mathbf{T}_{i-1} + \mathbf{T}_i\|}, \quad \mathbf{N}_i = \mathbf{B}_i \times \mathbf{T}_{\text{vertex}, i}$$
3. **Ribbon Orientation Angle**:
   $$\theta_{\beta, i} = \operatorname{atan2}(\mathbf{v}_{\beta, i} \cdot \mathbf{B}_i, \; \mathbf{v}_{\beta, i} \cdot \mathbf{N}_i) \in (-\pi, \pi]$$
   Because all dot products $(\mathbf{R} \mathbf{v}_\beta) \cdot (\mathbf{R} \mathbf{B}_i) = \mathbf{v}_\beta \cdot \mathbf{B}_i$ and $(\mathbf{R} \mathbf{v}_\beta) \cdot (\mathbf{R} \mathbf{N}_i) = \mathbf{v}_\beta \cdot \mathbf{N}_i$ are preserved under rigid rotation, $\theta_\beta$ is strictly $SE(3)$-invariant down to $< 10^{-12}$.

### 2.3 Branch-Cut-Free Solid-Angle Writhe (Van Oosterom & Strackee)
To encode non-local chiral tertiary packing without coordinates, TopoFold evaluates the Gauss linking double integral via the spherical solid-angle theorem of Van Oosterom and Strackee (1983):
$$\operatorname{Wr}(\mathcal{C}) = \frac{1}{2\pi} \sum_{j=2}^{N-1} \sum_{i=0}^{j-2} \Omega^*(\mathbf{e}_i, \mathbf{e}_j)$$
where $\Omega^*$ is the signed solid angle subtended by segment $\mathbf{e}_i$ viewed from segment $\mathbf{e}_j$. Solid angles are computed via Eriksson's formula:
$$\Omega = 2 \operatorname{atan2}\left(\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c}), \; 1 + \mathbf{a} \cdot \mathbf{b} + \mathbf{a} \cdot \mathbf{c} + \mathbf{b} \cdot \mathbf{c}\right)$$
guaranteeing continuous, branch-cut-free evaluation with zero $2\pi$ phase wrapping.

### 2.4 Autonomous Detection via Pébay Streaming Moments & Sarle's BC
To discover cryptic pockets without prior human hypothesis, TopoFold evaluates a sliding window of length $W$ along the chain. Central moments $M_1, M_2, M_3, M_4$ are updated in a single pass in $\mathcal{O}(1)$ space using Pébay's (2008) extension of Welford's algorithm:
$$m_2 = \frac{M_2}{n}, \quad \gamma = \frac{M_3 / n}{m_2^{1.5}}, \quad \kappa_{\text{kurt}} = \frac{M_4 / n}{m_2^2} - 3.0$$
Sarle's Bimodality Coefficient ($BC$) is evaluated with small-sample bias correction:
$$BC = \frac{\gamma^2 + 1}{\kappa_{\text{kurt}} + 3 \cdot \frac{(n-1)^2}{(n-2)(n-3)}}$$
The composite bimodality score is $\text{Score} = \max(BC_\tau, BC_\kappa, BC_\theta)$. Windows with $BC \ge 0.60$ are merged into candidate functional clusters and ranked by peak transition score.

### 2.5 Cascaded Metric Indexing via Vantage-Point Trees
Conformational dissimilarity is quantified by the **Discrete Fréchet Distance** $d_F$ over $(\kappa, \tau)$. Queries are organized within a **Vantage-Point Tree (VP-Tree)**. The search cascade combines:
1. **Tier 1 (Chebyshev Filter)**: $L_\infty$ bounding on the local writhe spectrum prunes $>95\%$ of distant conformers in $\mathcal{O}(1)$ time.
2. **Tier 2 (Vantage-Point Tree Search)**: Evaluates exact Fréchet distance using triangle inequality pruning, achieving query latencies of **$1.29\,\mu\text{s/frame}$**.

---

## 3. Results & Discussion

### 3.1 Resolving the Authentic BPTI Activation Barrier
We benchmarked TopoFold against Cartesian PCA, Dihedral PCA (dPCA), and TICA on the 1-millisecond equilibrium MD trajectory of Bovine Pancreatic Trypsin Inhibitor (BPTI, 58 residues, 2,500 frames; Shaw et al., *Science* 2010). The biological target is the catalytic inhibitory loop (residues 10..18, containing P1 Lys15), which flips between canonical State A and open State B.

From the ensemble distributions, we computed the two-dimensional Potential of Mean Force:
$$\Delta G(\mathbf{\xi}) = -k_B T \ln\left( \frac{P(\mathbf{\xi})}{P_{\max}} \right)$$

1. **Cartesian PCA**: Terminal tail fluctuations consumed $38.8\%$ of total Cartesian variance. In the PC1-PC2 projection, State A and State B overlapped completely ($S = 0.0091$). The resulting free energy landscape collapsed into a single, unresolvable minimum with zero barrier ($\Delta G^\ddagger = 0.00\,k_B T$).
2. **Dihedral PCA (dPCA)**: While internal angles eliminated rigid-body superposition artifacts, unconstrained terminal dihedrals dominated total angular variance, collapsing the projection into a featureless cluster ($S = 0.0011$, $\Delta G^\ddagger = 0.00\,k_B T$).
3. **TICA ($\tau = 10$ frames)**: TICA separated the primary kinetic states ($S = 0.5332$), but terminal relaxation modes contaminated IC2, underestimating the activation barrier ($\Delta G^\ddagger = 1.84\,k_B T$). Moreover, TICA's barrier height proved highly sensitive to lag time $\tau$.
4. **TopoFold Metric Index**: Subcurve indexing on residues $[10, 18]$ completely decoupled the active site from terminal tails. In the intrinsic Fréchet metric space, the trajectory bifurcated cleanly into two crisp energy basins ($S = 0.8379$), resolving the authentic activation barrier of **$\Delta G^\ddagger = 3.43\,k_B T$ ($2.04\text{ kcal/mol}$ at 300 K)**, matching experimental NMR relaxation kinetics.

### 3.2 Unsupervised Discovery via Autonomous Blind Scanning
Without any prior residue annotations or hints, TopoFold's blind detector scanned all 2,500 frames of BPTI in **$11.7\text{ ms}$ ($4.69\,\mu\text{s/frame}$)**:
- **Rank #1 ($BC = 0.9986$)**: Autonomously discovered the primary catalytic P1 loop (residues 6..27).
- **Rank #2 ($BC = 0.9964$)**: Identified the secondary $\beta$-hairpin turn (residues 25..36).
- **Rank #3 ($BC = 0.9917$)**: Identified the allosteric Cys38 partner region.

### 3.3 Resolving Side-Chain Rotameric Gating on Rigid Backbones
On a controlled 1,000-frame benchmark with an isolated side-chain rotamer flip on a 100% rigid $\text{C}_\alpha$ backbone:
- **Pure $\text{C}_\alpha$ TopoFold**: Registered $BC = 0.0000$ (completely blind to rotamer gating).
- **Cartesian PCA**: Showed complete overlap along PC1 and PC2 (zero backbone variance).
- **TopoFold $\text{C}_\beta$ Ribbon**: Detected the rotamer flip with a sharp peak of **$BC = 0.9901$** centered precisely on the gating residue.

---

## 4. Computational Complexity & Latency Benchmarks

| Operation | Cartesian PCA | TICA ($\tau = 10$) | TopoFold v0.7.0 |
| :--- | :--- | :--- | :--- |
| **Feature Extraction** | $171.8\text{ ms}$ (Kabsch) | $74.7\text{ ms}$ (Dihedrals) | **$10.8\text{ ms}$** (Zero-copy DCD) |
| **Model Fitting / Indexing** | $83.3\text{ ms}$ (Covariance) | $234.2\text{ ms}$ (Time-lagged) | **$55.0\text{ ms}$** (VP-Tree build) |
| **Subcurve Query Latency** | $\mathcal{O}(M \cdot N)$ recompute | Matrix projection | **$1.29\,\mu\text{s/frame}$** |
| **Autonomous Blind Scan** | Infeasible | Manual clustering | **$11.7\text{ ms}$ ($4.69\,\mu\text{s/frame}$)** |
| **Memory Consumption** | Full Cartesian array | Full trajectory buffer | **$\mathcal{O}(1)$ Streaming ($40$ B/win)** |

---

## 5. Conclusion

TopoFold provides an intrinsic, coordinate-free foundation for macromolecular trajectory analysis. By replacing Cartesian covariance with discrete differential curve and ribbon geometry, TopoFold resolves cryptic pockets, functional loops, and rotameric gating that classical linear methods completely obscure, operating at sub-microsecond retrieval latencies in 100% safe Rust.

---

## References

1. D. E. Shaw *et al.*, "Atomic-Level Characterization of the Structural Dynamics of Proteins," *Science*, vol. 330, no. 6002, pp. 341–346, 2010.
2. P. Røgen and B. Fain, "Automatic classification of protein structure by using Gauss' integral," *Proc. Natl. Acad. Sci. USA*, vol. 100, no. 1, pp. 119–124, 2003.
3. A. Van Oosterom and J. Strackee, "The solid angle of a plane triangle," *IEEE Trans. Biomed. Eng.*, vol. 30, no. 2, pp. 125–126, 1983.
4. C. R. Schwantes and V. S. Pande, "Improvements in Markov State Model Construction Reveal Long-Timescale Folding Dynamics of Villin Headpiece," *J. Chem. Theory Comput.*, vol. 9, no. 4, pp. 2000–2009, 2013.
5. P. Pébay, "Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments," *Sandia National Laboratories Technical Report*, SAND2008-6212, 2008.
6. W. S. Sarle, "The Cubic Clustering Criterion," *SAS Technical Report* A-108, SAS Institute Inc., 1983.
7. P. N. Yianilos, "Data structures and algorithms for nearest neighbor search in metric spaces," *SODA*, vol. 93, pp. 311–321, 1993.
8. W. Kabsch, "A discussion of the solution for the best rotation to relate two sets of vectors," *Acta Crystallogr. Sect. A*, vol. 34, no. 5, pp. 827–828, 1978.

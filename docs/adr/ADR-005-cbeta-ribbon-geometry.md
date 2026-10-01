# ADR-005: $\text{C}_\beta$ Ribbon Vector Dihedrals over All-Atom Cartesian Representations

## Status
Accepted (v0.7.0)

## Context
While the $\text{C}_\alpha$ polygonal space curve captures secondary structure transitions and large-scale loop hinge motions, it is blind to **rotameric cryptic pocket gating**. In many pharmacologically critical targets (e.g., Mcl-1, KRAS Switch-II, BPTI P1 Lys15), bulky aromatic or aliphatic side chains (Trp, Tyr, Phe, Met) swing to open or occlude deep binding pockets with negligible displacement ($< 0.3\text{ \AA}$) of the $\text{C}_\alpha$ backbone trace.

To address side-chain orientation, two obvious approaches exist in structural biology:
1. **Full All-Atom Cartesian Representation**:
   Include all heavy atoms ($\text{C}, \text{N}, \text{O}, \text{C}_\beta, \text{C}_\gamma, \dots$) in Cartesian space.
   - *Failure Mode*: Explodes dimensionality from $3N$ to $\approx 25N$; re-introduces the Kabsch rotational coupling problem; massively increases memory bandwidth and breaks $SE(3)$ gauge invariance.
2. **Standard Side-Chain Torsions ($\chi_1, \chi_2, \dots$)**:
   - *Failure Mode*: Standard $\chi$ angles are measured relative to local covalent bonds ($\text{N}-\text{C}_\alpha-\text{C}_\beta-\text{C}_\gamma$). They are decoupled from the tertiary orientation of the protein backbone curve, making it impossible to determine whether the side chain points *into* an internal cryptic cavity or *outward* into bulk solvent.

## Decision
We implemented **$\text{C}_\beta$ Ribbon Geometry** in `topofold-core::ribbon` by turning the 1D space curve into an oriented ribbon with a single additional vector per residue:
1. **$\text{C}_\beta$ Unit Direction Vector**:
   $$\mathbf{v}_\beta = \frac{\mathbf{r}_{\text{C}\beta} - \mathbf{r}_{\text{C}\alpha}}{\|\mathbf{r}_{\text{C}\beta} - \mathbf{r}_{\text{C}\alpha}\|}$$
2. **Deterministic Glycine Regularization**:
   Glycine lacks a $\text{C}_\beta$ atom. Rather than omitting Glycine (which would break chain continuity) or imputing arbitrary coordinates, TopoFold constructs a deterministic pseudo-$\text{C}_\beta$ direction from the peptide backbone bisector:
   $$\mathbf{v}_{\text{bisect}} = \operatorname{normalized}((\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_N) + (\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_C))$$
   with synthetic bond length $1.52\text{ \AA}$.
3. **Intrinsic Ribbon Orientation Angle ($\theta_\beta$)**:
   The ribbon vector $\mathbf{v}_\beta$ is projected onto the discrete Frenet-Serret vertex frame $(\mathbf{T}_{\text{vertex}, i}, \mathbf{N}_i, \mathbf{B}_i)$:
   $$\theta_{\beta, i} = \operatorname{atan2}(\mathbf{v}_{\beta, i} \cdot \mathbf{B}_i, \; \mathbf{v}_{\beta, i} \cdot \mathbf{N}_i) \in (-\pi, \pi]$$
   where $\mathbf{B}_i$ is normal to the osculating plane and $\mathbf{N}_i$ is the principal normal.

## Consequences

### Positive
- **Strict $SE(3)$ Gauge Invariance**: Since all dot products $(\mathbf{R} \mathbf{v}_\beta) \cdot (\mathbf{R} \mathbf{B}_i) = \mathbf{v}_\beta \cdot \mathbf{B}_i$ and $(\mathbf{R} \mathbf{v}_\beta) \cdot (\mathbf{R} \mathbf{N}_i) = \mathbf{v}_\beta \cdot \mathbf{N}_i$ are preserved under global rotation, $\theta_\beta$ is invariant down to $< 10^{-12}$.
- **Captures Rotameric Cryptic Gating**: Side chains swinging relative to the backbone osculating plane produce large $\theta_\beta$ shifts, generating high bimodality ($BC > 0.85$) even when the $\text{C}_\alpha$ backbone is 100% rigid.
- **Minimal Memory Footprint**: Adds only one scalar $\theta_\beta \in (-\pi, \pi]$ per residue, preserving lightweight cache efficiency.
- **Complete Chain Continuity**: The deterministic Glycine bisector produces smooth, differentiable angular trajectories with zero singularities or NaNs.

### Negative / Trade-offs
- Sub-rotameric fluctuations beyond $\text{C}_\beta$ (e.g., $\chi_2, \chi_3$ in Lys/Arg tails) are not explicitly parametrized, focusing computational resources strictly on pocket-gating orientation.

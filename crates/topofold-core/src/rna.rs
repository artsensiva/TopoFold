//! Ribonucleic Ribbon Geometry and Invariants for RNA 3D Dynamics.
//!
//! Extends TopoFold's discrete differential geometry engine to ribonucleic acids (RNA),
//! resolving tertiary folding, pseudoknots, and riboswitch conformational switching hinges.
//!
//! # Biophysical Motivation
//!
//! RNA backbones possess six rotatable dihedral angles per nucleotide ($\alpha, \beta, \gamma, \delta, \epsilon, \zeta$)
//! alongside ribose sugar pucker (C2'-endo / C3'-endo), creating extreme conformational flexibility.
//! CASP-RNA benchmarks demonstrated that deep learning structural predictors systematically struggle
//! to predict RNA tertiary dynamics and allosteric switching transitions between ligand-bound
//! (compact pseudoknot) and ligand-free (open expression platform) states.
//!
//! TopoFold formulates SE(3)-invariant **Ribonucleic Ribbon Geometry**:
//! 1. **Backbone Space Curve**: Modeled along Phosphorus ($\text{P}$) atoms.
//! 2. **Glycosidic Base Orientation Vector**:
//!    $$\mathbf{v}_{\text{base}} = \frac{\mathbf{r}_N - \mathbf{r}_{C1'}}{\|\mathbf{r}_N - \mathbf{r}_{C1'}\|}$$
//!    pointing from ribose $\text{C1'}$ to the glycosidic nitrogen ($\text{N9}$ for purines A/G, $\text{N1}$ for pyrimidines C/U).
//! 3. **Discrete Frenet-Serret Framing & Ribbon Invariants**:
//!    - Curvature $\kappa_{P, i} \in [0, \pi]$ (turning angles between consecutive $\text{P}-\text{P}$ segments).
//!    - Torsion $\tau_{P, i} \in (-\pi, \pi]$ (dihedral angle between consecutive osculating planes).
//!    - Local Writhe $\text{Wr}_{P, i}$ (discrete Gauss linking integral measuring tertiary non-planarity).
//!    - Glycosidic Ribbon Dihedral Angle $\theta_{\text{base}, i} \in (-\pi, \pi]$ (orientation of the nucleobase around the curve tangent).
//!
//! Under any rigid-body rotation $\mathbf{R} \in SO(3)$ and translation $\mathbf{t} \in \mathbb{R}^3$,
//! all four invariants are preserved down to $< 10^{-12}$.

#![forbid(unsafe_code)]

use nalgebra::{Point3, UnitVector3, Vector3};

use crate::bimodality::StreamingMoments;
use crate::error::GeometryError;
use crate::geometry::GEOMETRY_EPSILON;
use crate::types::BackboneTrace;

/// An RNA ribonucleic ribbon trace tracking Phosphorus backbone, ribose C4', C1',
/// and glycosidic nitrogen atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct RnaRibbonTrace {
    /// 3D coordinates of Phosphorus atoms (in Ångströms).
    pub(crate) p_coords: Vec<Point3<f64>>,
    /// 3D coordinates of ribose C4' atoms (in Ångströms).
    pub(crate) c4_coords: Vec<Point3<f64>>,
    /// 3D coordinates of ribose C1' atoms (in Ångströms).
    pub(crate) c1_coords: Vec<Point3<f64>>,
    /// 3D coordinates of glycosidic nitrogen atoms (N9 for purines, N1 for pyrimidines, in Ångströms).
    pub(crate) n_coords: Vec<Point3<f64>>,
    /// Nucleotide 1-letter or 3-letter names (e.g. "A", "G", "C", "U").
    pub(crate) names: Vec<String>,
    /// Residue sequence numbers (PDB resSeq).
    pub(crate) seq_ids: Vec<i32>,
    /// Chain identifiers.
    pub(crate) chain_ids: Vec<char>,
}

impl RnaRibbonTrace {
    /// Creates a new `RnaRibbonTrace` from matching atomic coordinate sets.
    ///
    /// # Errors
    /// Returns [`GeometryError::InsufficientPoints`] if fewer than 2 nucleotides are provided,
    /// or if coordinate slice lengths do not match.
    pub fn new(
        p_coords: Vec<Point3<f64>>,
        c4_coords: Vec<Point3<f64>>,
        c1_coords: Vec<Point3<f64>>,
        n_coords: Vec<Point3<f64>>,
    ) -> Result<Self, GeometryError> {
        let n = p_coords.len();
        if c4_coords.len() != n || c1_coords.len() != n || n_coords.len() != n {
            return Err(GeometryError::InsufficientPoints {
                required: n,
                actual: c4_coords.len().min(c1_coords.len()).min(n_coords.len()),
            });
        }
        if n < 2 {
            return Err(GeometryError::InsufficientPoints {
                required: 2,
                actual: n,
            });
        }
        let names = vec!["N".to_string(); n];
        let seq_ids = (1..=n as i32).collect();
        let chain_ids = vec!['A'; n];

        Ok(Self {
            p_coords,
            c4_coords,
            c1_coords,
            n_coords,
            names,
            seq_ids,
            chain_ids,
        })
    }

    /// Creates an `RnaRibbonTrace` with full nucleotide metadata.
    ///
    /// # Errors
    /// Returns [`GeometryError`] if coordinate or metadata lengths do not match.
    pub fn with_meta(
        p_coords: Vec<Point3<f64>>,
        c4_coords: Vec<Point3<f64>>,
        c1_coords: Vec<Point3<f64>>,
        n_coords: Vec<Point3<f64>>,
        names: Vec<String>,
        seq_ids: Vec<i32>,
        chain_ids: Vec<char>,
    ) -> Result<Self, GeometryError> {
        let n = p_coords.len();
        if c4_coords.len() != n
            || c1_coords.len() != n
            || n_coords.len() != n
            || names.len() != n
            || seq_ids.len() != n
            || chain_ids.len() != n
        {
            return Err(GeometryError::InsufficientPoints {
                required: n,
                actual: c4_coords.len(),
            });
        }
        if n < 2 {
            return Err(GeometryError::InsufficientPoints {
                required: 2,
                actual: n,
            });
        }
        Ok(Self {
            p_coords,
            c4_coords,
            c1_coords,
            n_coords,
            names,
            seq_ids,
            chain_ids,
        })
    }

    /// Creates an `RnaRibbonTrace` from coordinate slices `&[[f64; 3]]`.
    ///
    /// # Errors
    /// Returns [`GeometryError`] if lengths do not match or are insufficient.
    pub fn from_arrays(
        p: &[[f64; 3]],
        c4: &[[f64; 3]],
        c1: &[[f64; 3]],
        n: &[[f64; 3]],
    ) -> Result<Self, GeometryError> {
        let p_pts: Vec<Point3<f64>> = p.iter().map(|&[x, y, z]| Point3::new(x, y, z)).collect();
        let c4_pts: Vec<Point3<f64>> = c4.iter().map(|&[x, y, z]| Point3::new(x, y, z)).collect();
        let c1_pts: Vec<Point3<f64>> = c1.iter().map(|&[x, y, z]| Point3::new(x, y, z)).collect();
        let n_pts: Vec<Point3<f64>> = n.iter().map(|&[x, y, z]| Point3::new(x, y, z)).collect();
        Self::new(p_pts, c4_pts, c1_pts, n_pts)
    }

    /// Returns the number of nucleotides in the RNA trace.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.p_coords.len()
    }

    /// Returns `true` if the trace contains zero nucleotides.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.p_coords.is_empty()
    }

    /// Slice access to Phosphorus coordinates.
    #[inline]
    #[must_use]
    pub fn p_coords(&self) -> &[Point3<f64>] {
        &self.p_coords
    }

    /// Slice access to ribose C4' coordinates.
    #[inline]
    #[must_use]
    pub fn c4_coords(&self) -> &[Point3<f64>] {
        &self.c4_coords
    }

    /// Slice access to ribose C1' coordinates.
    #[inline]
    #[must_use]
    pub fn c1_coords(&self) -> &[Point3<f64>] {
        &self.c1_coords
    }

    /// Slice access to glycosidic nitrogen coordinates.
    #[inline]
    #[must_use]
    pub fn n_coords(&self) -> &[Point3<f64>] {
        &self.n_coords
    }

    /// Slice access to nucleotide names.
    #[inline]
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Slice access to sequence residue numbers.
    #[inline]
    #[must_use]
    pub fn seq_ids(&self) -> &[i32] {
        &self.seq_ids
    }

    /// Slice access to chain identifiers.
    #[inline]
    #[must_use]
    pub fn chain_ids(&self) -> &[char] {
        &self.chain_ids
    }

    /// Converts Phosphorus coordinates into a standard [`BackboneTrace`].
    #[must_use]
    pub fn to_p_backbone(&self) -> BackboneTrace {
        BackboneTrace::new(self.p_coords.clone())
    }

    /// Computes the unit vectors pointing from C1' to the glycosidic nitrogen
    /// for each nucleotide:
    ///
    /// $$\mathbf{v}_{\text{base}, i} = \frac{\mathbf{r}_{N, i} - \mathbf{r}_{C1', i}}{\|\mathbf{r}_{N, i} - \mathbf{r}_{C1', i}\|}$$
    ///
    /// # Errors
    /// Returns [`GeometryError::DegenerateSegment`] if any C1' to N distance is $< \text{GEOMETRY\_EPSILON}$.
    pub fn base_unit_vectors(&self) -> Result<Vec<UnitVector3<f64>>, GeometryError> {
        let mut vectors = Vec::with_capacity(self.len());
        for i in 0..self.len() {
            let diff: Vector3<f64> = self.n_coords[i] - self.c1_coords[i];
            let len = diff.norm();
            if len < GEOMETRY_EPSILON {
                return Err(GeometryError::DegenerateSegment {
                    index: i,
                    next_index: i,
                    length: len,
                });
            }
            vectors.push(UnitVector3::new_unchecked(diff / len));
        }
        Ok(vectors)
    }
}

/// Discrete geometric invariants of an RNA ribonucleic ribbon space curve.
#[derive(Debug, Clone, PartialEq)]
pub struct RnaCurveInvariants {
    /// Virtual bond lengths between consecutive Phosphorus atoms: $\|P_{i+1} - P_i\|$, length $N - 1$.
    pub segment_lengths: Vec<f64>,
    /// Discrete curvature of the Phosphorus backbone: turning angle $\kappa_{P, i} \in [0, \pi]$, length $N - 2$.
    pub curvatures: Vec<f64>,
    /// Discrete torsion of the Phosphorus backbone: dihedral angle $\tau_{P, i} \in (-\pi, \pi]$, length $N - 3$.
    pub torsions: Vec<f64>,
    /// Glycosidic base ribbon orientation dihedral angles $\theta_{\text{base}, i} \in (-\pi, \pi]$, length $N - 2$.
    pub base_dihedrals: Vec<f64>,
    /// Local writhe spectrum along Phosphorus backbone, length $N - 2 \times \text{window\_radius}$ (or empty).
    pub local_writhes: Vec<f64>,
}

impl RnaCurveInvariants {
    /// Total integrated curvature along the Phosphorus backbone trace.
    #[must_use]
    pub fn total_curvature(&self) -> f64 {
        self.curvatures.iter().sum()
    }

    /// Average Phosphorus-Phosphorus virtual bond length.
    #[must_use]
    pub fn mean_segment_length(&self) -> f64 {
        if self.segment_lengths.is_empty() {
            0.0
        } else {
            self.segment_lengths.iter().sum::<f64>() / self.segment_lengths.len() as f64
        }
    }
}

/// Computes the glycosidic base ribbon orientation dihedral angles $\theta_{\text{base}}$
/// for all internal nucleotides $i \in \{1, \dots, N-2\}$.
///
/// At internal Phosphorus vertex $i$:
/// - Evaluates Frenet-Serret frame $(\mathbf{T}_{\text{vertex}, i}, \mathbf{N}_i, \mathbf{B}_i)$ from Phosphorus curve.
/// - Evaluates base unit vector $\mathbf{v}_{\text{base}, i} = \frac{\mathbf{r}_{N, i} - \mathbf{r}_{C1', i}}{\|\mathbf{r}_{N, i} - \mathbf{r}_{C1', i}\|}$.
/// - Returns signed dihedral angle around $\mathbf{T}_{\text{vertex}, i}$:
///   $$\theta_{\text{base}, i} = \operatorname{atan2}(\mathbf{v}_{\text{base}, i} \cdot \mathbf{B}_i, \; \mathbf{v}_{\text{base}, i} \cdot \mathbf{N}_i) \in (-\pi, \pi]$$
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if $N < 3$.
/// Returns [`GeometryError::DegenerateSegment`] if any Phosphorus bond or base vector is degenerate.
/// Returns [`GeometryError::CollinearSegments`] if any consecutive tangents are collinear.
pub fn compute_base_ribbon_dihedrals(trace: &RnaRibbonTrace) -> Result<Vec<f64>, GeometryError> {
    let n = trace.len();
    if n < 3 {
        return Err(GeometryError::InsufficientPoints {
            required: 3,
            actual: n,
        });
    }

    let p_coords = trace.p_coords();

    // 1. Tangents along Phosphorus backbone (length N - 1)
    let mut tangents = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let diff: Vector3<f64> = p_coords[i + 1] - p_coords[i];
        let len = diff.norm();
        if len < GEOMETRY_EPSILON {
            return Err(GeometryError::DegenerateSegment {
                index: i,
                next_index: i + 1,
                length: len,
            });
        }
        tangents.push(UnitVector3::new_unchecked(diff / len));
    }

    // 2. Base unit vectors (length N)
    let base_vecs = trace.base_unit_vectors()?;

    // 3. Ribbon dihedrals at each internal vertex i = 1 .. N - 2
    let mut dihedrals = Vec::with_capacity(n - 2);
    for i in 1..n - 1 {
        let t_prev = &tangents[i - 1];
        let t_curr = &tangents[i];

        // Binormal B_i = normalized(T_{i-1} x T_i)
        let cross = t_prev.cross(t_curr);
        let cross_norm = cross.norm();
        if cross_norm < GEOMETRY_EPSILON {
            let dot = t_prev.dot(t_curr).clamp(-1.0, 1.0);
            return Err(GeometryError::CollinearSegments {
                vertex_index: i,
                turning_angle: dot.acos(),
            });
        }
        let b_i = UnitVector3::new_unchecked(cross / cross_norm);

        // Vertex tangent T_vertex = normalized(T_{i-1} + T_i)
        let t_sum = t_prev.as_ref() + t_curr.as_ref();
        let t_sum_norm = t_sum.norm();
        let t_vertex = if t_sum_norm >= GEOMETRY_EPSILON {
            UnitVector3::new_unchecked(t_sum / t_sum_norm)
        } else {
            *t_curr
        };

        // In-plane principal normal N_i = B_i x T_vertex
        let n_i = UnitVector3::new_unchecked(b_i.cross(&t_vertex));

        let v_base = &base_vecs[i];
        let proj_normal = v_base.dot(&n_i);
        let proj_binormal = v_base.dot(&b_i);

        let theta = proj_binormal.atan2(proj_normal);
        dihedrals.push(theta);
    }

    Ok(dihedrals)
}

/// Extracts complete SE(3)-invariant ribonucleic curve and ribbon invariants for an [`RnaRibbonTrace`].
///
/// Uses default local writhe window radius $w = 2$.
///
/// # Errors
/// Returns [`GeometryError`] if trace has fewer than 4 nucleotides or degenerate geometry.
pub fn extract_rna_ribbon_invariants(
    trace: &RnaRibbonTrace,
) -> Result<RnaCurveInvariants, GeometryError> {
    extract_rna_ribbon_invariants_with_window(trace, 2)
}

/// Extracts complete SE(3)-invariant ribonucleic curve and ribbon invariants with a custom
/// writhe sliding window radius.
///
/// # Errors
/// Returns [`GeometryError`] if trace has fewer than 4 nucleotides or degenerate geometry.
pub fn extract_rna_ribbon_invariants_with_window(
    trace: &RnaRibbonTrace,
    writhe_window_radius: usize,
) -> Result<RnaCurveInvariants, GeometryError> {
    if trace.len() < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: trace.len(),
        });
    }

    let p_bb = trace.to_p_backbone();
    let (segment_lengths, tangents) = crate::geometry::compute_tangents(&p_bb)?;
    let curvatures = crate::geometry::compute_curvatures(&tangents)?;
    let torsions = crate::geometry::compute_torsions(&tangents)?;
    let base_dihedrals = compute_base_ribbon_dihedrals(trace)?;

    let local_writhes = if trace.len() >= 2 * writhe_window_radius + 2 {
        crate::writhe::compute_local_writhe(&p_bb, writhe_window_radius)?
    } else {
        Vec::new()
    };

    Ok(RnaCurveInvariants {
        segment_lengths,
        curvatures,
        torsions,
        base_dihedrals,
        local_writhes,
    })
}

/// A detected candidate RNA conformational switching hinge or pseudoknot junction.
#[derive(Debug, Clone, PartialEq)]
pub struct RnaHingeCandidate {
    /// Starting nucleotide index of the detected candidate cluster (0-based, inclusive).
    pub start_res: usize,
    /// Ending nucleotide index of the detected candidate cluster (0-based, inclusive).
    pub end_res: usize,
    /// Peak composite bimodality score within this candidate cluster.
    pub score: f64,
    /// Sarle's bimodality coefficient for discrete torsion at the peak window.
    pub bc_tau: f64,
    /// Sarle's bimodality coefficient for discrete curvature at the peak window.
    pub bc_kappa: f64,
    /// Sarle's bimodality coefficient for glycosidic base ribbon dihedral at the peak window.
    pub bc_theta: f64,
    /// Nucleotide index of the peak window within the cluster.
    pub peak_res: usize,
}

/// Sequence-wide bimodality profile tuple: `(scores, bc_theta, bc_kappa, bc_tau)`.
pub type RnaBimodalityProfile = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

/// Computes the sequence-wide bimodality profile for an ensemble of RNA ribbon traces.
///
/// Evaluates Sarle's Bimodality Coefficient across sliding windows of size `window_size`
/// for composite score, base ribbon dihedral $\theta_{\text{base}}$, curvature $\kappa_P$, and torsion $\tau_P$.
///
/// # Returns
/// A 4-tuple of vectors: `(scores, bc_theta, bc_kappa, bc_tau)`.
///
/// # Errors
/// Returns [`GeometryError`] if trajectory has fewer than 4 frames or mismatched lengths.
pub fn compute_rna_bimodality_profile(
    traces: &[RnaRibbonTrace],
    window_size: usize,
) -> Result<RnaBimodalityProfile, GeometryError> {
    if traces.len() < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: traces.len(),
        });
    }

    let n_nt = traces[0].len();
    if n_nt < window_size + 2 {
        return Err(GeometryError::InsufficientPoints {
            required: window_size + 2,
            actual: n_nt,
        });
    }

    // Extract invariants for all frames
    let mut ensemble_invariants = Vec::with_capacity(traces.len());
    for trace in traces {
        if trace.len() != n_nt {
            return Err(GeometryError::InsufficientPoints {
                required: n_nt,
                actual: trace.len(),
            });
        }
        ensemble_invariants.push(extract_rna_ribbon_invariants(trace)?);
    }

    // Evaluated internal vertices count: n_nt - 2 (curvatures and base dihedrals)
    let n_internal = n_nt - 2;
    if n_internal < window_size {
        return Ok((Vec::new(), Vec::new(), Vec::new(), Vec::new()));
    }

    let n_windows = n_internal - window_size + 1;
    let mut scores = Vec::with_capacity(n_windows);
    let mut bc_theta = Vec::with_capacity(n_windows);
    let mut bc_kappa = Vec::with_capacity(n_windows);
    let mut bc_tau = Vec::with_capacity(n_windows);

    let ref_inv = &ensemble_invariants[0];

    for w_idx in 0..n_windows {
        let mut moments_kappa = StreamingMoments::new();
        let mut moments_tau = StreamingMoments::new();
        let mut moments_theta = StreamingMoments::new();

        for inv in &ensemble_invariants {
            // Mean curvature in window (curvature in [0, pi] has no branch cut)
            let k_slice = &inv.curvatures[w_idx..w_idx + window_size];
            let mean_k = k_slice.iter().sum::<f64>() / window_size as f64;
            moments_kappa.update(mean_k);

            // Mean unwrapped base dihedral in window
            let mut th_sum = 0.0;
            for j in 0..window_size {
                let th = inv.base_dihedrals[w_idx + j];
                let th_ref = ref_inv.base_dihedrals[w_idx + j];
                let diff = (th - th_ref).sin().atan2((th - th_ref).cos());
                th_sum += th_ref + diff;
            }
            let mean_th = th_sum / window_size as f64;
            moments_theta.update(mean_th);

            // Mean unwrapped torsion in window (clamped to available torsions)
            let tau_start = w_idx.min(inv.torsions.len().saturating_sub(1));
            let tau_end = (w_idx + window_size).min(inv.torsions.len());
            if tau_end > tau_start {
                let mut tau_sum = 0.0;
                for j in tau_start..tau_end {
                    let tau = inv.torsions[j];
                    let tau_ref = ref_inv.torsions[j];
                    let diff = (tau - tau_ref).sin().atan2((tau - tau_ref).cos());
                    tau_sum += tau_ref + diff;
                }
                let mean_t = tau_sum / (tau_end - tau_start) as f64;
                moments_tau.update(mean_t);
            }
        }

        let bc_k = moments_kappa.bimodality_coefficient();
        let bc_t = moments_tau.bimodality_coefficient();
        let bc_th = moments_theta.bimodality_coefficient();
        let score = bc_k.max(bc_t).max(bc_th);

        scores.push(score);
        bc_theta.push(bc_th);
        bc_kappa.push(bc_k);
        bc_tau.push(bc_t);
    }

    Ok((scores, bc_theta, bc_kappa, bc_tau))
}

/// Scans an ensemble of RNA ribbon traces across trajectory frames to detect bistable
/// conformational switching hinges (such as the regulatory P1 stem junction in riboswitches).
///
/// Evaluates Sarle's Bimodality Coefficient over sliding windows of size `window_size`
/// across frames for curvature $\kappa_P$, torsion $\tau_P$, and base ribbon dihedral $\theta_{\text{base}}$.
///
/// # Errors
/// Returns [`GeometryError`] if trajectory has fewer than 4 frames or mismatched lengths.
pub fn scan_rna_switching_hinges(
    traces: &[RnaRibbonTrace],
    window_size: usize,
    threshold: f64,
) -> Result<Vec<RnaHingeCandidate>, GeometryError> {
    let (scores, bc_theta, bc_kappa, bc_tau) = compute_rna_bimodality_profile(traces, window_size)?;

    // Cluster contiguous windows with score >= threshold
    let mut candidates = Vec::new();
    let mut in_cluster = false;
    let mut cluster_start = 0;
    let mut cluster_end = 0;
    let mut peak_score = 0.0;
    let mut peak_res = 0;
    let mut peak_bc_t = 0.0;
    let mut peak_bc_k = 0.0;
    let mut peak_bc_th = 0.0;

    for (w_idx, (((&score, &bc_th), &bc_k), &bc_t)) in scores
        .iter()
        .zip(bc_theta.iter())
        .zip(bc_kappa.iter())
        .zip(bc_tau.iter())
        .enumerate()
    {
        if score >= threshold {
            if !in_cluster {
                in_cluster = true;
                cluster_start = w_idx + 1; // +1 to convert internal vertex to 0-based nt index
                cluster_end = w_idx + window_size;
                peak_score = score;
                peak_res = w_idx + 1;
                peak_bc_t = bc_t;
                peak_bc_k = bc_k;
                peak_bc_th = bc_th;
            } else {
                cluster_end = w_idx + window_size;
                if score > peak_score {
                    peak_score = score;
                    peak_res = w_idx + 1;
                    peak_bc_t = bc_t;
                    peak_bc_k = bc_k;
                    peak_bc_th = bc_th;
                }
            }
        } else if in_cluster {
            in_cluster = false;
            candidates.push(RnaHingeCandidate {
                start_res: cluster_start,
                end_res: cluster_end,
                score: peak_score,
                bc_tau: peak_bc_t,
                bc_kappa: peak_bc_k,
                bc_theta: peak_bc_th,
                peak_res,
            });
        }
    }

    if in_cluster {
        candidates.push(RnaHingeCandidate {
            start_res: cluster_start,
            end_res: cluster_end,
            score: peak_score,
            bc_tau: peak_bc_t,
            bc_kappa: peak_bc_k,
            bc_theta: peak_bc_th,
            peak_res,
        });
    }

    // Rank candidates by peak score descending
    candidates.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(candidates)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Isometry3, Translation3, UnitQuaternion};

    #[test]
    fn test_se3_strict_invariance_rna() {
        // Construct an idealized RNA helix segment (8 nucleotides)
        let n = 8;
        let mut p_coords = Vec::with_capacity(n);
        let mut c4_coords = Vec::with_capacity(n);
        let mut c1_coords = Vec::with_capacity(n);
        let mut n_coords = Vec::with_capacity(n);

        let radius_p = 8.8; // Å
        let pitch_per_nt = 2.8; // Å
        let angle_step = 0.55; // rad (~32.7 deg, standard A-form RNA helix)

        for i in 0..n {
            let theta = i as f64 * angle_step;
            let z = i as f64 * pitch_per_nt;

            let p = Point3::new(radius_p * theta.cos(), radius_p * theta.sin(), z);
            let c4 = Point3::new(
                (radius_p - 2.5) * (theta + 0.1).cos(),
                (radius_p - 2.5) * (theta + 0.1).sin(),
                z + 1.2,
            );
            let c1 = Point3::new(
                (radius_p - 4.5) * (theta + 0.2).cos(),
                (radius_p - 4.5) * (theta + 0.2).sin(),
                z + 1.8,
            );
            let n_pt = Point3::new(
                (radius_p - 6.0) * (theta + 0.3).cos(),
                (radius_p - 6.0) * (theta + 0.3).sin(),
                z + 2.0,
            );

            p_coords.push(p);
            c4_coords.push(c4);
            c1_coords.push(c1);
            n_coords.push(n_pt);
        }

        let trace = RnaRibbonTrace::new(
            p_coords.clone(),
            c4_coords.clone(),
            c1_coords.clone(),
            n_coords.clone(),
        )
        .unwrap();
        let inv_original = extract_rna_ribbon_invariants(&trace).unwrap();

        // Apply arbitrary 3D rigid rotation and translation (SE(3))
        let rotation = UnitQuaternion::from_euler_angles(0.785, -1.234, 2.456);
        let translation = Translation3::new(42.5, -99.3, 137.8);
        let isometry = Isometry3::from_parts(translation, rotation);

        let p_trans: Vec<Point3<f64>> = p_coords
            .iter()
            .map(|&p| isometry.transform_point(&p))
            .collect();
        let c4_trans: Vec<Point3<f64>> = c4_coords
            .iter()
            .map(|&p| isometry.transform_point(&p))
            .collect();
        let c1_trans: Vec<Point3<f64>> = c1_coords
            .iter()
            .map(|&p| isometry.transform_point(&p))
            .collect();
        let n_trans: Vec<Point3<f64>> = n_coords
            .iter()
            .map(|&p| isometry.transform_point(&p))
            .collect();

        let trace_trans = RnaRibbonTrace::new(p_trans, c4_trans, c1_trans, n_trans).unwrap();
        let inv_transformed = extract_rna_ribbon_invariants(&trace_trans).unwrap();

        // Verify strict SE(3) invariance down to < 10^-12
        for i in 0..inv_original.segment_lengths.len() {
            assert!(
                (inv_original.segment_lengths[i] - inv_transformed.segment_lengths[i]).abs()
                    < 1e-12,
                "Segment length mismatch at {}: {} vs {}",
                i,
                inv_original.segment_lengths[i],
                inv_transformed.segment_lengths[i]
            );
        }

        for i in 0..inv_original.curvatures.len() {
            assert!(
                (inv_original.curvatures[i] - inv_transformed.curvatures[i]).abs() < 1e-12,
                "Curvature mismatch at {}: {} vs {}",
                i,
                inv_original.curvatures[i],
                inv_transformed.curvatures[i]
            );
        }

        for i in 0..inv_original.torsions.len() {
            assert!(
                (inv_original.torsions[i] - inv_transformed.torsions[i]).abs() < 1e-12,
                "Torsion mismatch at {}: {} vs {}",
                i,
                inv_original.torsions[i],
                inv_transformed.torsions[i]
            );
        }

        for i in 0..inv_original.base_dihedrals.len() {
            assert!(
                (inv_original.base_dihedrals[i] - inv_transformed.base_dihedrals[i]).abs() < 1e-12,
                "Base dihedral mismatch at {}: {} vs {}",
                i,
                inv_original.base_dihedrals[i],
                inv_transformed.base_dihedrals[i]
            );
        }

        for i in 0..inv_original.local_writhes.len() {
            assert!(
                (inv_original.local_writhes[i] - inv_transformed.local_writhes[i]).abs() < 1e-12,
                "Local writhe mismatch at {}: {} vs {}",
                i,
                inv_original.local_writhes[i],
                inv_transformed.local_writhes[i]
            );
        }
    }

    #[test]
    fn test_rna_switching_hinge_detection() {
        // Construct an ensemble of 100 frames across 16 nucleotides with a bistable switching hinge
        // at nucleotides 5..8 (simulating riboswitch stem unpairing/bending)
        let n_frames = 100;
        let n_nt = 16;
        let mut ensemble = Vec::with_capacity(n_frames);

        for f in 0..n_frames {
            let is_state_b = f >= 50;
            let mut p_coords = Vec::with_capacity(n_nt);
            let mut c4_coords = Vec::with_capacity(n_nt);
            let mut c1_coords = Vec::with_capacity(n_nt);
            let mut n_coords = Vec::with_capacity(n_nt);

            for i in 0..n_nt {
                let hinge_offset = if is_state_b && (5..=8).contains(&i) {
                    1.4 // Distinct conformation in State B
                } else {
                    0.0
                };

                let theta = i as f64 * 0.55 + hinge_offset;
                let z = i as f64 * 2.8;

                let p = Point3::new(8.8 * theta.cos(), 8.8 * theta.sin(), z);
                let c4 = Point3::new(
                    6.3 * (theta + 0.1).cos(),
                    6.3 * (theta + 0.1).sin(),
                    z + 1.2,
                );
                let c1 = Point3::new(
                    4.3 * (theta + 0.2).cos(),
                    4.3 * (theta + 0.2).sin(),
                    z + 1.8,
                );
                let n_pt = Point3::new(
                    2.8 * (theta + 0.3).cos(),
                    2.8 * (theta + 0.3).sin(),
                    z + 2.0,
                );

                p_coords.push(p);
                c4_coords.push(c4);
                c1_coords.push(c1);
                n_coords.push(n_pt);
            }

            ensemble.push(RnaRibbonTrace::new(p_coords, c4_coords, c1_coords, n_coords).unwrap());
        }

        let hinges = scan_rna_switching_hinges(&ensemble, 3, 0.70).unwrap();
        assert!(
            !hinges.is_empty(),
            "Expected candidate switching hinges to be detected"
        );
        let top_hinge = &hinges[0];
        assert!(
            top_hinge.score > 0.85,
            "Expected top hinge score > 0.85, got {}",
            top_hinge.score
        );
        // Hinge cluster should encompass the modified hinge residues 5..8
        assert!(
            top_hinge.start_res <= 5 && top_hinge.end_res >= 7,
            "Hinge cluster [{}, {}] does not overlap expected hinge 5..8",
            top_hinge.start_res,
            top_hinge.end_res
        );
    }
}

//! Spectral Topological Density and Transient Motif Detection for Disordered Ensembles (IDPs).
//!
//! # Biophysical Motivation
//!
//! Intrinsically Disordered Proteins (IDPs) and intrinsically disordered regions (IDRs)—such
//! as human Alpha-Synuclein (140 residues), Tau, and p53 transactivation domain—lack a fixed,
//! static tertiary structure. Machine-learning structure prediction tools like AlphaFold output
//! arbitrary low-confidence conformations (pLDDT < 50) that fail to capture the true thermodynamic
//! ensemble.
//!
//! In Cartesian space, conformational ensembles of IDPs exhibit enormous RMSD fluctuations
//! (> 15 Å), causing Cartesian PCA and structural alignment methods to collapse into uninformative
//! isotropic Gaussian clouds (Silhouette score ~ 0).
//!
//! However, pathological self-assembly (e.g. amyloid aggregation in Parkinson's disease) is
//! driven by *transient pre-structured motifs*—metastable conformations where localized segments
//! (such as the hydrophobic non-amyloid component NACore, residues 61..95 in Alpha-Synuclein)
//! intermittently nucleate high-curvature turns, $\beta$-hairpins, or transient helical turns.
//!
//! # Mathematical Foundations
//!
//! TopoFold formulates the **Spectral Topological Density** by tracking the localized coupling
//! between solid-angle writhe ($\text{Wr}$) and discrete curvature ($\kappa$) along the space curve:
//!
//! ## 1. Localized Topological Compactness
//!
//! For frame $f$ and residue $c \in [0, N-1]$ with sliding window radius $R$:
//! - Let the subchain window be $[w_L, w_R]$ where $w_L = \max(0, c - R)$ and $w_R = \min(N - 1, c + R)$.
//! - The local discrete writhe $\text{Wr}_{\text{local}}(c, f)$ is computed via the Gauss linking integral
//!   over all non-adjacent segment pairs within the subchain window:
//!   $$\text{Wr}_{\text{local}}(c, f) = \sum_{i = w_L}^{w_R - 2} \sum_{j = i + 2}^{w_R - 1} w(\mathbf{r}_i, \mathbf{r}_{i+1}, \mathbf{r}_j, \mathbf{r}_{j+1})$$
//! - The local discrete curvature $\bar{\kappa}_{\text{local}}(c, f)$ is the mean turning angle
//!   across internal vertices in $[w_L, w_R]$.
//! - The instantaneous topological compactness is the product:
//!   $$\mathcal{C}(c, f) = |\text{Wr}_{\text{local}}(c, f)| \cdot \bar{\kappa}_{\text{local}}(c, f)$$
//!
//! Extended random coils and Brownian tails have negligible local writhe ($|\text{Wr}| \to 0$),
//! suppressing false-positive thermal fluctuations. Only when a segment simultaneously achieves
//! non-planar chiral self-wrapping and high backbone curvature does $\mathcal{C}(c, f)$ peak.
//!
//! ## 2. Ensemble Spectral Topological Density
//!
//! Across an ensemble of $F$ conformations:
//! - **Mean Topological Density**:
//!   $$S_{\text{topo}}(c) = \mathbb{E}_F[\mathcal{C}(c, f)] = \frac{1}{F} \sum_{f=1}^F \mathcal{C}(c, f)$$
//! - **Ensemble Variance**:
//!   $$\sigma^2_{\text{topo}}(c) = \frac{1}{F - 1} \sum_{f=1}^F (\mathcal{C}(c, f) - S_{\text{topo}}(c))^2$$
//! - **Sequence-Normalized Z-score**:
//!   $$Z(c) = \frac{S_{\text{topo}}(c) - \mu_S}{\max(\sigma_S, 10^{-12})}$$
//!   where $\mu_S$ and $\sigma_S$ are the mean and standard deviation of $S_{\text{topo}}$ across the sequence.
//!
//! ## 3. Autonomous Transient Motif Discovery
//!
//! Contiguous or near-contiguous regions with $Z(c) \ge Z_{\text{threshold}}$ (default: 1.0)
//! are grouped into [`IdpMotif`] candidate islands, discovering transiently ordered nucleation hubs.

#![forbid(unsafe_code)]

use rayon::prelude::*;

use crate::error::GeometryError;
use crate::geometry::compute_tangents;
use crate::ribbon::RibbonTrace;
use crate::types::BackboneTrace;
use crate::writhe::segment_pair_writhe;

/// A detected transiently structured pre-nucleation motif in an intrinsically disordered ensemble.
#[derive(Debug, Clone, PartialEq)]
pub struct IdpMotif {
    /// Starting residue index (0-based, inclusive).
    pub start_res: usize,
    /// Ending residue index (0-based, inclusive).
    pub end_res: usize,
    /// Residue index of the peak topological compactness within this motif.
    pub peak_res: usize,
    /// Mean topological compactness across the motif residues.
    pub mean_compactness: f64,
    /// Peak topological compactness within this motif.
    pub peak_compactness: f64,
    /// Peak sequence Z-score of this motif.
    pub z_score: f64,
}

/// Complete Spectral Topological Density profile across all residues of an IDP ensemble.
#[derive(Debug, Clone, PartialEq)]
pub struct IdpDensityProfile {
    /// Mean topological compactness spectrum $S_{\text{topo}}(i)$ across all frames.
    pub mean_density: Vec<f64>,
    /// Ensemble variance of topological compactness across frames $\sigma^2_{\text{topo}}(i)$.
    pub variance_density: Vec<f64>,
    /// Ensemble standard deviation of topological compactness $\sigma_{\text{topo}}(i)$.
    pub std_density: Vec<f64>,
    /// Sequence-normalized Z-scores of topological compactness $Z(i)$.
    pub z_scores: Vec<f64>,
}

/// Computes the instantaneous localized topological compactness $\mathcal{C}(i)$ for each residue
/// of a single backbone conformation.
///
/// Output length: `trace.len()`.
///
/// # Arguments
/// - `trace`: The C-alpha backbone trace.
/// - `window_radius`: Half-width of the sliding window in residues (e.g. 4 for a 9-residue window).
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `trace.len() < 4`.
pub fn compute_frame_topological_compactness(
    trace: &BackboneTrace,
    window_radius: usize,
) -> Result<Vec<f64>, GeometryError> {
    let coords = trace.coordinates();
    let n = coords.len();

    if n < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: n,
        });
    }

    // Precompute tangents and curvature for each vertex
    let (_, tangents) = compute_tangents(trace)?;
    let mut vertex_kappas = Vec::with_capacity(n);
    vertex_kappas.push(0.0); // vertex 0 (boundary)

    for k in 1..n - 1 {
        let dot = tangents[k - 1].dot(&tangents[k]).clamp(-1.0, 1.0);
        vertex_kappas.push(dot.acos());
    }
    vertex_kappas.push(0.0); // vertex n - 1 (boundary)

    let mut compactness = vec![0.0; n];

    for c in 0..n {
        let start = c.saturating_sub(window_radius);
        let end = (c + window_radius).min(n - 1);
        let window_len = end - start + 1;

        if window_len < 4 {
            compactness[c] = 0.0;
            continue;
        }

        // 1. Local discrete writhe of subsegment [start ..= end]
        let mut local_wr = 0.0;
        for i in start..end - 1 {
            let p_i = &coords[i];
            let p_i1 = &coords[i + 1];
            for j in (i + 2)..end {
                let p_j = &coords[j];
                let p_j1 = &coords[j + 1];
                local_wr += segment_pair_writhe(p_i, p_i1, p_j, p_j1);
            }
        }

        // 2. Local mean curvature across vertices in the window
        let k_start = start + 1;
        let k_end = end; // exclusive (internal vertices k_start .. k_end)
        let local_kappa = if k_end > k_start {
            let sum_kappa: f64 = vertex_kappas[k_start..k_end].iter().sum();
            sum_kappa / ((k_end - k_start) as f64)
        } else {
            vertex_kappas[c]
        };

        // 3. Compactness product: |Wr_local| * kappa_local
        compactness[c] = local_wr.abs() * local_kappa;
    }

    Ok(compactness)
}

/// Computes the complete Spectral Topological Density profile across an ensemble of backbone traces.
///
/// Evaluates:
/// - Mean Topological Density: $S_{\text{topo}}(i) = \mathbb{E}_F[\mathcal{C}(i, f)]$
/// - Ensemble Variance: $\sigma^2_{\text{topo}}(i) = \operatorname{Var}_F[\mathcal{C}(i, f)]$
/// - Sequence-wide Z-scores: $Z(i) = (S_{\text{topo}}(i) - \mu_S) / \sigma_S$
///
/// # Arguments
/// - `traces`: Ensemble of backbone traces across trajectory frames.
/// - `window_radius`: Half-width of sliding subchain window (e.g. 4 for 9-residue subchains).
///
/// # Errors
/// Returns [`GeometryError`] if traces have mismatched lengths or fewer than 4 residues.
pub fn compute_idp_density_profile(
    traces: &[BackboneTrace],
    window_radius: usize,
) -> Result<IdpDensityProfile, GeometryError> {
    if traces.is_empty() {
        return Err(GeometryError::InsufficientPoints {
            required: 1,
            actual: 0,
        });
    }

    let n = traces[0].len();
    if n < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: n,
        });
    }

    for (idx, trace) in traces.iter().enumerate() {
        if trace.len() != n {
            return Err(GeometryError::InsufficientPoints {
                required: n,
                actual: trace.len(),
            });
        }
        if idx == 0 && trace.len() < 4 {
            return Err(GeometryError::InsufficientPoints {
                required: 4,
                actual: trace.len(),
            });
        }
    }

    // Parallel evaluation of compactness across all frames via Rayon
    let frame_compactness: Vec<Vec<f64>> = traces
        .par_iter()
        .map(|trace| compute_frame_topological_compactness(trace, window_radius))
        .collect::<Result<Vec<_>, _>>()?;

    let f_count = frame_compactness.len();
    let mut mean_density = vec![0.0; n];
    let mut variance_density = vec![0.0; n];
    let mut std_density = vec![0.0; n];

    // Compute mean
    for frame in &frame_compactness {
        for (i, &val) in frame.iter().enumerate() {
            mean_density[i] += val;
        }
    }
    for m in &mut mean_density {
        *m /= f_count as f64;
    }

    // Compute sample variance and standard deviation
    if f_count > 1 {
        let denom = (f_count - 1) as f64;
        for frame in &frame_compactness {
            for (i, &val) in frame.iter().enumerate() {
                let diff = val - mean_density[i];
                variance_density[i] += diff * diff;
            }
        }
        for (var, std) in variance_density.iter_mut().zip(std_density.iter_mut()) {
            *var /= denom;
            *std = var.sqrt();
        }
    }

    // Sequence-wide Z-score
    let sum_mean: f64 = mean_density.iter().sum();
    let seq_mean = sum_mean / (n as f64);
    let seq_var: f64 = mean_density
        .iter()
        .map(|&x| (x - seq_mean).powi(2))
        .sum::<f64>()
        / ((n - 1).max(1) as f64);
    let seq_std = seq_var.sqrt().max(1e-12);

    let z_scores: Vec<f64> = mean_density
        .iter()
        .map(|&x| (x - seq_mean) / seq_std)
        .collect();

    Ok(IdpDensityProfile {
        mean_density,
        variance_density,
        std_density,
        z_scores,
    })
}

/// Computes the Spectral Topological Density profile for an ensemble of [`RibbonTrace`] structures.
///
/// # Errors
/// Returns [`GeometryError`] if coordinate traces are invalid.
pub fn compute_idp_density_profile_ribbon(
    traces: &[RibbonTrace],
    window_radius: usize,
) -> Result<IdpDensityProfile, GeometryError> {
    let bb_traces: Vec<BackboneTrace> = traces.iter().map(|r| r.to_backbone_trace()).collect();
    compute_idp_density_profile(&bb_traces, window_radius)
}

/// Detects statistically significant transiently structured motifs within an ensemble of [`RibbonTrace`] structures.
///
/// Uses an autonomous clustering of sequence positions with $Z(i) \ge 1.0$ and bridges small gaps ($\le 2$ residues).
///
/// # Arguments
/// - `traces`: Ensemble of ribbon traces.
/// - `window_size`: Width of sliding subchain window (e.g. 8 or 9 residues).
///
/// # Errors
/// Returns [`GeometryError`] if trajectory data is invalid.
pub fn detect_transient_motifs(
    traces: &[RibbonTrace],
    window_size: usize,
) -> Result<Vec<IdpMotif>, GeometryError> {
    let bb_traces: Vec<BackboneTrace> = traces.iter().map(|r| r.to_backbone_trace()).collect();
    detect_transient_motifs_backbone(&bb_traces, window_size)
}

/// Detects statistically significant transiently structured motifs within an ensemble of [`BackboneTrace`] structures.
///
/// # Arguments
/// - `traces`: Ensemble of backbone traces.
/// - `window_size`: Width of sliding subchain window (e.g. 8 or 9 residues).
///
/// # Errors
/// Returns [`GeometryError`] if trajectory data is invalid.
pub fn detect_transient_motifs_backbone(
    traces: &[BackboneTrace],
    window_size: usize,
) -> Result<Vec<IdpMotif>, GeometryError> {
    let window_radius = (window_size / 2).max(1);
    detect_transient_motifs_with_params(traces, window_radius, 1.0, 2)
}

/// Detects transiently structured motifs with configurable parameters.
///
/// # Arguments
/// - `traces`: Ensemble of backbone traces.
/// - `window_radius`: Half-width of sliding window.
/// - `z_threshold`: Minimum sequence Z-score to include a residue in a motif (default: 1.0).
/// - `max_gap`: Maximum allowed residue gap between adjacent residues in the same cluster (default: 2).
///
/// # Errors
/// Returns [`GeometryError`] if trajectory data is invalid.
pub fn detect_transient_motifs_with_params(
    traces: &[BackboneTrace],
    window_radius: usize,
    z_threshold: f64,
    max_gap: usize,
) -> Result<Vec<IdpMotif>, GeometryError> {
    let profile = compute_idp_density_profile(traces, window_radius)?;
    let n = profile.mean_density.len();

    let above_threshold: Vec<usize> = (0..n)
        .filter(|&i| profile.z_scores[i] >= z_threshold)
        .collect();

    if above_threshold.is_empty() {
        return Ok(Vec::new());
    }

    // Cluster consecutive or near-consecutive residues bridging small gaps
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut current_cluster = vec![above_threshold[0]];

    for &idx in &above_threshold[1..] {
        let last_idx = *current_cluster.last().unwrap();
        if idx - last_idx <= max_gap + 1 {
            current_cluster.push(idx);
        } else {
            clusters.push(current_cluster);
            current_cluster = vec![idx];
        }
    }
    if !current_cluster.is_empty() {
        clusters.push(current_cluster);
    }

    let mut motifs = Vec::with_capacity(clusters.len());

    for cluster in clusters {
        let start_res = cluster[0];
        let end_res = *cluster.last().unwrap();

        let mut peak_val = 0.0;
        let mut peak_res = start_res;
        let mut peak_z = 0.0;
        let mut sum_val = 0.0;

        for &res_idx in &cluster {
            let val = profile.mean_density[res_idx];
            let z = profile.z_scores[res_idx];
            sum_val += val;
            if val > peak_val {
                peak_val = val;
                peak_res = res_idx;
                peak_z = z;
            }
        }

        let mean_compactness = sum_val / (cluster.len() as f64);

        motifs.push(IdpMotif {
            start_res,
            end_res,
            peak_res,
            mean_compactness,
            peak_compactness: peak_val,
            z_score: peak_z,
        });
    }

    // Sort by peak compactness in descending order
    motifs.sort_by(|a, b| {
        b.peak_compactness
            .partial_cmp(&a.peak_compactness)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(motifs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    #[test]
    fn test_straight_line_zero_topological_compactness() {
        // A straight line of 20 points
        let pts: Vec<Point3<f64>> = (0..20)
            .map(|i| Point3::new(i as f64 * 3.8, 0.0, 0.0))
            .collect();
        let trace = BackboneTrace::new(pts);

        let compactness = compute_frame_topological_compactness(&trace, 4).unwrap();
        assert_eq!(compactness.len(), 20);
        for (i, &val) in compactness.iter().enumerate() {
            assert!(
                val.abs() < 1e-10,
                "Residue {} in straight line should have 0 compactness, got {}",
                i,
                val
            );
        }
    }

    #[test]
    fn test_planar_curve_zero_topological_compactness() {
        // A planar zigzag/circle in the XY plane: local writhe must be zero
        let pts: Vec<Point3<f64>> = (0..20)
            .map(|i| {
                let theta = (i as f64) * std::f64::consts::PI / 10.0;
                Point3::new(10.0 * theta.cos(), 10.0 * theta.sin(), 0.0)
            })
            .collect();
        let trace = BackboneTrace::new(pts);

        let compactness = compute_frame_topological_compactness(&trace, 4).unwrap();
        assert_eq!(compactness.len(), 20);
        for (i, &val) in compactness.iter().enumerate() {
            assert!(
                val.abs() < 1e-10,
                "Planar curve residue {} should have 0 compactness due to planar writhe=0, got {}",
                i,
                val
            );
        }
    }

    #[test]
    fn test_helix_motif_detection() {
        // Build a synthetic 40-residue chain:
        // Residues 0..15: straight coil
        // Residues 15..28: chiral 3D helix (high writhe, high curvature)
        // Residues 28..40: straight coil
        let mut pts = Vec::with_capacity(40);
        for i in 0..15 {
            pts.push(Point3::new(i as f64 * 3.8, 0.0, 0.0));
        }
        for i in 15..28 {
            let t = (i - 15) as f64;
            let angle = t * 1.5;
            let x = 15.0 * 3.8 + t * 1.5;
            let y = 3.0 * angle.cos();
            let z = 3.0 * angle.sin();
            pts.push(Point3::new(x, y, z));
        }
        for i in 28..40 {
            let t = (i - 28) as f64;
            let x = 15.0 * 3.8 + 13.0 * 1.5 + t * 3.8;
            pts.push(Point3::new(x, 0.0, 0.0));
        }

        let trace = BackboneTrace::new(pts);
        let traces = vec![trace.clone(), trace];

        let profile = compute_idp_density_profile(&traces, 4).unwrap();
        assert_eq!(profile.mean_density.len(), 40);

        // Residues in the straight tails (0..10 and 32..40) must have 0 compactness
        assert!(profile.mean_density[2] < 1e-6);
        assert!(profile.mean_density[38] < 1e-6);

        // Helical region must have significant topological density
        let max_helix_compactness = profile.mean_density[15..28]
            .iter()
            .cloned()
            .fold(0.0, f64::max);
        assert!(
            max_helix_compactness > 0.01,
            "Helical segment should have positive topological compactness, got {}",
            max_helix_compactness
        );

        let motifs = detect_transient_motifs_backbone(&traces, 8).unwrap();
        assert!(!motifs.is_empty(), "Should detect at least 1 motif");
        let top = &motifs[0];
        assert!(
            top.start_res <= 20 && top.end_res >= 20,
            "Top motif should span helical region around residue 20, got {}..{}",
            top.start_res,
            top.end_res
        );
    }
}

//! Intrinsic Allosteric Networks via Mutual Information on Differential Curve Invariants.
//!
//! Maps conformational trajectories into an $N \times N$ communication adjacency matrix
//! by computing the Generalized Correlation Coefficient $r_{\text{MI}}$ from the
//! Mutual Information $I(\mathbf{V}_i; \mathbf{V}_j)$ of discrete curve invariants:
//! $\mathbf{V}_i(t) = [\kappa_i(t), \tau_i(t), \theta_{\beta, i}(t)] \in \mathbb{R}^3$.
//!
//! Under a Gaussian copula approximation, the Mutual Information between two multivariate
//! invariant random vectors is:
//! $$I(\mathbf{V}_i; \mathbf{V}_j) = -\frac{1}{2} \ln \left( \frac{\det \mathbf{C}_{ij}}{\det \mathbf{C}_i \det \mathbf{C}_j} \right)$$
//! where $\mathbf{C}_i, \mathbf{C}_j \in \mathbb{R}^{d \times d}$ are the sample covariance matrices,
//! and $\mathbf{C}_{ij} \in \mathbb{R}^{2d \times 2d}$ is the joint covariance matrix regularized
//! with diagonal jitter $\epsilon \mathbf{I}$.
//!
//! The generalized correlation coefficient:
//! $$r_{\text{MI}}(i, j) = \sqrt{1 - \exp(-2 I(\mathbf{V}_i; \mathbf{V}_j))} \in [0, 1]$$
//! provides a strictly $SE(3)$-invariant, superposition-free allosteric network metric.

#![forbid(unsafe_code)]

use nalgebra::DMatrix;
use ndarray::Array2;
use rayon::prelude::*;

use crate::error::GeometryError;
use crate::geometry::extract_curve_invariants;
use crate::ribbon::RibbonTrace;
use crate::types::BackboneTrace;

/// Computes the log determinant of a symmetric positive-definite matrix using Cholesky decomposition.
/// Falls back gracefully to regularized determinant if Cholesky fails.
fn log_det_symm(m: &DMatrix<f64>) -> f64 {
    let dim = m.nrows();
    if let Some(chol) = m.clone().cholesky() {
        let l = chol.l();
        let mut sum_ln = 0.0;
        let mut valid = true;
        for k in 0..dim {
            let diag = l[(k, k)];
            if diag > 0.0 {
                sum_ln += diag.ln();
            } else {
                valid = false;
                break;
            }
        }
        if valid && sum_ln.is_finite() {
            return 2.0 * sum_ln;
        }
    }

    // Fallback: LU decomposition or jittered determinant
    let det = m.determinant();
    if det > 1e-30 {
        det.ln()
    } else {
        let jitter = DMatrix::<f64>::identity(dim, dim) * 1e-4;
        let alt = m + jitter;
        let alt_det = alt.determinant();
        if alt_det > 1e-30 {
            alt_det.ln()
        } else {
            0.0
        }
    }
}

/// Single-residue covariance matrix and regularized log-determinant.
struct SingleResidueCov {
    cov: DMatrix<f64>,
    log_det: f64,
}

/// Extracted and centered differential curve invariants for a conformational chain ensemble.
struct CenteredChainFeatures {
    n_residues: usize,
    n_internal: usize,
    feat_dim: usize,
    centered: Vec<Vec<Vec<f64>>>, // [res_internal][frame][feat]
    single_covs: Vec<SingleResidueCov>,
}

/// Extracts centered differential curve invariant features for a single chain across all frames.
fn extract_centered_features(
    traces: &[BackboneTrace],
    regularizer_eps: f64,
) -> Result<CenteredChainFeatures, GeometryError> {
    let f_count = traces.len();
    if f_count < 3 {
        return Err(GeometryError::InsufficientPoints {
            required: 3,
            actual: f_count,
        });
    }

    let n_residues = traces[0].len();
    if n_residues < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: n_residues,
        });
    }

    for trace in traces {
        if trace.len() != n_residues {
            return Err(GeometryError::InsufficientPoints {
                required: n_residues,
                actual: trace.len(),
            });
        }
    }

    let has_cbeta = traces[0].cb_coordinates().is_some();
    let feat_dim = if has_cbeta { 3 } else { 2 };

    // Step 1: Extract invariants for all frames in parallel
    let invariants: Vec<_> = traces
        .par_iter()
        .map(extract_curve_invariants)
        .collect::<Result<Vec<_>, _>>()?;

    let n_internal = n_residues - 2;

    // Raw feature extraction: [residue_internal_idx][frame_idx][dim]
    let mut raw_features = vec![vec![vec![0.0; feat_dim]; f_count]; n_internal];

    for (f, inv) in invariants.iter().enumerate() {
        let n_curv = inv.curvatures.len();
        let n_tors = inv.torsions.len();
        let has_sc = !inv.sidechain_dihedrals.is_empty();

        for (i, feat_i) in raw_features.iter_mut().enumerate() {
            let k = if i < n_curv { inv.curvatures[i] } else { 0.0 };
            let t = if i < n_tors {
                inv.torsions[i]
            } else if n_tors > 0 {
                inv.torsions[n_tors - 1]
            } else {
                0.0
            };

            feat_i[f][0] = k;
            feat_i[f][1] = t;

            if has_cbeta {
                let th = if has_sc && i < inv.sidechain_dihedrals.len() {
                    inv.sidechain_dihedrals[i]
                } else {
                    0.0
                };
                feat_i[f][2] = th;
            }
        }
    }

    // Step 2: Mean-center features with circular wrapping on torsion and ribbon angles
    let mut centered = vec![vec![vec![0.0; feat_dim]; f_count]; n_internal];

    for i in 0..n_internal {
        let mean_k: f64 = raw_features[i].iter().map(|v| v[0]).sum::<f64>() / f_count as f64;

        let sum_sin_t: f64 = raw_features[i].iter().map(|v| v[1].sin()).sum();
        let sum_cos_t: f64 = raw_features[i].iter().map(|v| v[1].cos()).sum();
        let mean_t = sum_sin_t.atan2(sum_cos_t);

        let mean_th = if has_cbeta {
            let sum_sin_th: f64 = raw_features[i].iter().map(|v| v[2].sin()).sum();
            let sum_cos_th: f64 = raw_features[i].iter().map(|v| v[2].cos()).sum();
            sum_sin_th.atan2(sum_cos_th)
        } else {
            0.0
        };

        for f in 0..f_count {
            centered[i][f][0] = raw_features[i][f][0] - mean_k;

            let raw_t = raw_features[i][f][1];
            let diff_t = (raw_t - mean_t).sin().atan2((raw_t - mean_t).cos());
            centered[i][f][1] = diff_t;

            if has_cbeta {
                let raw_th = raw_features[i][f][2];
                let diff_th = (raw_th - mean_th).sin().atan2((raw_th - mean_th).cos());
                centered[i][f][2] = diff_th;
            }
        }
    }

    // Step 3: Precompute single-residue covariance matrices and their regularized log determinants
    let eps = if regularizer_eps > 0.0 { regularizer_eps } else { 1e-6 };
    let inv_df = 1.0 / (f_count - 1) as f64;

    let single_covs: Vec<SingleResidueCov> = (0..n_internal)
        .into_par_iter()
        .map(|i| {
            let mut cov = DMatrix::<f64>::zeros(feat_dim, feat_dim);
            for u in &centered[i] {
                for r in 0..feat_dim {
                    for c in 0..feat_dim {
                        cov[(r, c)] += u[r] * u[c];
                    }
                }
            }
            cov *= inv_df;
            for d in 0..feat_dim {
                cov[(d, d)] += eps;
            }
            let log_det = log_det_symm(&cov);
            SingleResidueCov { cov, log_det }
        })
        .collect();

    Ok(CenteredChainFeatures {
        n_residues,
        n_internal,
        feat_dim,
        centered,
        single_covs,
    })
}

/// Computes the $N \times N$ Intrinsic Allosteric Network Matrix across a trajectory of conformers.
///
/// Parameters
/// ----------
/// traces : &[BackboneTrace]
///     Ensemble of $F$ conformational frames, each with $N$ residues.
/// regularizer_eps : f64
///     Diagonal jitter $\epsilon$ added to covariance matrices (default: $10^{-6}$).
///
/// Returns
/// -------
/// Array2<f64>
///     Symmetric $N \times N$ matrix of generalized correlation coefficients $r_{\text{MI}} \in [0, 1]$.
pub fn compute_intrinsic_allosteric_network(
    traces: &[BackboneTrace],
    regularizer_eps: f64,
) -> Result<Array2<f64>, GeometryError> {
    let f_count = traces.len();
    let feats = extract_centered_features(traces, regularizer_eps)?;

    let n_residues = feats.n_residues;
    let n_internal = feats.n_internal;
    let feat_dim = feats.feat_dim;
    let joint_dim = feat_dim * 2;
    let inv_df = 1.0 / (f_count - 1) as f64;
    let eps = if regularizer_eps > 0.0 { regularizer_eps } else { 1e-6 };

    // Parallelize outer loop over rows i in 0..n_internal
    let row_results: Vec<Vec<(usize, f64)>> = (0..n_internal)
        .into_par_iter()
        .map(|i| {
            let mut row = Vec::with_capacity(n_internal - i);
            let cov_i = &feats.single_covs[i].cov;
            let log_det_i = feats.single_covs[i].log_det;

            for j in (i + 1)..n_internal {
                let cov_j = &feats.single_covs[j].cov;
                let log_det_j = feats.single_covs[j].log_det;

                // Compute cross-covariance C_ij_cross of shape (feat_dim, feat_dim)
                let mut cross = DMatrix::<f64>::zeros(feat_dim, feat_dim);
                for f in 0..f_count {
                    let u_i = &feats.centered[i][f];
                    let u_j = &feats.centered[j][f];
                    for r in 0..feat_dim {
                        for c in 0..feat_dim {
                            cross[(r, c)] += u_i[r] * u_j[c];
                        }
                    }
                }
                cross *= inv_df;

                // Assemble 2d x 2d joint covariance matrix
                let mut joint_cov = DMatrix::<f64>::zeros(joint_dim, joint_dim);
                // Top-left: cov_i
                for r in 0..feat_dim {
                    for c in 0..feat_dim {
                        joint_cov[(r, c)] = cov_i[(r, c)];
                    }
                }
                // Bottom-right: cov_j
                for r in 0..feat_dim {
                    for c in 0..feat_dim {
                        joint_cov[(r + feat_dim, c + feat_dim)] = cov_j[(r, c)];
                    }
                }
                // Top-right: cross
                for r in 0..feat_dim {
                    for c in 0..feat_dim {
                        joint_cov[(r, c + feat_dim)] = cross[(r, c)];
                    }
                }
                // Bottom-left: cross.T
                for r in 0..feat_dim {
                    for c in 0..feat_dim {
                        joint_cov[(r + feat_dim, c)] = cross[(c, r)];
                    }
                }

                // Add regularizer to joint covariance
                for d in 0..joint_dim {
                    joint_cov[(d, d)] += eps;
                }

                let log_det_joint = log_det_symm(&joint_cov);
                let mi = -0.5 * (log_det_joint - log_det_i - log_det_j);
                let clamped_mi = if mi > 0.0 { mi } else { 0.0 };

                // Generalized Correlation Coefficient r_MI = sqrt(1 - exp(-2 * MI))
                let r_mi = (1.0 - (-2.0 * clamped_mi).exp()).clamp(0.0, 1.0).sqrt();
                row.push((j, r_mi));
            }
            row
        })
        .collect();

    // Step 5: Fill into N x N symmetric output array
    let mut network = Array2::<f64>::zeros((n_residues, n_residues));

    // Diagonal elements: r_MI(i, i) = 1.0
    for i in 0..n_residues {
        network[[i, i]] = 1.0;
    }

    // Fill internal residues (residues 1 .. n_residues - 2)
    for (i_idx, row) in row_results.into_iter().enumerate() {
        let res_i = i_idx + 1;
        for (j_idx, r_mi) in row {
            let res_j = j_idx + 1;
            network[[res_i, res_j]] = r_mi;
            network[[res_j, res_i]] = r_mi;
        }
    }

    Ok(network)
}

/// Computes the $N_A \times N_B$ Inter-Molecular Allosteric Network Matrix between two protein chains.
///
/// Parameters
/// ----------
/// chain_a : &[BackboneTrace]
///     Ensemble of $F$ conformational frames for Chain A ($N_A$ residues).
/// chain_b : &[BackboneTrace]
///     Ensemble of $F$ conformational frames for Chain B ($N_B$ residues).
/// regularizer_eps : f64
///     Diagonal jitter $\epsilon$ added to covariance matrices (default: $10^{-6}$).
///
/// Returns
/// -------
/// Array2<f64>
///     Rectangular $N_A \times N_B$ matrix of generalized correlation coefficients $r_{\text{MI}}(i, j) \in [0, 1]$.
pub fn compute_intermolecular_allosteric_network(
    chain_a: &[BackboneTrace],
    chain_b: &[BackboneTrace],
    regularizer_eps: f64,
) -> Result<Array2<f64>, GeometryError> {
    let f_count = chain_a.len();
    if f_count != chain_b.len() {
        return Err(GeometryError::InsufficientPoints {
            required: f_count,
            actual: chain_b.len(),
        });
    }
    if f_count < 3 {
        return Err(GeometryError::InsufficientPoints {
            required: 3,
            actual: f_count,
        });
    }

    let feats_a = extract_centered_features(chain_a, regularizer_eps)?;
    let feats_b = extract_centered_features(chain_b, regularizer_eps)?;

    let d_a = feats_a.feat_dim;
    let d_b = feats_b.feat_dim;
    let joint_dim = d_a + d_b;
    let inv_df = 1.0 / (f_count - 1) as f64;
    let eps = if regularizer_eps > 0.0 { regularizer_eps } else { 1e-6 };

    // Parallelize over rows i in 0..feats_a.n_internal
    let row_results: Vec<Vec<f64>> = (0..feats_a.n_internal)
        .into_par_iter()
        .map(|i| {
            let mut row = Vec::with_capacity(feats_b.n_internal);
            let cov_i = &feats_a.single_covs[i].cov;
            let log_det_i = feats_a.single_covs[i].log_det;

            for j in 0..feats_b.n_internal {
                let cov_j = &feats_b.single_covs[j].cov;
                let log_det_j = feats_b.single_covs[j].log_det;

                // Compute cross-covariance C_ij of shape (d_a, d_b)
                let mut cross = DMatrix::<f64>::zeros(d_a, d_b);
                for f in 0..f_count {
                    let u_i = &feats_a.centered[i][f];
                    let u_j = &feats_b.centered[j][f];
                    for r in 0..d_a {
                        for c in 0..d_b {
                            cross[(r, c)] += u_i[r] * u_j[c];
                        }
                    }
                }
                cross *= inv_df;

                // Assemble (d_a + d_b) x (d_a + d_b) joint covariance matrix
                let mut joint_cov = DMatrix::<f64>::zeros(joint_dim, joint_dim);
                // Top-left: cov_i (d_a x d_a)
                for r in 0..d_a {
                    for c in 0..d_a {
                        joint_cov[(r, c)] = cov_i[(r, c)];
                    }
                }
                // Bottom-right: cov_j (d_b x d_b)
                for r in 0..d_b {
                    for c in 0..d_b {
                        joint_cov[(r + d_a, c + d_a)] = cov_j[(r, c)];
                    }
                }
                // Top-right: cross (d_a x d_b)
                for r in 0..d_a {
                    for c in 0..d_b {
                        joint_cov[(r, c + d_a)] = cross[(r, c)];
                    }
                }
                // Bottom-left: cross.T (d_b x d_a)
                for r in 0..d_b {
                    for c in 0..d_a {
                        joint_cov[(r + d_a, c)] = cross[(c, r)];
                    }
                }

                // Add regularizer to joint covariance
                for d in 0..joint_dim {
                    joint_cov[(d, d)] += eps;
                }

                let log_det_joint = log_det_symm(&joint_cov);
                let mi = -0.5 * (log_det_joint - log_det_i - log_det_j);
                let clamped_mi = if mi > 0.0 { mi } else { 0.0 };

                let r_mi = (1.0 - (-2.0 * clamped_mi).exp()).clamp(0.0, 1.0).sqrt();
                row.push(r_mi);
            }
            row
        })
        .collect();

    // Fill into N_A x N_B rectangular output array
    let mut network = Array2::<f64>::zeros((feats_a.n_residues, feats_b.n_residues));

    // Internal residues: chain_a index i_idx + 1, chain_b index j_idx + 1
    for (i_idx, row) in row_results.into_iter().enumerate() {
        let res_i = i_idx + 1;
        for (j_idx, r_mi) in row.into_iter().enumerate() {
            let res_j = j_idx + 1;
            network[[res_i, res_j]] = r_mi;
        }
    }

    Ok(network)
}

/// Computes the $N_A \times N_B$ Inter-Molecular Allosteric Network for ribbon traces.
///
/// # Errors
/// Returns [`GeometryError`] if trajectory traces have invalid geometry or mismatched lengths.
pub fn compute_intermolecular_allosteric_network_ribbon(
    chain_a: &[RibbonTrace],
    chain_b: &[RibbonTrace],
    regularizer_eps: f64,
) -> Result<Array2<f64>, GeometryError> {
    let bb_a: Vec<BackboneTrace> = chain_a.iter().map(|r| r.to_backbone_trace()).collect();
    let bb_b: Vec<BackboneTrace> = chain_b.iter().map(|r| r.to_backbone_trace()).collect();
    compute_intermolecular_allosteric_network(&bb_a, &bb_b, regularizer_eps)
}

/// Computes the average dynamic allosteric cooperativity index across an inter-protein contact interface.
///
/// Parameters
/// ----------
/// inter_matrix : &Array2<f64>
///     The $N_A \times N_B$ cross-chain generalized correlation matrix $r_{\text{MI}}$.
/// interface_pairs : &[(usize, usize)]
///     List of 0-based residue pairs $(i_A, j_B)$ located within the physical binding interface (e.g. within 5 Å).
///
/// Returns
/// -------
/// f64
///     The mean dynamic cooperativity score across the interface contacts, $\mathcal{I}_{\text{coop}} \in [0, 1]$.
#[must_use]
pub fn compute_ternary_cooperativity_index(
    inter_matrix: &Array2<f64>,
    interface_pairs: &[(usize, usize)],
) -> f64 {
    if interface_pairs.is_empty() {
        return 0.0;
    }
    let n_rows = inter_matrix.nrows();
    let n_cols = inter_matrix.ncols();

    let mut sum = 0.0;
    let mut count = 0;

    for &(i, j) in interface_pairs {
        if i < n_rows && j < n_cols {
            sum += inter_matrix[[i, j]];
            count += 1;
        }
    }

    if count > 0 {
        sum / (count as f64)
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    #[test]
    fn test_allosteric_network_diagonal() {
        // Construct a simple 10-residue synthetic trajectory with 5 frames
        let mut traces = Vec::new();
        for f in 0..5 {
            let mut pts = Vec::new();
            for i in 0..10 {
                let theta = (i as f64) * 0.5 + (f as f64) * 0.05;
                pts.push(Point3::new(
                    (i as f64) * 3.8,
                    3.0 * theta.cos(),
                    3.0 * theta.sin(),
                ));
            }
            traces.push(BackboneTrace::new(pts));
        }

        let net = compute_intrinsic_allosteric_network(&traces, 1e-6).unwrap();
        assert_eq!(net.shape(), &[10, 10]);

        // Diagonal elements should be 1.0
        for i in 0..10 {
            assert!((net[[i, i]] - 1.0).abs() < 1e-6);
        }

        // Symmetry: net[i, j] == net[j, i]
        for i in 0..10 {
            for j in 0..10 {
                assert!((net[[i, j]] - net[[j, i]]).abs() < 1e-6);
                assert!(net[[i, j]] >= 0.0 && net[[i, j]] <= 1.0);
            }
        }
    }

    #[test]
    fn test_intermolecular_allosteric_network_and_cooperativity() {
        // Chain A has 8 residues, Chain B has 12 residues across 10 frames
        let mut chain_a = Vec::new();
        let mut chain_b = Vec::new();

        for f in 0..10 {
            let f_flt = f as f64;
            // Chain A: internal conformational breathing (radius and pitch vary with f)
            let radius_a = 2.0 + f_flt * 0.15;
            let pitch_a = 3.8 + f_flt * 0.08;
            let mut pts_a = Vec::new();
            for i in 0..8 {
                let ang = (i as f64) * 0.6;
                pts_a.push(Point3::new((i as f64) * pitch_a, ang.cos() * radius_a, ang.sin() * radius_a));
            }
            chain_a.push(BackboneTrace::new(pts_a));

            // Chain B: coupled conformational breathing with Chain A
            let radius_b = 2.5 + f_flt * 0.20;
            let pitch_b = 3.8 + f_flt * 0.08;
            let mut pts_b = Vec::new();
            for j in 0..12 {
                let ang = (j as f64) * 0.4;
                pts_b.push(Point3::new((j as f64) * pitch_b, ang.cos() * radius_b, ang.sin() * radius_b));
            }
            chain_b.push(BackboneTrace::new(pts_b));
        }

        let inter_net = compute_intermolecular_allosteric_network(&chain_a, &chain_b, 1e-6).unwrap();
        assert_eq!(inter_net.shape(), &[8, 12]);

        // Check values bounded in [0, 1]
        for i in 0..8 {
            for j in 0..12 {
                assert!(inter_net[[i, j]] >= 0.0 && inter_net[[i, j]] <= 1.0);
            }
        }

        // Test interface cooperativity index
        let interface_pairs = vec![(2, 3), (3, 4), (4, 5)];
        let coop_idx = compute_ternary_cooperativity_index(&inter_net, &interface_pairs);
        assert!(coop_idx > 0.0 && coop_idx <= 1.0);
    }
}

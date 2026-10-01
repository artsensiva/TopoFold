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
#[allow(clippy::needless_range_loop)]
pub fn compute_intrinsic_allosteric_network(
    traces: &[BackboneTrace],
    regularizer_eps: f64,
) -> Result<Array2<f64>, GeometryError> {
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

    // Internal residue indices that possess invariant features:
    // Vertices 1 .. n_residues - 2 (indices 0 .. n_residues - 3 in curvatures)
    let n_internal = n_residues - 2;

    // Feature matrix: [residue_internal_idx][frame_idx] -> Vec<f64> of length feat_dim
    // Raw feature extraction:
    let mut raw_features = vec![vec![vec![0.0; feat_dim]; f_count]; n_internal];

    for (f, inv) in invariants.iter().enumerate() {
        let n_curv = inv.curvatures.len();
        let n_tors = inv.torsions.len();
        let has_sc = !inv.sidechain_dihedrals.is_empty();

        for i in 0..n_internal {
            let k = if i < n_curv { inv.curvatures[i] } else { 0.0 };
            let t = if i < n_tors {
                inv.torsions[i]
            } else if n_tors > 0 {
                inv.torsions[n_tors - 1]
            } else {
                0.0
            };

            raw_features[i][f][0] = k;
            raw_features[i][f][1] = t;

            if has_cbeta {
                let th = if has_sc && i < inv.sidechain_dihedrals.len() {
                    inv.sidechain_dihedrals[i]
                } else {
                    0.0
                };
                raw_features[i][f][2] = th;
            }
        }
    }

    // Step 2: Mean-center features with circular wrapping on torsion and ribbon angles
    let mut centered_features = vec![vec![vec![0.0; feat_dim]; f_count]; n_internal];

    for i in 0..n_internal {
        // Curvature linear mean
        let mean_k: f64 = raw_features[i].iter().map(|v| v[0]).sum::<f64>() / f_count as f64;

        // Torsion circular mean
        let sum_sin_t: f64 = raw_features[i].iter().map(|v| v[1].sin()).sum();
        let sum_cos_t: f64 = raw_features[i].iter().map(|v| v[1].cos()).sum();
        let mean_t = sum_sin_t.atan2(sum_cos_t);

        // Ribbon circular mean
        let mean_th = if has_cbeta {
            let sum_sin_th: f64 = raw_features[i].iter().map(|v| v[2].sin()).sum();
            let sum_cos_th: f64 = raw_features[i].iter().map(|v| v[2].cos()).sum();
            sum_sin_th.atan2(sum_cos_th)
        } else {
            0.0
        };

        for f in 0..f_count {
            // Delta curvature
            centered_features[i][f][0] = raw_features[i][f][0] - mean_k;

            // Delta torsion on S^1
            let raw_t = raw_features[i][f][1];
            let diff_t = (raw_t - mean_t).sin().atan2((raw_t - mean_t).cos());
            centered_features[i][f][1] = diff_t;

            if has_cbeta {
                let raw_th = raw_features[i][f][2];
                let diff_th = (raw_th - mean_th).sin().atan2((raw_th - mean_th).cos());
                centered_features[i][f][2] = diff_th;
            }
        }
    }

    // Step 3: Precompute single-residue covariance matrices and their regularized log determinants
    let eps = if regularizer_eps > 0.0 { regularizer_eps } else { 1e-6 };
    let inv_df = 1.0 / (f_count - 1) as f64;

    struct SingleResidueCov {
        cov: DMatrix<f64>,
        log_det: f64,
    }

    let single_covs: Vec<SingleResidueCov> = (0..n_internal)
        .into_par_iter()
        .map(|i| {
            let mut cov = DMatrix::<f64>::zeros(feat_dim, feat_dim);
            for f in 0..f_count {
                let u = &centered_features[i][f];
                for r in 0..feat_dim {
                    for c in 0..feat_dim {
                        cov[(r, c)] += u[r] * u[c];
                    }
                }
            }
            cov *= inv_df;
            // Add diagonal jitter
            for d in 0..feat_dim {
                cov[(d, d)] += eps;
            }
            let log_det = log_det_symm(&cov);
            SingleResidueCov { cov, log_det }
        })
        .collect();

    // Step 4: Compute pair-wise mutual information and generalized correlation in parallel
    let joint_dim = feat_dim * 2;

    // Parallelize outer loop over rows i in 0..n_internal
    let row_results: Vec<Vec<(usize, f64)>> = (0..n_internal)
        .into_par_iter()
        .map(|i| {
            let mut row = Vec::with_capacity(n_internal - i);
            let cov_i = &single_covs[i].cov;
            let log_det_i = single_covs[i].log_det;

            for j in (i + 1)..n_internal {
                let cov_j = &single_covs[j].cov;
                let log_det_j = single_covs[j].log_det;

                // Compute cross-covariance C_ij_cross of shape (feat_dim, feat_dim)
                let mut cross = DMatrix::<f64>::zeros(feat_dim, feat_dim);
                for f in 0..f_count {
                    let u_i = &centered_features[i][f];
                    let u_j = &centered_features[j][f];
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

    // Fill internal residues (residues 1 .. n_residues - 2, corresponding to PDB res 2 .. N-1)
    for (i_idx, row) in row_results.into_iter().enumerate() {
        let res_i = i_idx + 1; // 1-based index in full chain
        for (j_idx, r_mi) in row {
            let res_j = j_idx + 1;
            network[[res_i, res_j]] = r_mi;
            network[[res_j, res_i]] = r_mi;
        }
    }

    Ok(network)
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
                    10.0 * theta.cos(),
                    10.0 * theta.sin(),
                    (i as f64) * 3.8,
                ));
            }
            traces.push(BackboneTrace::new(pts));
        }

        let net = compute_intrinsic_allosteric_network(&traces, 1e-6).expect("Network computation failed");
        assert_eq!(net.shape(), &[10, 10]);

        for i in 0..10 {
            assert!((net[[i, i]] - 1.0).abs() < 1e-6, "Diagonal must be 1.0");
        }

        // Test symmetry
        for i in 0..10 {
            for j in 0..10 {
                assert!(
                    (net[[i, j]] - net[[j, i]]).abs() < 1e-12,
                    "Matrix must be strictly symmetric"
                );
                assert!(net[[i, j]] >= 0.0 && net[[i, j]] <= 1.0, "Values must be in [0, 1]");
            }
        }
    }
}

//! Intrinsic Allosteric Networks via Mutual Information on Differential Curve Invariants.
//!
//! Maps conformational trajectories into an $N \times N$ communication adjacency matrix
//! by computing the Generalized Correlation Coefficient $r_{\text{MI}}$ from the
//! Mutual Information $I(\mathbf{V}_i; \mathbf{V}_j)$ of discrete curve invariants:
//! $\mathbf{V}_i(t) = [\kappa_i(t), \tau_i(t), \theta_{\beta, i}(t)] \in \mathbb{R}^3$.
//!
//! The estimator is a **regularized Gaussian (covariance) MI** on the centered per-residue features: the formula
//! below is the MI of the Gaussian model fitted with the regularized sample covariances. It is a finite-sample
//! plug-in estimate, not the population MI, even when the features are jointly Gaussian. No rank (copula)
//! transform, no circular model and no significance test is applied. Torsion and $\theta_\beta$ are centered on their circular means and wrapped to $(-\pi, \pi]$,
//! then treated as linear variables. With unregularized sample covariances $\mathbf{S}_i, \mathbf{S}_j$ and joint
//! $\mathbf{S}_{ij}$ (denominator $F - 1$), and a single ridge $\epsilon$ added once to every diagonal:
//! $$I_\epsilon(\mathbf{V}_i; \mathbf{V}_j) = \frac{1}{2} \ln \left( \frac{\det(\mathbf{S}_i + \epsilon \mathbf{I}) \det(\mathbf{S}_j + \epsilon \mathbf{I})}{\det(\mathbf{S}_{ij} + \epsilon \mathbf{I})} \right) \ge 0.$$
//! $I_\epsilon$ is not the unregularized Gaussian MI $I_0$ and can differ from it when $\epsilon$ is not negligible
//! against the feature variances. The ridge need not reduce the value in every case; features with zero sample
//! cross-covariance give $I_\epsilon = I_0 = 0$. $I_\epsilon$ stays finite when $\mathbf{S}_{ij}$ is singular.
//! The $\ge 0$ bound holds in exact arithmetic for positive-definite matrices; computed values can differ (see
//! `log_det_symm` for the determinant fallback).
//! *Change 2026-10-22 (KNOWN_ISSUES M7):* earlier versions added $\epsilon$ to the joint diagonal a second time,
//! which biased $I$ towards zero and, after clipping, returned 0 for correlated low-variance features. Outputs of
//! both network functions can therefore differ from earlier versions.
//!
//! The generalized correlation coefficient:
//! $$r_{\text{MI}}(i, j) = \sqrt{1 - \exp(-2 I_\epsilon(\mathbf{V}_i; \mathbf{V}_j))} \in [0, 1]$$
//! is $SE(3)$-invariant and superposition-free because the features are. It measures statistical association within
//! one ensemble only. Frames are treated as independent samples (temporal dependence is not handled), and no
//! causal or biological (allosteric) interpretation follows from it.

#![forbid(unsafe_code)]

use nalgebra::DMatrix;
use ndarray::Array2;
use rayon::prelude::*;

use crate::error::GeometryError;
use crate::geometry::extract_curve_invariants;
use crate::ribbon::RibbonTrace;
use crate::types::BackboneTrace;

/// Computes the log determinant of a symmetric positive-definite matrix using Cholesky decomposition.
///
/// Known limitation (not changed in the 2026-10-22 M7 repair): if Cholesky fails, e.g. when $\epsilon$ is negligible
/// against the feature scale and the matrix is numerically singular, the fallback uses the LU determinant only if it
/// exceeds 1e-30, otherwise adds an extra 1e-4 jitter to this matrix alone, and otherwise returns 0. In those cases
/// the marginal and joint terms are no longer regularized consistently, so the resulting MI is not $I_\epsilon$ and
/// may be negative (then clipped to 0) or otherwise wrong.
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
    let eps = if regularizer_eps > 0.0 {
        regularizer_eps
    } else {
        1e-6
    };
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

/// Regularized Gaussian MI of residues i and j:
/// $I_\epsilon = \tfrac{1}{2}[\ln\det(S_i + \epsilon I) + \ln\det(S_j + \epsilon I) - \ln\det(S_{ij} + \epsilon I)]$.
///
/// `cov_i` and `cov_j` are the marginal covariances **already regularized** ($S + \epsilon I$) with their
/// log-determinants. The joint matrix is assembled from these blocks plus the unregularized cross-covariance, so
/// its diagonal carries $\epsilon$ exactly once, consistent with the marginals. In exact arithmetic
/// $I_\epsilon \ge 0$ by Fischer's inequality. Negative computed values, from rounding or from the determinant
/// fallback in `log_det_symm`, are clipped to 0.
fn regularized_gaussian_mi(
    cov_i: &DMatrix<f64>,
    log_det_i: f64,
    cov_j: &DMatrix<f64>,
    log_det_j: f64,
    cross: &DMatrix<f64>,
) -> f64 {
    let d_a = cov_i.nrows();
    let d_b = cov_j.nrows();
    let joint_dim = d_a + d_b;

    let mut joint_cov = DMatrix::<f64>::zeros(joint_dim, joint_dim);
    joint_cov.view_mut((0, 0), (d_a, d_a)).copy_from(cov_i);
    joint_cov.view_mut((d_a, d_a), (d_b, d_b)).copy_from(cov_j);
    joint_cov.view_mut((0, d_a), (d_a, d_b)).copy_from(cross);
    joint_cov
        .view_mut((d_a, 0), (d_b, d_a))
        .copy_from(&cross.transpose());

    let log_det_joint = log_det_symm(&joint_cov);
    let mi = -0.5 * (log_det_joint - log_det_i - log_det_j);
    if mi > 0.0 {
        mi
    } else {
        0.0
    }
}

/// Computes the $N \times N$ Intrinsic Allosteric Network Matrix across a trajectory of conformers.
///
/// Parameters
/// ----------
/// traces : &[BackboneTrace]
///     Ensemble of $F$ conformational frames, each with $N$ residues.
/// regularizer_eps : f64
///     Ridge $\epsilon$ added once to the diagonal of every marginal and joint sample covariance; non-positive
///     values select $10^{-6}$. The result is the regularized quantity $I_\epsilon$ (module docs), not $I_0$.
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
    let inv_df = 1.0 / (f_count - 1) as f64;

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

                let clamped_mi =
                    regularized_gaussian_mi(cov_i, log_det_i, cov_j, log_det_j, &cross);

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
///     Ridge $\epsilon$ added once to the diagonal of every marginal and joint sample covariance; non-positive
///     values select $10^{-6}$. The result is the regularized quantity $I_\epsilon$ (module docs), not $I_0$.
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
    let inv_df = 1.0 / (f_count - 1) as f64;

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

                let clamped_mi =
                    regularized_gaussian_mi(cov_i, log_det_i, cov_j, log_det_j, &cross);

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
                pts_a.push(Point3::new(
                    (i as f64) * pitch_a,
                    ang.cos() * radius_a,
                    ang.sin() * radius_a,
                ));
            }
            chain_a.push(BackboneTrace::new(pts_a));

            // Chain B: coupled conformational breathing with Chain A
            let radius_b = 2.5 + f_flt * 0.20;
            let pitch_b = 3.8 + f_flt * 0.08;
            let mut pts_b = Vec::new();
            for j in 0..12 {
                let ang = (j as f64) * 0.4;
                pts_b.push(Point3::new(
                    (j as f64) * pitch_b,
                    ang.cos() * radius_b,
                    ang.sin() * radius_b,
                ));
            }
            chain_b.push(BackboneTrace::new(pts_b));
        }

        let inter_net =
            compute_intermolecular_allosteric_network(&chain_a, &chain_b, 1e-6).unwrap();
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

    // --- M7 regression tests: regularization consistency (2026-10-22) ---
    //
    // Estimand. For unregularized sample covariances S_i, S_j and joint S_ij (blocks S_i, S_j and the
    // cross-covariance), the regularized Gaussian MI with a single, consistent ridge ε is
    //   I_ε = ½ [ln det(S_i + εI) + ln det(S_j + εI) − ln det(S_ij + εI)].
    // It is ≥ 0 (Fischer's inequality) and differs from the unregularized Gaussian MI
    //   I_0 = ½ [ln det S_i + ln det S_j − ln det S_ij]
    // whenever ε is not negligible against the covariances. The references below are computed
    // independently of `log_det_symm` (closed forms, or symmetric eigenvalues).

    fn log_det_eig(m: &DMatrix<f64>) -> f64 {
        m.clone()
            .symmetric_eigen()
            .eigenvalues
            .iter()
            .map(|v| v.ln())
            .sum()
    }

    fn ridge(m: &DMatrix<f64>, eps: f64) -> DMatrix<f64> {
        m + DMatrix::<f64>::identity(m.nrows(), m.ncols()) * eps
    }

    /// Reference I_ε from unregularized blocks (independent of the code under test).
    fn reference_mi(si: &DMatrix<f64>, sj: &DMatrix<f64>, c: &DMatrix<f64>, eps: f64) -> f64 {
        let (da, db) = (si.nrows(), sj.nrows());
        let mut joint = DMatrix::<f64>::zeros(da + db, da + db);
        joint.view_mut((0, 0), (da, da)).copy_from(si);
        joint.view_mut((da, da), (db, db)).copy_from(sj);
        joint.view_mut((0, da), (da, db)).copy_from(c);
        joint.view_mut((da, 0), (db, da)).copy_from(&c.transpose());
        0.5 * (log_det_eig(&ridge(si, eps)) + log_det_eig(&ridge(sj, eps))
            - log_det_eig(&ridge(&joint, eps)))
    }

    /// Calls the helper exactly as the public paths do: marginals arrive already regularized.
    fn helper_mi(si: &DMatrix<f64>, sj: &DMatrix<f64>, c: &DMatrix<f64>, eps: f64) -> f64 {
        let (ri, rj) = (ridge(si, eps), ridge(sj, eps));
        let (li, lj) = (log_det_symm(&ri), log_det_symm(&rj));
        regularized_gaussian_mi(&ri, li, &rj, lj, c)
    }

    fn iso(d: usize, v: f64) -> DMatrix<f64> {
        DMatrix::<f64>::identity(d, d) * v
    }

    #[test]
    fn test_mi_helper_matches_closed_form_correlated_low_variance() {
        // S_i = S_j = aI_d, cross = cI_d: det(S_ij + εI) = ((a+ε)² − c²)^d, so
        // I_ε = (d/2) ln[(a+ε)² / ((a+ε)² − c²)]; I_0 = −(d/2) ln(1 − (c/a)²).
        let eps = 1e-6;
        for &d in &[2usize, 3] {
            for &(a, c) in &[
                (1e-6f64, 0.9e-6f64),
                (2e-6, 1.5e-6),
                (1e-3, 0.5e-3),
                (1e-6, 0.0),
            ] {
                let df = d as f64;
                let expected = 0.5 * df * ((a + eps).powi(2) / ((a + eps).powi(2) - c * c)).ln();
                let got = helper_mi(&iso(d, a), &iso(d, a), &iso(d, c), eps);
                assert!(
                    (got - expected).abs() <= 1e-9 * expected.abs().max(1e-3),
                    "d={d} a={a} c={c}: got {got}, closed form {expected}"
                );
                if c > 0.0 && a <= 1e-5 {
                    // ε is material here: the regularized quantity is not the unregularized MI.
                    let unregularized = -0.5 * df * (1.0 - (c / a).powi(2)).ln();
                    assert!(
                        got < 0.5 * unregularized,
                        "I_ε should be well below I_0 when a ~ ε"
                    );
                }
            }
        }
        // The published-analysis example: Σ_i = Σ_j = εI₂, C = 0.9εI₂ gives I_ε = ln(4/3.19) > 0.
        let got = helper_mi(&iso(2, eps), &iso(2, eps), &iso(2, 0.9 * eps), eps);
        assert!((got - (4.0f64 / 3.19).ln()).abs() < 1e-9, "got {got}");
    }

    #[test]
    fn test_mi_helper_matches_eigen_reference_general_blocks() {
        // Non-isotropic, unequal dimensions (2 × 3), deterministic positive semi-definite blocks.
        let eps = 1e-4;
        let x = DMatrix::<f64>::from_fn(5, 7, |r, c| ((r * 7 + c) as f64 * 0.37).sin() * 1e-2);
        let s = &x * x.transpose() / 6.0; // 5×5 sample-covariance-like PSD matrix
        let si = s.view((0, 0), (2, 2)).into_owned();
        let sj = s.view((2, 2), (3, 3)).into_owned();
        let c = s.view((0, 2), (2, 3)).into_owned();
        let expected = reference_mi(&si, &sj, &c, eps);
        assert!(expected > 1e-3, "fixture should be informative: {expected}");
        let got = helper_mi(&si, &sj, &c, eps);
        assert!(
            (got - expected).abs() < 1e-10,
            "got {got}, reference {expected}"
        );
        // Symmetry under exchange of i and j.
        let swapped = helper_mi(&sj, &si, &c.transpose(), eps);
        assert!((got - swapped).abs() < 1e-12);
    }

    #[test]
    fn test_mi_helper_independence_and_degenerate_inputs() {
        let eps = 1e-6;
        // Independence (zero cross-covariance): I_ε = 0 for any marginal scale.
        for &a in &[0.0, 1e-8, 1e-6, 1.0] {
            let got = helper_mi(&iso(3, a), &iso(2, 2.0 * a), &DMatrix::zeros(3, 2), eps);
            assert!(got.abs() < 1e-12, "a={a}: {got}");
        }
        // Constant features (all covariances zero): regularized blocks are εI, MI is 0.
        assert!(helper_mi(&iso(2, 0.0), &iso(2, 0.0), &iso(2, 0.0), eps).abs() < 1e-12);
        // Perfect correlation (S_ij singular): I_0 is infinite, I_ε is finite and matches the closed form.
        let a: f64 = 1e-6;
        let expected = ((a + eps).powi(2) / ((a + eps).powi(2) - a * a)).ln();
        let got = helper_mi(&iso(2, a), &iso(2, a), &iso(2, a), eps);
        assert!(
            got.is_finite() && (got - expected).abs() < 1e-9,
            "got {got}, closed form {expected}"
        );
    }

    /// One coupled degree of freedom `amp·sin(0.7f + phase)` deforms every residue of a helix-like chain.
    fn coupled_chain(n: usize, frame: usize, amp: f64, phase: f64, with_cb: bool) -> BackboneTrace {
        let s = amp * ((frame as f64) * 0.7 + phase).sin();
        let pt = |i: usize, r: f64| {
            let a = i as f64 * 1.745 + s * (i as f64 * 0.9).cos();
            Point3::new(
                r * a.cos(),
                r * a.sin(),
                1.5 * i as f64 + s * (i as f64 * 0.5).sin(),
            )
        };
        let ca: Vec<_> = (0..n).map(|i| pt(i, 2.3)).collect();
        if with_cb {
            let cb: Vec<_> = (0..n)
                .map(|i| pt(i, 3.6) + nalgebra::Vector3::new(0.0, 0.0, 0.4))
                .collect();
            BackboneTrace::with_cbeta(ca, cb)
        } else {
            BackboneTrace::new(ca)
        }
    }

    fn sample_blocks(
        fa: &CenteredChainFeatures,
        i: usize,
        fb: &CenteredChainFeatures,
        j: usize,
    ) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
        let n = fa.centered[i].len();
        let xa = DMatrix::<f64>::from_fn(n, fa.feat_dim, |f, k| fa.centered[i][f][k]);
        let xb = DMatrix::<f64>::from_fn(n, fb.feat_dim, |f, k| fb.centered[j][f][k]);
        let df = (n - 1) as f64;
        (
            xa.transpose() * &xa / df,
            xb.transpose() * &xb / df,
            xa.transpose() * &xb / df,
        )
    }

    fn mi_from_r(r: f64) -> f64 {
        -0.5 * (1.0 - r * r).ln()
    }

    /// Mean diagonal variance of the centered features; used to set ε on the scale of the data.
    fn feature_scale(f: &CenteredChainFeatures) -> f64 {
        let n = f.centered[0].len() as f64;
        let mut sum = 0.0;
        for res in &f.centered {
            for frame in res {
                sum += frame.iter().map(|v| v * v).sum::<f64>();
            }
        }
        sum / ((n - 1.0) * (f.centered.len() * f.feat_dim) as f64)
    }

    #[test]
    fn test_intrinsic_network_matches_reference_correlated_low_variance() {
        for &with_cb in &[false, true] {
            let traces: Vec<_> = (0..16)
                .map(|f| coupled_chain(10, f, 0.01, 0.0, with_cb))
                .collect();
            let feats = extract_centered_features(&traces, 1.0).unwrap();
            let eps = feature_scale(&feats); // ε comparable to the feature variances
            let net = compute_intrinsic_allosteric_network(&traces, eps).unwrap();
            let mut informative = 0;
            for i in 0..feats.n_internal {
                for j in (i + 1)..feats.n_internal {
                    let (si, sj, c) = sample_blocks(&feats, i, &feats, j);
                    let expected = reference_mi(&si, &sj, &c, eps);
                    let got = mi_from_r(net[[i + 1, j + 1]]);
                    assert!(
                        (got - expected).abs() < 1e-8 * expected.max(1.0),
                        "cb={with_cb} ({i},{j}): got {got}, reference {expected}"
                    );
                    assert_eq!(net[[i + 1, j + 1]].to_bits(), net[[j + 1, i + 1]].to_bits());
                    if expected > 0.05 {
                        informative += 1;
                    }
                }
            }
            assert!(informative > 0, "fixture must contain correlated pairs");
        }
    }

    #[test]
    fn test_intermolecular_network_matches_reference_and_intrinsic() {
        for &with_cb in &[false, true] {
            let a: Vec<_> = (0..16)
                .map(|f| coupled_chain(9, f, 0.01, 0.0, with_cb))
                .collect();
            let b: Vec<_> = (0..16)
                .map(|f| coupled_chain(7, f, 0.01, 0.4, with_cb))
                .collect();
            let fa = extract_centered_features(&a, 1.0).unwrap();
            let fb = extract_centered_features(&b, 1.0).unwrap();
            let eps = 0.5 * (feature_scale(&fa) + feature_scale(&fb));
            let ab = compute_intermolecular_allosteric_network(&a, &b, eps).unwrap();
            let ba = compute_intermolecular_allosteric_network(&b, &a, eps).unwrap();
            let mut informative = 0;
            for i in 0..fa.n_internal {
                for j in 0..fb.n_internal {
                    let (si, sj, c) = sample_blocks(&fa, i, &fb, j);
                    let expected = reference_mi(&si, &sj, &c, eps);
                    let got = mi_from_r(ab[[i + 1, j + 1]]);
                    assert!(
                        (got - expected).abs() < 1e-8 * expected.max(1.0),
                        "cb={with_cb} ({i},{j}): got {got}, reference {expected}"
                    );
                    assert!((ab[[i + 1, j + 1]] - ba[[j + 1, i + 1]]).abs() < 1e-12);
                    if expected > 0.05 {
                        informative += 1;
                    }
                }
            }
            assert!(informative > 0, "fixture must contain correlated pairs");

            // Same chain through both paths: off-diagonal entries agree.
            let intra = compute_intrinsic_allosteric_network(&a, eps).unwrap();
            let inter = compute_intermolecular_allosteric_network(&a, &a, eps).unwrap();
            for i in 1..=fa.n_internal {
                for j in 1..=fa.n_internal {
                    if i != j {
                        assert!((intra[[i, j]] - inter[[i, j]]).abs() < 1e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn test_network_degenerate_inputs() {
        // Identical frames: every feature is constant, so every off-diagonal coefficient is 0.
        let still: Vec<_> = (0..6)
            .map(|_| coupled_chain(8, 0, 0.0, 0.0, false))
            .collect();
        let net = compute_intrinsic_allosteric_network(&still, 1e-6).unwrap();
        let inter = compute_intermolecular_allosteric_network(&still, &still, 1e-6).unwrap();
        for i in 0..8 {
            for j in 0..8 {
                if i != j {
                    assert_eq!(net[[i, j]], 0.0);
                }
                assert_eq!(inter[[i, j]], 0.0);
            }
        }
        // Non-positive ε selects the documented default 1e-6.
        let traces: Vec<_> = (0..10)
            .map(|f| coupled_chain(8, f, 0.05, 0.0, false))
            .collect();
        let reference = compute_intrinsic_allosteric_network(&traces, 1e-6).unwrap();
        for eps in [0.0, -1.0] {
            assert_eq!(
                compute_intrinsic_allosteric_network(&traces, eps).unwrap(),
                reference
            );
        }
        // Too few frames is an error, not a silent result.
        assert!(compute_intrinsic_allosteric_network(&traces[..2], 1e-6).is_err());
    }
}

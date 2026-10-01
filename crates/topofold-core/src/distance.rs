//! Metric distance functions between protein backbone curve signatures.
//!
//! Provides two tiers of metric comparison:
//! 1. Coarse Pass: $L_2$ Euclidean metric on localized Writhe spectra.
//! 2. Fine Pass: Discrete Fréchet distance and Dynamic Time Warping (DTW) on $(\kappa, \tau)$
//!    discrete differential invariants, incorporating the $S^1$ branch-cut for torsion.

#![forbid(unsafe_code)]

use std::f64::consts::PI;

use crate::invariants::CurveInvariants;

/// Computes the minimal angular distance on the circle $S^1$ between two dihedral angles in $(-\pi, \pi]$.
#[inline]
pub fn angular_distance_s1(a: f64, b: f64) -> f64 {
    let mut diff = (a - b).abs();
    if diff > PI {
        diff = (2.0 * PI - diff).abs();
    }
    diff
}

/// Pointwise distance between two local differential invariants $(\kappa_1, \tau_1)$ and $(\kappa_2, \tau_2)$.
///
/// Curvature $\kappa \in [0, \pi]$ is Euclidean, while torsion $\tau \in (-\pi, \pi]$ wraps periodically.
#[inline]
pub fn invariant_point_distance(
    kappa1: f64,
    tau1: f64,
    kappa2: f64,
    tau2: f64,
    weight_kappa: f64,
    weight_tau: f64,
) -> f64 {
    let d_kappa = (kappa1 - kappa2).abs();
    let d_tau = angular_distance_s1(tau1, tau2);
    (weight_kappa * d_kappa * d_kappa + weight_tau * d_tau * d_tau).sqrt()
}

/// Computes the $L_2$ Euclidean distance between two Writhe spectrum vectors.
///
/// Used as a rapid coarse filter prior to non-linear curve alignment.
#[must_use]
pub fn writhe_spectrum_distance(a: &[f64], b: &[f64]) -> f64 {
    let min_len = a.len().min(b.len());
    if min_len == 0 {
        return 0.0;
    }

    let mut sum_sq = 0.0;
    for i in 0..min_len {
        let diff = a[i] - b[i];
        sum_sq += diff * diff;
    }

    (sum_sq / min_len as f64).sqrt()
}

/// Computes the Discrete Fréchet Distance between two $(\kappa, \tau)$ curve invariant sequences.
///
/// Evaluates the minimax coupling metric:
/// $$\delta_F(P, Q) = \min_{\alpha, \beta} \max_t d(P(\alpha(t)), Q(\beta(t)))$$
///
/// Implemented using dynamic programming with an $O(\min(m, n))$ space-optimized rolling buffer.
#[must_use]
pub fn discrete_frechet_invariants(
    inv_a: &CurveInvariants,
    inv_b: &CurveInvariants,
    weight_kappa: f64,
    weight_tau: f64,
) -> f64 {
    let (k_a, t_a) = (&inv_a.curvatures, &inv_a.torsions);
    let (k_b, t_b) = (&inv_b.curvatures, &inv_b.torsions);

    let len_a = k_a.len().min(t_a.len());
    let len_b = k_b.len().min(t_b.len());

    if len_a == 0 || len_b == 0 {
        return 0.0;
    }

    // Dynamic programming matrix with 2 rows for space optimization
    let mut prev_row = vec![f64::INFINITY; len_b];
    let mut curr_row = vec![f64::INFINITY; len_b];

    // Base case: (0, 0)
    prev_row[0] = invariant_point_distance(
        k_a[0],
        t_a[0],
        k_b[0],
        t_b[0],
        weight_kappa,
        weight_tau,
    );

    // First row: coupling with a[0]
    for j in 1..len_b {
        let d = invariant_point_distance(
            k_a[0],
            t_a[0],
            k_b[j],
            t_b[j],
            weight_kappa,
            weight_tau,
        );
        prev_row[j] = prev_row[j - 1].max(d);
    }

    // Iterate through subsequent rows
    for i in 1..len_a {
        // First column: coupling with b[0]
        let d_first = invariant_point_distance(
            k_a[i],
            t_a[i],
            k_b[0],
            t_b[0],
            weight_kappa,
            weight_tau,
        );
        curr_row[0] = prev_row[0].max(d_first);

        for j in 1..len_b {
            let d = invariant_point_distance(
                k_a[i],
                t_a[i],
                k_b[j],
                t_b[j],
                weight_kappa,
                weight_tau,
            );
            let min_pred = prev_row[j].min(curr_row[j - 1]).min(prev_row[j - 1]);
            curr_row[j] = min_pred.max(d);
        }

        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    prev_row[len_b - 1]
}

/// Computes the Dynamic Time Warping (DTW) distance between two $(\kappa, \tau)$ sequences.
///
/// Accumulates the minimum alignment cost across both curves.
#[must_use]
pub fn dtw_invariants(
    inv_a: &CurveInvariants,
    inv_b: &CurveInvariants,
    weight_kappa: f64,
    weight_tau: f64,
) -> f64 {
    let (k_a, t_a) = (&inv_a.curvatures, &inv_a.torsions);
    let (k_b, t_b) = (&inv_b.curvatures, &inv_b.torsions);

    let len_a = k_a.len().min(t_a.len());
    let len_b = k_b.len().min(t_b.len());

    if len_a == 0 || len_b == 0 {
        return 0.0;
    }

    let mut prev_row = vec![f64::INFINITY; len_b];
    let mut curr_row = vec![f64::INFINITY; len_b];

    prev_row[0] = invariant_point_distance(
        k_a[0],
        t_a[0],
        k_b[0],
        t_b[0],
        weight_kappa,
        weight_tau,
    );

    for j in 1..len_b {
        let d = invariant_point_distance(
            k_a[0],
            t_a[0],
            k_b[j],
            t_b[j],
            weight_kappa,
            weight_tau,
        );
        prev_row[j] = prev_row[j - 1] + d;
    }

    for i in 1..len_a {
        let d_first = invariant_point_distance(
            k_a[i],
            t_a[i],
            k_b[0],
            t_b[0],
            weight_kappa,
            weight_tau,
        );
        curr_row[0] = prev_row[0] + d_first;

        for j in 1..len_b {
            let d = invariant_point_distance(
                k_a[i],
                t_a[i],
                k_b[j],
                t_b[j],
                weight_kappa,
                weight_tau,
            );
            let min_pred = prev_row[j].min(curr_row[j - 1]).min(prev_row[j - 1]);
            curr_row[j] = min_pred + d;
        }

        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    prev_row[len_b - 1] / (len_a + len_b) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_angular_distance() {
        assert!((angular_distance_s1(0.0, 0.5) - 0.5).abs() < 1e-12);
        assert!((angular_distance_s1(3.0, -3.0) - (2.0 * PI - 6.0)).abs() < 1e-12);
        assert!((angular_distance_s1(PI, -PI) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_frechet_identity() {
        let inv = CurveInvariants {
            segment_lengths: vec![3.8, 3.8, 3.8],
            curvatures: vec![1.5, 1.5, 1.5],
            torsions: vec![0.8, 0.8],
            sidechain_dihedrals: vec![],
        };

        let dist = discrete_frechet_invariants(&inv, &inv, 1.0, 1.0);
        assert!(dist < 1e-12, "Fréchet distance to itself must be 0, got {dist}");
    }
}

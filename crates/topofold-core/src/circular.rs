//! Circular statistics for angular (periodic) data.
//!
//! Angular quantities like torsion angles τ ∈ (−π, π] and side-chain orientation
//! angles θβ are periodic and require circular statistics. Linear averaging fails
//! near the branch cut at ±π, producing false bimodality signals for actually
//! unimodal distributions (e.g., β-strand torsions near −170°).
//!
//! # Mathematical Background
//!
//! For a set of angles {θ₁, ..., θₙ}, the circular mean θ̄ is:
//!
//! ```text
//! θ̄ = atan2(mean(sin(θᵢ)), mean(cos(θᵢ)))
//! ```
//!
//! This is equivalent to computing the argument of the mean unit vector:
//!
//! ```text
//! θ̄ = arg(mean(exp(iθᵢ)))
//! ```
//!
//! The circular variance is:
//!
//! ```text
//! V = 1 − R
//! ```
//!
//! where R is the mean resultant length:
//!
//! ```text
//! R = |mean(exp(iθᵢ))|
//! ```
//!
//! # References
//!
//! - Mardia, K. V. & Jupp, P. E. (2000). *Directional Statistics*. Wiley.
//! - Fisher, N. I. (1993). *Statistical Analysis of Circular Data*. Cambridge University Press.

use std::f64::consts::PI;

/// Computes the circular mean of angles in radians.
///
/// For angles θ ∈ (−π, π], the circular mean is defined as:
///
/// ```text
/// θ̄ = atan2(Σ sin(θᵢ)/n, Σ cos(θᵢ)/n)
/// ```
///
/// This correctly handles angles near the ±π branch cut.
///
/// # Arguments
/// * `angles` - Slice of angles in radians, typically in (−π, π]
///
/// # Returns
/// Circular mean angle in radians, in (−π, π]
///
/// # Examples
/// ```
/// use topofold_core::circular::circular_mean;
/// use std::f64::consts::PI;
///
/// // Angles clustered near +π
/// let angles = vec![PI - 0.1, PI - 0.05, PI, -PI + 0.1];
/// let mean = circular_mean(&angles);
/// assert!((mean - PI).abs() < 0.2);
///
/// // Angles clustered near −π
/// let angles = vec![-PI + 0.1, -PI + 0.05, PI - 0.05];
/// let mean = circular_mean(&angles);
/// assert!((mean.abs() - PI).abs() < 0.2);
/// ```
#[must_use]
pub fn circular_mean(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }

    let n = angles.len() as f64;
    let sum_sin: f64 = angles.iter().map(|&a| a.sin()).sum();
    let sum_cos: f64 = angles.iter().map(|&a| a.cos()).sum();

    let mean_sin = sum_sin / n;
    let mean_cos = sum_cos / n;

    mean_sin.atan2(mean_cos)
}

/// Computes the mean resultant length R for circular data.
///
/// R measures the concentration of angular data around the mean direction:
///
/// ```text
/// R = √((Σ sin(θᵢ)/n)² + (Σ cos(θᵢ)/n)²)
/// ```
///
/// # Properties
/// - R = 1: All angles identical (perfect concentration)
/// - R ≈ 0: Uniform distribution around the circle
///
/// # Arguments
/// * `angles` - Slice of angles in radians
///
/// # Returns
/// Mean resultant length R ∈ [0, 1]
#[must_use]
pub fn mean_resultant_length(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }

    let n = angles.len() as f64;
    let sum_sin: f64 = angles.iter().map(|&a| a.sin()).sum();
    let sum_cos: f64 = angles.iter().map(|&a| a.cos()).sum();

    let mean_sin = sum_sin / n;
    let mean_cos = sum_cos / n;

    (mean_sin * mean_sin + mean_cos * mean_cos).sqrt()
}

/// Computes the circular variance V = 1 − R.
///
/// Circular variance ranges from 0 (perfect concentration) to 1 (uniform).
///
/// # Arguments
/// * `angles` - Slice of angles in radians
///
/// # Returns
/// Circular variance V ∈ [0, 1]
#[must_use]
pub fn circular_variance(angles: &[f64]) -> f64 {
    1.0 - mean_resultant_length(angles)
}

/// Computes centered angular deviations from the circular mean.
///
/// For each angle θᵢ, computes the signed angular distance to the circular mean θ̄,
/// wrapped to (−π, π]. Useful for computing moment-based statistics on centered data.
///
/// # Arguments
/// * `angles` - Slice of angles in radians
///
/// # Returns
/// Vector of deviations dᵢ = wrap(θᵢ − θ̄), each in (−π, π]
///
/// # Examples
/// ```
/// use topofold_core::circular::centered_deviations;
/// use std::f64::consts::PI;
///
/// // Symmetric angles around 0 - mean will be exactly 0
/// let angles = vec![0.1, -0.1, 0.15, -0.15];
/// let devs = centered_deviations(&angles);
/// // Mean is 0, so deviations = original angles
/// assert!((devs[0] - 0.1).abs() < 1e-10);
/// assert!((devs[1] + 0.1).abs() < 1e-10);
/// ```
#[must_use]
pub fn centered_deviations(angles: &[f64]) -> Vec<f64> {
    if angles.is_empty() {
        return Vec::new();
    }

    let mean = circular_mean(angles);
    angles
        .iter()
        .map(|&a| {
            let diff = a - mean;
            // Wrap to (−π, π]
            let wrapped = ((diff + PI) % (2.0 * PI)) - PI;
            // Handle negative modulo edge case
            if wrapped <= -PI {
                wrapped + 2.0 * PI
            } else {
                wrapped
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_circular_mean_basic() {
        // Mean of angles around 0°
        let angles = vec![0.0, 0.1, -0.1, 0.05];
        let mean = circular_mean(&angles);
        assert_relative_eq!(mean, 0.0125, epsilon = 0.01);
    }

    #[test]
    fn test_circular_mean_near_positive_pi() {
        // Angles clustered near +π
        let angles = vec![PI - 0.1, PI - 0.05, PI, -PI + 0.1];
        let mean = circular_mean(&angles);
        // Should be near +π (or equivalently -π)
        assert!((mean.abs() - PI).abs() < 0.1);
    }

    #[test]
    fn test_circular_mean_near_negative_pi() {
        // Angles clustered near −π
        let angles = vec![-PI + 0.1, -PI + 0.05, -PI + 0.15];
        let mean = circular_mean(&angles);
        // Should be near −π (≈ -3.04 or equivalently +3.24)
        // The mean should be around -3.04, so abs(mean) ≈ 3.04, close to π
        assert!(
            mean < -3.0 && mean > -PI,
            "Mean should be near -π, got {mean}"
        );
        assert_relative_eq!(mean, -PI + 0.1, epsilon = 0.05);
    }

    #[test]
    fn test_circular_mean_branch_cut() {
        // Half near +π, half near −π (which are the same direction)
        let angles = vec![PI - 0.01, -PI + 0.01, PI - 0.02, -PI + 0.02];
        let mean = circular_mean(&angles);
        // Linear mean would be ~0, circular mean should be near ±π
        assert!((mean.abs() - PI).abs() < 0.05);
    }

    #[test]
    fn test_mean_resultant_length_concentrated() {
        // All angles identical → R = 1
        let angles = vec![0.5; 100];
        let r = mean_resultant_length(&angles);
        assert_relative_eq!(r, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mean_resultant_length_uniform() {
        // Approximately uniform distribution → R ≈ 0
        let n = 360;
        let angles: Vec<f64> = (0..n)
            .map(|i| 2.0 * PI * (i as f64) / (n as f64) - PI)
            .collect();
        let r = mean_resultant_length(&angles);
        assert!(r < 0.01, "Uniform distribution should have R ≈ 0, got {r}");
    }

    #[test]
    fn test_circular_variance() {
        // Concentrated: V ≈ 0
        let angles = vec![0.0, 0.01, -0.01, 0.005];
        let v = circular_variance(&angles);
        assert!(v < 0.01);

        // Dispersed: V ≈ 1
        let angles_uniform: Vec<f64> = (0..360)
            .map(|i| 2.0 * PI * (i as f64) / 360.0 - PI)
            .collect();
        let v_uniform = circular_variance(&angles_uniform);
        assert!(v_uniform > 0.99);
    }

    #[test]
    fn test_centered_deviations() {
        let angles = vec![0.1, 0.2, -0.1, 0.0];
        let devs = centered_deviations(&angles);

        // Deviations should have circular mean ≈ 0
        let dev_mean = circular_mean(&devs);
        assert_relative_eq!(dev_mean, 0.0, epsilon = 1e-10);

        // Standard deviation of deviations should match input spread
        assert_eq!(devs.len(), angles.len());
    }

    #[test]
    fn test_centered_deviations_branch_cut() {
        // Angles near +π and −π
        let angles = vec![PI - 0.1, -PI + 0.1, PI - 0.05, -PI + 0.05];
        let devs = centered_deviations(&angles);

        // All deviations should be small (near the mean ±π)
        for &d in &devs {
            assert!(
                d.abs() < 0.15,
                "Deviation {d} too large for branch-cut case"
            );
        }
    }

    #[test]
    fn test_empty_input() {
        assert_eq!(circular_mean(&[]), 0.0);
        assert_eq!(mean_resultant_length(&[]), 0.0);
        assert_eq!(circular_variance(&[]), 1.0);
        assert!(centered_deviations(&[]).is_empty());
    }
}

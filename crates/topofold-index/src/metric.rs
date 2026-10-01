//! Metric space abstractions for curve comparison.

use topofold_core::distance::{
    discrete_frechet_invariants, dtw_invariants, writhe_spectrum_distance,
};
use topofold_core::invariants::CurveInvariants;

/// Abstraction of a metric space $(X, d)$ satisfying the triangle inequality:
/// $d(x, z) \le d(x, y) + d(y, z)$.
pub trait Metric<T> {
    /// Computes the distance between elements `a` and `b`.
    fn distance(&self, a: &T, b: &T) -> f64;
}

/// $L_2$ Euclidean metric on localized Writhe spectrum vectors.
#[derive(Debug, Clone, Default)]
pub struct WritheSpectrumMetric;

impl Metric<Vec<f64>> for WritheSpectrumMetric {
    #[inline]
    fn distance(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        writhe_spectrum_distance(a, b)
    }
}

/// Discrete Fréchet metric on $(\kappa, \tau)$ curve invariants.
#[derive(Debug, Clone)]
pub struct FrechetInvariantMetric {
    /// Relative weight for curvature difference.
    pub weight_kappa: f64,
    /// Relative weight for torsion difference.
    pub weight_tau: f64,
}

impl Default for FrechetInvariantMetric {
    fn default() -> Self {
        Self {
            weight_kappa: 1.0,
            weight_tau: 1.0,
        }
    }
}

impl Metric<CurveInvariants> for FrechetInvariantMetric {
    #[inline]
    fn distance(&self, a: &CurveInvariants, b: &CurveInvariants) -> f64 {
        discrete_frechet_invariants(a, b, self.weight_kappa, self.weight_tau)
    }
}

/// Dynamic Time Warping (DTW) alignment metric on $(\kappa, \tau)$ curve invariants.
#[derive(Debug, Clone)]
pub struct DtwInvariantMetric {
    /// Relative weight for curvature difference.
    pub weight_kappa: f64,
    /// Relative weight for torsion difference.
    pub weight_tau: f64,
}

impl Default for DtwInvariantMetric {
    fn default() -> Self {
        Self {
            weight_kappa: 1.0,
            weight_tau: 1.0,
        }
    }
}

impl Metric<CurveInvariants> for DtwInvariantMetric {
    #[inline]
    fn distance(&self, a: &CurveInvariants, b: &CurveInvariants) -> f64 {
        dtw_invariants(a, b, self.weight_kappa, self.weight_tau)
    }
}

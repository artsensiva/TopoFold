//! Autonomous Detection of Cryptic Pockets and Bistable Loops via Intrinsic Geometry and Sarle's Bimodality.
//!
//! Provides unsupervised, alignment-free identification of functionally mobile segments
//! along a protein backbone by evaluating Sarle's Bimodality Coefficient (BC) over sliding-window
//! discrete differential geometry invariants ($\kappa, \tau$).
//!
//! # Mathematical Foundations
//!
//! For an ensemble of $F$ trajectory frames, a subcurve window of length $W$ is slid along
//! the backbone C-alpha trace ($i \dots i + W - 1$). For each window position, the distribution
//! of the mean discrete curvature $\bar{\kappa}$ and torsion $\bar{\tau}$ across all $F$ frames is evaluated.
//!
//! ## Sarle's Bimodality Coefficient (BC)
//!
//! For a 1D sample $X = \{x_1, \dots, x_F\}$ of local curve invariant values:
//!
//! $$BC = \frac{\gamma^2 + 1}{\kappa_{\text{kurt}} + 3 \cdot \frac{(n - 1)^2}{(n - 2)(n - 3)}}$$
//!
//! where:
//! - $\gamma$ is the sample skewness: $\gamma = m_3 / m_2^{1.5}$
//! - $\kappa_{\text{kurt}}$ is the sample excess kurtosis: $\kappa_{\text{kurt}} = m_4 / m_2^2 - 3.0$
//! - $m_r = \frac{1}{n} \sum_{i=1}^n (x_i - \bar{x})^r$ are the sample central moments.
//!
//! ### Biophysical Interpretation:
//! - **Unimodal regimes ($BC < 0.555$)**: Rigid secondary structures (e.g. stable $\alpha$-helices,
//!   $\beta$-sheets) or Gaussian harmonic thermal vibrations yield $BC \approx 0.333$ (normal) to $0.555$ (uniform).
//! - **Bistable regimes ($BC \gg 0.555$, up to $1.0$)**: Functional hinge-bending loops, cryptic pocket openings,
//!   or allosteric rotamer flips switch between discrete metastable states, producing bimodal distributions.
//!
//! ## Streaming Calculation
//!
//! Moments $M_1, M_2, M_3, M_4$ are updated in a single pass using Pébay's (2008) extension of
//! Welford's algorithm in $O(1)$ memory.

use rayon::prelude::*;

use crate::circular::circular_mean;
use crate::error::GeometryError;
use crate::geometry::extract_curve_invariants;
use crate::invariants::CurveInvariants;
use crate::types::BackboneTrace;

/// Streaming online calculation of sample central moments ($M_1, M_2, M_3, M_4$)
/// using Pébay's (2008) extension of Welford's algorithm.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamingMoments {
    n: usize,
    m1: f64, // mean
    m2: f64, // sum of (x - mean)^2
    m3: f64, // sum of (x - mean)^3
    m4: f64, // sum of (x - mean)^4
}

impl Default for StreamingMoments {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingMoments {
    /// Creates a new empty `StreamingMoments` accumulator.
    #[must_use]
    pub fn new() -> Self {
        Self {
            n: 0,
            m1: 0.0,
            m2: 0.0,
            m3: 0.0,
            m4: 0.0,
        }
    }

    /// Feeds a new observation into the accumulator in $O(1)$ time and $O(1)$ space.
    pub fn update(&mut self, x: f64) {
        let n1 = self.n as f64;
        self.n += 1;
        let n = self.n as f64;

        let delta = x - self.m1;
        let delta_n = delta / n;
        let delta_n2 = delta_n * delta_n;
        let term1 = delta * delta_n * n1;

        self.m4 += term1 * delta_n2 * (n * n - 3.0 * n + 3.0) + 6.0 * delta_n2 * self.m2
            - 4.0 * delta_n * self.m3;
        self.m3 += term1 * delta_n * (n - 2.0) - 3.0 * delta_n * self.m2;
        self.m2 += term1;
        self.m1 += delta_n;
    }

    /// Returns the number of observations processed so far.
    #[inline]
    #[must_use]
    pub fn count(&self) -> usize {
        self.n
    }

    /// Returns the sample mean $\bar{x}$.
    #[inline]
    #[must_use]
    pub fn mean(&self) -> f64 {
        self.m1
    }

    /// Returns the sample central variance $m_2 = \frac{1}{n} \sum (x_i - \bar{x})^2$.
    #[inline]
    #[must_use]
    pub fn variance(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.m2 / (self.n as f64)
        }
    }

    /// Returns the sample skewness $\gamma = m_3 / m_2^{1.5}$.
    #[must_use]
    pub fn skewness(&self) -> f64 {
        if self.n < 3 {
            return 0.0;
        }
        let m2 = self.variance();
        if m2 <= 1e-14 {
            return 0.0;
        }
        let m3 = self.m3 / (self.n as f64);
        m3 / m2.powf(1.5)
    }

    /// Returns the sample excess kurtosis $\kappa_{\text{kurt}} = m_4 / m_2^2 - 3.0$.
    #[must_use]
    pub fn excess_kurtosis(&self) -> f64 {
        if self.n < 4 {
            return 0.0;
        }
        let m2 = self.variance();
        if m2 <= 1e-14 {
            return 0.0;
        }
        let m4 = self.m4 / (self.n as f64);
        (m4 / (m2 * m2)) - 3.0
    }

    /// Computes Sarle's Bimodality Coefficient (BC):
    ///
    /// $$BC = \frac{\gamma^2 + 1}{\kappa_{\text{kurt}} + 3 \cdot \frac{(n - 1)^2}{(n - 2)(n - 3)}}$$
    ///
    /// Values range in $[0.0, 1.0]$:
    /// - Normal distribution: $BC \approx 0.333$
    /// - Uniform distribution: $BC \approx 0.555$
    /// - Symmetric bimodal distribution: $BC \to 1.0$
    #[must_use]
    pub fn bimodality_coefficient(&self) -> f64 {
        if self.n < 4 {
            return 0.0;
        }
        let m2 = self.variance();
        if m2 <= 1e-14 {
            return 0.0;
        }

        let gamma = self.skewness();
        let kappa_kurt = self.excess_kurtosis();

        let n = self.n as f64;
        let c = 3.0 * ((n - 1.0) * (n - 1.0)) / ((n - 2.0) * (n - 3.0));
        let denom = kappa_kurt + c;

        if denom <= 1e-12 {
            return 0.0;
        }

        let bc = (gamma * gamma + 1.0) / denom;
        bc.clamp(0.0, 1.0)
    }
}

/// Computes Sarle's Bimodality Coefficient for an arbitrary slice of `f64` values.
///
/// # Examples
/// ```
/// use topofold_core::sarles_bimodality_coefficient;
///
/// // Clean bimodal distribution with sufficient sample size
/// let mut bimodal = Vec::new();
/// for _ in 0..50 {
///     bimodal.push(-2.0);
///     bimodal.push(2.0);
/// }
/// let bc = sarles_bimodality_coefficient(&bimodal);
/// assert!(bc > 0.85);
/// ```
#[must_use]
pub fn sarles_bimodality_coefficient(values: &[f64]) -> f64 {
    let mut moments = StreamingMoments::new();
    for &val in values {
        moments.update(val);
    }
    moments.bimodality_coefficient()
}

/// Evaluated bimodality metrics for a single sliding subcurve window.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowBimodality {
    /// Starting residue index of the window (0-based, inclusive).
    pub window_start: usize,
    /// Ending residue index of the window (0-based, inclusive).
    pub window_end: usize,
    /// Sarle's bimodality coefficient for mean discrete torsion $\tau$.
    pub bc_tau: f64,
    /// Sarle's bimodality coefficient for mean discrete curvature $\kappa$.
    pub bc_kappa: f64,
    /// Sarle's bimodality coefficient for mean side-chain orientation dihedral $\theta_\beta$.
    pub bc_theta: f64,
    /// Composite bimodality score: $\max(BC_\tau, BC_\kappa, BC_\theta)$.
    pub score: f64,
}

/// A detected candidate bistable cryptic pocket or functional loop segment.
#[derive(Debug, Clone, PartialEq)]
pub struct PocketCandidate {
    /// Starting residue index of the detected candidate cluster (0-based, inclusive).
    pub start_res: usize,
    /// Ending residue index of the detected candidate cluster (0-based, inclusive).
    pub end_res: usize,
    /// Peak bimodality score within this candidate cluster.
    pub score: f64,
    /// Sarle's bimodality coefficient for discrete torsion at the peak window.
    pub bc_tau: f64,
    /// Sarle's bimodality coefficient for discrete curvature at the peak window.
    pub bc_kappa: f64,
    /// Sarle's bimodality coefficient for side-chain orientation dihedral at the peak window.
    pub bc_theta: f64,
    /// Starting residue index of the peak window within the cluster.
    pub peak_res: usize,
}

/// Computes the sliding-window bimodality coefficient profile across an ensemble of traces.
///
/// For each window position $i$ of length `window_size` (0-indexed $i \dots i + W - 1$),
/// computes the distribution of the mean discrete curvature and torsion across all frames,
/// evaluating Sarle's Bimodality Coefficient.
///
/// # Errors
/// Returns [`GeometryError`] if trajectory traces have invalid geometry or mismatched lengths.
pub fn compute_bimodality_profile(
    traces: &[BackboneTrace],
    window_size: usize,
) -> Result<Vec<WindowBimodality>, GeometryError> {
    if traces.len() < 4 {
        return Ok(Vec::new());
    }

    let n_residues = traces[0].len();
    if n_residues < window_size || window_size < 4 {
        return Ok(Vec::new());
    }

    for t in traces {
        if t.len() != n_residues {
            return Err(GeometryError::InsufficientPoints {
                required: n_residues,
                actual: t.len(),
            });
        }
    }

    // Precompute curve invariants for all frames in parallel
    let invariants: Result<Vec<CurveInvariants>, GeometryError> =
        traces.par_iter().map(extract_curve_invariants).collect();
    let invariants = invariants?;

    let n_windows = n_residues - window_size + 1;
    let mut profile = Vec::with_capacity(n_windows);

    let n_kappas = window_size - 2;
    let n_torsions = window_size - 3;

    for i in 0..n_windows {
        let window_start = i;
        let window_end = i + window_size - 1;

        let has_theta = !invariants.is_empty() && !invariants[0].sidechain_dihedrals.is_empty();

        // Collect circular means across all frames for this window
        let mut tau_means = Vec::with_capacity(invariants.len());
        let mut kappa_means = Vec::with_capacity(invariants.len());
        let mut theta_means = Vec::with_capacity(invariants.len());

        for inv in &invariants {
            // Mean discrete curvature across the window
            let k_sum: f64 = inv.curvatures[i..i + n_kappas].iter().sum();
            let k_mean = k_sum / (n_kappas as f64);
            kappa_means.push(k_mean);

            // Mean discrete torsion across the window (circular mean for angular data)
            let t_slice = &inv.torsions[i..i + n_torsions];
            let t_mean = circular_mean(t_slice);
            tau_means.push(t_mean);

            // Mean side-chain orientation dihedral across the window (if available, circular mean)
            if has_theta && inv.sidechain_dihedrals.len() >= i + n_kappas {
                let th_slice = &inv.sidechain_dihedrals[i..i + n_kappas];
                let th_mean = circular_mean(th_slice);
                theta_means.push(th_mean);
            }
        }

        // For angular quantities (tau, theta), compute centered deviations from circular mean
        // before computing moments (prevents branch-cut artifacts)
        use crate::circular::centered_deviations;
        let tau_deviations = centered_deviations(&tau_means);
        let theta_deviations = if has_theta {
            centered_deviations(&theta_means)
        } else {
            Vec::new()
        };

        // Compute moments from the centered deviations
        let mut moments_tau = StreamingMoments::new();
        let mut moments_kappa = StreamingMoments::new();
        let mut moments_theta = StreamingMoments::new();

        for &k in &kappa_means {
            moments_kappa.update(k);
        }
        for &t_dev in &tau_deviations {
            moments_tau.update(t_dev);
        }
        if has_theta {
            for &th_dev in &theta_deviations {
                moments_theta.update(th_dev);
            }
        }

        let bc_tau = moments_tau.bimodality_coefficient();
        let bc_kappa = moments_kappa.bimodality_coefficient();
        let bc_theta = if has_theta {
            moments_theta.bimodality_coefficient()
        } else {
            0.0
        };
        let score = bc_tau.max(bc_kappa).max(bc_theta);

        profile.push(WindowBimodality {
            window_start,
            window_end,
            bc_tau,
            bc_kappa,
            bc_theta,
            score,
        });
    }

    Ok(profile)
}

/// Detects bistable functional loop segments and cryptic pockets without prior structural knowledge,
/// allowing a configurable gap tolerance between adjacent windows.
///
/// # Arguments
/// - `traces`: Ensemble of backbone coordinate traces across the trajectory.
/// - `window_size`: Length of the sliding subcurve window (default: 8 residues).
/// - `bc_threshold`: Minimum Sarle's Bimodality Coefficient threshold (e.g. 0.60).
/// - `max_gap`: Maximum allowed residue gap between windows to merge into a continuous cluster (default: 1).
///
/// # Returns
/// A list of [`PocketCandidate`] segments sorted in descending order of transition score.
pub fn detect_bistable_segments_with_gap(
    traces: &[BackboneTrace],
    window_size: usize,
    bc_threshold: f64,
    max_gap: usize,
) -> Vec<PocketCandidate> {
    let profile = match compute_bimodality_profile(traces, window_size) {
        Ok(p) => p,
        Err(_) => return Vec::new(),
    };

    if profile.is_empty() {
        return Vec::new();
    }

    // Collect all window indices where bimodality exceeds threshold
    let above_threshold_indices: Vec<usize> = profile
        .iter()
        .enumerate()
        .filter(|(_, w)| w.score >= bc_threshold)
        .map(|(i, _)| i)
        .collect();

    if above_threshold_indices.is_empty() {
        return Vec::new();
    }

    // Cluster consecutive or near-consecutive windows (bridging small gaps <= max_gap)
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut current_cluster = vec![above_threshold_indices[0]];

    for &idx in &above_threshold_indices[1..] {
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

    let mut candidates = Vec::with_capacity(clusters.len());

    for cluster in clusters {
        let first_win = cluster[0];
        let last_win = *cluster.last().unwrap();

        let start_res = profile[first_win].window_start;
        let end_res = profile[last_win].window_end;

        let mut peak_score = 0.0;
        let mut peak_tau = 0.0;
        let mut peak_kappa = 0.0;
        let mut peak_theta = 0.0;
        let mut peak_res = start_res;

        for &w_idx in &cluster {
            let w = &profile[w_idx];
            if w.score > peak_score {
                peak_score = w.score;
                peak_tau = w.bc_tau;
                peak_kappa = w.bc_kappa;
                peak_theta = w.bc_theta;
                peak_res = w.window_start;
            }
        }

        candidates.push(PocketCandidate {
            start_res,
            end_res,
            score: peak_score,
            bc_tau: peak_tau,
            bc_kappa: peak_kappa,
            bc_theta: peak_theta,
            peak_res,
        });
    }

    // Rank candidate segments by peak bimodality score descending
    candidates.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    candidates
}

/// Detects bistable functional loop segments and cryptic pockets without prior structural knowledge.
///
/// Slides a window of length `window_size` (default: 8 residues) along the backbone trace,
/// evaluates Sarle's Bimodality Coefficient on discrete curve invariants across the trajectory,
/// identifies continuous residue clusters where $BC \ge \text{bc\_threshold}$ (e.g. 0.60),
/// and ranks candidate segments by peak bimodality score.
///
/// # Arguments
/// - `traces`: Ensemble of backbone coordinate traces across the trajectory.
/// - `window_size`: Length of the sliding subcurve window in residues (typically 8).
/// - `bc_threshold`: Minimum Sarle's Bimodality Coefficient (typically 0.60).
///
/// # Returns
/// A list of [`PocketCandidate`] segments sorted in descending order of transition score.
#[must_use]
pub fn detect_bistable_segments(
    traces: &[BackboneTrace],
    window_size: usize,
    bc_threshold: f64,
) -> Vec<PocketCandidate> {
    detect_bistable_segments_with_gap(traces, window_size, bc_threshold, 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    #[test]
    fn test_streaming_moments_normal_distribution() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut moments = StreamingMoments::new();

        // Standard normal distribution: skewness ~ 0, excess kurtosis ~ 0, BC ~ 0.333
        for _ in 0..10_000 {
            // Box-Muller transform
            let u1: f64 = rng.gen();
            let u2: f64 = rng.gen();
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            moments.update(z);
        }

        assert_eq!(moments.count(), 10_000);
        assert_relative_eq!(moments.mean(), 0.0, epsilon = 0.05);
        assert_relative_eq!(moments.variance(), 1.0, epsilon = 0.05);
        assert_relative_eq!(moments.skewness(), 0.0, epsilon = 0.05);
        assert_relative_eq!(moments.excess_kurtosis(), 0.0, epsilon = 0.15);

        let bc = moments.bimodality_coefficient();
        assert!(
            bc < 0.555,
            "Normal distribution must be unimodal (< 0.555), got {bc}"
        );
        assert_relative_eq!(bc, 0.333, epsilon = 0.05);
    }

    #[test]
    fn test_streaming_moments_uniform_distribution() {
        let mut rng = StdRng::seed_from_u64(12345);
        let mut moments = StreamingMoments::new();

        // Uniform distribution: skewness ~ 0, excess kurtosis ~ -1.2, BC ~ 0.555
        for _ in 0..10_000 {
            let u: f64 = rng.gen_range(-1.0..1.0);
            moments.update(u);
        }

        let bc = moments.bimodality_coefficient();
        assert_relative_eq!(bc, 0.555, epsilon = 0.05);
    }

    #[test]
    fn test_streaming_moments_bimodal_distribution() {
        let mut moments = StreamingMoments::new();

        // 50% at -3.0 and 50% at +3.0
        for _ in 0..500 {
            moments.update(-3.0);
            moments.update(3.0);
        }

        let bc = moments.bimodality_coefficient();
        assert!(
            bc > 0.95,
            "Clean bimodal switch must have BC approaching 1.0, got {bc}"
        );
    }

    #[test]
    fn test_edge_cases_and_zero_variance() {
        assert_eq!(sarles_bimodality_coefficient(&[]), 0.0);
        assert_eq!(sarles_bimodality_coefficient(&[1.0, 2.0]), 0.0);
        assert_eq!(sarles_bimodality_coefficient(&[5.0; 100]), 0.0);
    }

    #[test]
    fn test_detect_bistable_synthetic_loop() {
        // Construct a 20-residue chain across 100 frames.
        // Residues 6..12 undergo a bistable flip:
        // Frames 0..49: State A (straight line in x)
        // Frames 50..99: State B (arch in y)
        let n_frames = 100;
        let n_residues = 20;
        let mut traces = Vec::with_capacity(n_frames);

        for f in 0..n_frames {
            let mut coords = Vec::with_capacity(n_residues);
            for i in 0..n_residues {
                let angle = i as f64 * 1.5;
                let x = i as f64 * 2.5;
                let mut y = angle.sin() * 2.0;
                let mut z = angle.cos() * 2.0;

                if f >= 50 && (6..=12).contains(&i) {
                    // Bistable flip in residues 6..12
                    y += 4.0 * ((i - 6) as f64 * std::f64::consts::PI / 6.0).sin();
                    z += 3.0 * ((i - 6) as f64 * std::f64::consts::PI / 6.0).sin();
                }

                coords.push([x, y, z]);
            }
            traces.push(BackboneTrace::from_arrays(&coords));
        }

        let candidates = detect_bistable_segments(&traces, 8, 0.6);
        assert!(
            !candidates.is_empty(),
            "Should detect at least 1 bistable candidate"
        );

        let best = &candidates[0];
        assert!(best.score > 0.6, "Candidate score should exceed 0.6");
        // Ensure detected range overlaps with the synthetic flip region [6, 12]
        assert!(
            best.start_res <= 12 && best.end_res >= 6,
            "Detected candidate [{}, {}] must overlap flip region [6, 12]",
            best.start_res,
            best.end_res
        );
    }

    #[test]
    fn test_detect_bistable_sidechain_rotamer_gating() {
        // Backbone C-alpha coordinates are completely STATIC across 100 frames (zero Cartesian/kappa/tau shift).
        // Residues 7..12 undergo a pure rotameric side-chain gating flip:
        // Frames 0..49: C-beta oriented at theta ~ 0
        // Frames 50..99: C-beta oriented at theta ~ pi/2
        let n_frames = 100;
        let n_residues = 20;
        let mut traces = Vec::with_capacity(n_frames);

        for f in 0..n_frames {
            let mut ca_coords = Vec::with_capacity(n_residues);
            let mut cb_coords = Vec::with_capacity(n_residues);

            for i in 0..n_residues {
                let angle = i as f64 * 0.8;
                let ca = [i as f64 * 2.0, angle.sin() * 2.0, angle.cos() * 2.0];
                ca_coords.push(ca);

                // Default side-chain orientation
                let mut cb = [ca[0], ca[1] + 1.5, ca[2]];
                if f >= 50 && (7..=12).contains(&i) {
                    // Rotameric gating flip: side-chain swings out of plane
                    cb = [ca[0], ca[1], ca[2] + 1.5];
                }
                cb_coords.push(cb);
            }
            traces.push(BackboneTrace::from_arrays_with_cbeta(
                &ca_coords, &cb_coords,
            ));
        }

        let candidates = detect_bistable_segments(&traces, 8, 0.6);
        assert!(
            !candidates.is_empty(),
            "Rotameric side-chain gating must be detected even with static backbone"
        );

        let best = &candidates[0];
        assert!(
            best.bc_theta > 0.8,
            "bc_theta should reflect sharp bimodal rotamer switch, got {}",
            best.bc_theta
        );
        assert!(
            best.start_res <= 12 && best.end_res >= 7,
            "Detected segment [{}, {}] must capture rotamer gating region [7, 12]",
            best.start_res,
            best.end_res
        );
    }

    #[test]
    fn test_circular_mean_branch_cut_regression() {
        // M2 regression test: torsion angles fluctuating unimodally around ±π
        // should NOT produce false bimodality (BC > 0.9) with circular mean.
        // This was the bug in v0.8.3 where linear averaging caused BC 0.96-0.98
        // for β-strand regions with τ ≈ -170° (see benchmarks/audit/bc_check.py).
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let mut rng = StdRng::seed_from_u64(999);
        let n_frames = 100;
        let n_residues = 15;
        let mut traces = Vec::with_capacity(n_frames);

        for _ in 0..n_frames {
            let mut coords = Vec::with_capacity(n_residues);
            for i in 0..n_residues {
                // Build a β-strand-like structure with torsions near ±180°
                let x = i as f64 * 3.8; // Extended β-strand spacing
                let y = (i as f64 * std::f64::consts::PI).sin() * 0.5; // Slight pleating
                let z = (i as f64 * std::f64::consts::PI).cos() * 0.5;

                // Add small unimodal thermal fluctuations (±0.1 radians ≈ ±6°)
                let noise_y = rng.gen_range(-0.1..0.1);
                let noise_z = rng.gen_range(-0.1..0.1);

                coords.push([x, y + noise_y, z + noise_z]);
            }
            traces.push(BackboneTrace::from_arrays(&coords));
        }

        // Extract invariants to verify that torsions are indeed near ±π
        let invariants: Vec<CurveInvariants> = traces
            .iter()
            .map(|trace| extract_curve_invariants(trace).unwrap())
            .collect();

        // Verify that middle segment torsions are near ±π
        let mid_frame_torsions = &invariants[0].torsions[5..10];
        for &tau in mid_frame_torsions {
            assert!(
                tau.abs() > 2.5,
                "β-strand torsions should be near ±π, got {tau}"
            );
        }

        // Run bimodality scan with window size 8
        let candidates = detect_bistable_segments(&traces, 8, 0.6);

        // With circular mean, unimodal fluctuations near ±π should NOT trigger
        // false bimodality. BC should be < 0.555 (unimodal threshold).
        // With linear mean (the bug), this would produce BC > 0.9.
        if let Some(candidate) = candidates.first() {
            assert!(
                candidate.bc_tau < 0.555,
                "Unimodal fluctuations around ±π must not produce false bimodality. \
                 BC_tau = {}, expected < 0.555 (unimodal). \
                 Circular mean should prevent this false positive.",
                candidate.bc_tau
            );
        }
        // If no candidates found, that's also correct (BC was below threshold)
    }

    #[test]
    fn test_true_bimodal_angular_distribution_preserved() {
        // Positive control: TRUE bimodal angular distribution should still yield high BC.
        // This test ensures we don't suppress legitimate bimodal signals after M2 fix.
        //
        // Simulates two rotamer states (e.g., side-chain χ₁ dihedral):
        // State A (frames 0-49): torsions ~ −60° (gauche-)
        // State B (frames 50-99): torsions ~ +60° (gauche+)
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let mut rng = StdRng::seed_from_u64(777);
        let n_frames = 100;
        let n_residues = 15;
        let mut traces = Vec::with_capacity(n_frames);

        for f in 0..n_frames {
            let mut coords = Vec::with_capacity(n_residues);
            for i in 0..n_residues {
                // Build a structure where residues 6-11 undergo bimodal rotamer switching
                let x = i as f64 * 3.5;
                let base_angle = i as f64 * 0.5;

                // Rotamer state: -60° for first half, +60° for second half
                let rotamer_shift = if f < 50 {
                    std::f64::consts::PI / 3.0 // +60° (state A)
                } else {
                    -std::f64::consts::PI / 3.0 // -60° (state B)
                };

                // Apply rotamer shift to middle residues 6-11
                let angle = if (6..=11).contains(&i) {
                    rotamer_shift + base_angle
                } else {
                    base_angle
                };

                // Small thermal noise
                let noise = rng.gen_range(-0.05..0.05);
                let y = (angle + noise).sin() * 2.5;
                let z = (angle + noise).cos() * 2.5;

                coords.push([x, y, z]);
            }
            traces.push(BackboneTrace::from_arrays(&coords));
        }

        // Run bimodality scan
        let candidates = detect_bistable_segments(&traces, 8, 0.6);

        // Should detect the bimodal rotamer switching in residues 6-11
        assert!(
            !candidates.is_empty(),
            "True bimodal angular distribution (rotamer switching) must be detected"
        );

        let best = &candidates[0];
        assert!(
            best.bc_tau > 0.7,
            "True bimodal rotamer switch should yield high BC_tau. Got {}, expected > 0.7. \
             Circular statistics should preserve legitimate bimodal signals.",
            best.bc_tau
        );

        // Verify detected region overlaps with the true switching region [6, 11]
        assert!(
            best.start_res <= 11 && best.end_res >= 6,
            "Detected segment [{}, {}] must overlap true bimodal region [6, 11]",
            best.start_res,
            best.end_res
        );
    }
}

//! Two-tier conformational search engine for protein trajectories.
//!
//! Tier 1: Rapid coarse filtering using a Vantage-Point Tree on localized Writhe spectra ($L_2$).
//! Tier 2: Precise fine verification using the Discrete Fréchet distance on $(\kappa, \tau)$ curve invariants.

#![forbid(unsafe_code)]

use topofold_core::distance::discrete_frechet_invariants;
use topofold_core::invariants::CurveInvariants;

use crate::metric::{Metric, WritheSpectrumMetric};
use crate::vptree::VpTree;

/// An indexed conformational frame from an MD trajectory or ensemble.
#[derive(Debug, Clone)]
pub struct ConformerFrame {
    /// Frame sequence index or unique identifier.
    pub frame_id: u64,
    /// Simulation timestamp in picoseconds.
    pub time_ps: f64,
    /// Complete discrete differential curve invariants $(\kappa, \tau, l)$.
    pub invariants: CurveInvariants,
    /// Localized Writhe spectrum vector for coarse metric filtering.
    pub writhe_spectrum: Vec<f64>,
}

/// Helper struct for VP-Tree indexing over writhe spectrum.
#[derive(Debug, Clone)]
struct IndexedWrithe {
    frame_idx: usize,
    spectrum: Vec<f64>,
}

impl Metric<IndexedWrithe> for WritheSpectrumMetric {
    #[inline]
    fn distance(&self, a: &IndexedWrithe, b: &IndexedWrithe) -> f64 {
        self.distance(&a.spectrum, &b.spectrum)
    }
}

/// Result returned from a conformational search query.
#[derive(Debug, Clone, PartialEq)]
pub struct ConformerHit {
    /// Frame identifier.
    pub frame_id: u64,
    /// Timestamp in picoseconds.
    pub time_ps: f64,
    /// Measured Discrete Fréchet distance on the $(\kappa, \tau)$ manifold.
    pub frechet_distance: f64,
    /// Measured coarse Writhe spectrum $L_2$ distance.
    pub writhe_distance: f64,
}

/// High-performance two-tier conformational index.
#[derive(Debug)]
pub struct ConformationalIndex {
    frames: Vec<ConformerFrame>,
    vp_tree: VpTree<IndexedWrithe, WritheSpectrumMetric>,
}

impl ConformationalIndex {
    /// Constructs a `ConformationalIndex` from a collection of conformer frames.
    #[must_use]
    pub fn build(frames: Vec<ConformerFrame>) -> Self {
        let indexed_items: Vec<IndexedWrithe> = frames
            .iter()
            .enumerate()
            .map(|(frame_idx, f)| IndexedWrithe {
                frame_idx,
                spectrum: f.writhe_spectrum.clone(),
            })
            .collect();

        let vp_tree = VpTree::new(indexed_items, WritheSpectrumMetric);
        Self { frames, vp_tree }
    }

    /// Number of frames indexed.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Returns `true` if the index contains no frames.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Access to an indexed frame by internal index.
    #[inline]
    #[must_use]
    pub fn get_frame(&self, index: usize) -> Option<&ConformerFrame> {
        self.frames.get(index)
    }

    /// Executes a two-tier search for conformers matching `target_invariants` within `frechet_tol`.
    ///
    /// - Tier 1: Filters candidates whose Writhe spectrum distance is within `writhe_tol`.
    /// - Tier 2: Validates candidates against exact Discrete Fréchet distance on $(\kappa, \tau)$.
    #[must_use]
    pub fn query_conformations(
        &self,
        target_invariants: &CurveInvariants,
        target_writhe: &[f64],
        frechet_tol: f64,
        writhe_tol: f64,
    ) -> Vec<ConformerHit> {
        let dummy_query = IndexedWrithe {
            frame_idx: 0,
            spectrum: target_writhe.to_vec(),
        };

        // Tier 1: Sub-linear VP-tree range search on Writhe spectrum
        let coarse_hits = self.vp_tree.range_search(&dummy_query, writhe_tol);

        let mut hits = Vec::new();

        // Tier 2: Exact Discrete Fréchet verification on filtered candidate conformers
        for (item, w_dist) in coarse_hits {
            let frame = &self.frames[item.frame_idx];
            let f_dist = discrete_frechet_invariants(
                &frame.invariants,
                target_invariants,
                1.0,
                1.0,
            );

            if f_dist <= frechet_tol {
                hits.push(ConformerHit {
                    frame_id: frame.frame_id,
                    time_ps: frame.time_ps,
                    frechet_distance: f_dist,
                    writhe_distance: w_dist,
                });
            }
        }

        // Sort by finest metric (Fréchet distance)
        hits.sort_by(|a, b| {
            a.frechet_distance
                .partial_cmp(&b.frechet_distance)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        hits
    }
}

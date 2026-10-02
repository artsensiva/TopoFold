//! # TopoFold Index
//!
//! Sub-linear metric space indexing and search engine for protein conformational
//! trajectories, using Vantage-Point Trees (VP-Trees) over discrete differential
//! curve invariants and Writhe spectra.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod engine;
pub mod metric;
pub mod vptree;

pub use engine::{ConformationalIndex, ConformerFrame, ConformerHit};
pub use metric::{DtwInvariantMetric, FrechetInvariantMetric, Metric, WritheSpectrumMetric};
pub use vptree::{VpNode, VpTree};

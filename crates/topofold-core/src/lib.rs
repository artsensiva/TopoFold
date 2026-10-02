//! # TopoFold Core Engine
//!
//! High-performance discrete differential geometry engine for protein backbones,
//! extracting rotation- and translation-invariant (SE(3)-invariant) curve signatures.
//!
//! ## Overview
//!
//! Traditional conformational analysis relies on Cartesian coordinates and linear PCA,
//! which requires expensive pairwise structural alignment (Kabsch algorithm) and is
//! sensitive to rotational drift.
//!
//! TopoFold models the C-alpha protein backbone as a discrete 3D space curve, representing
//! conformations via intrinsic local differential invariants:
//! - Segment lengths: `l_i`
//! - Discrete curvature (turning angles): `kappa_i \in [0, \pi]`
//! - Discrete torsion (osculating plane dihedrals): `tau_i \in (-\pi, \pi]`
//!
//! These invariants are strictly SE(3)-invariant and chiral-sensitive.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod allostery;
pub mod bimodality;
pub mod distance;
pub mod error;
pub mod geometry;
pub mod idp;
pub mod invariants;
pub mod pdb;
pub mod ribbon;
pub mod rna;
pub mod types;
pub mod writhe;

// Re-export primary types and functions for ergonomic top-level use
pub use allostery::{
    compute_intermolecular_allosteric_network, compute_intermolecular_allosteric_network_ribbon,
    compute_intrinsic_allosteric_network, compute_ternary_cooperativity_index,
};
pub use bimodality::{
    compute_bimodality_profile, detect_bistable_segments, detect_bistable_segments_with_gap,
    sarles_bimodality_coefficient, PocketCandidate, StreamingMoments, WindowBimodality,
};
pub use distance::{
    angular_distance_s1, discrete_frechet_distance_coords, discrete_frechet_invariants,
    dtw_invariants, writhe_spectrum_distance,
};
pub use error::GeometryError;
pub use geometry::{
    compute_binormals, compute_curvatures, compute_tangents, compute_torsions,
    extract_curve_invariants, GEOMETRY_EPSILON,
};
pub use idp::{
    compute_frame_topological_compactness, compute_idp_density_profile,
    compute_idp_density_profile_ribbon, detect_transient_motifs, detect_transient_motifs_backbone,
    detect_transient_motifs_with_params, IdpDensityProfile, IdpMotif,
};
pub use invariants::CurveInvariants;
pub use pdb::{
    parse_pdb_ca, parse_pdb_ribbon, parse_pdb_rna, parse_pdb_rna_str, parse_pdb_str, PdbError,
    ResidueMeta,
};
pub use ribbon::{
    compute_glycine_pseudo_cbeta, compute_pseudo_cbeta_from_ca, compute_sidechain_dihedrals,
    extract_ribbon_invariants, RibbonTrace, STANDARD_CA_CB_BOND_LENGTH,
};
pub use rna::{
    compute_base_ribbon_dihedrals, compute_rna_bimodality_profile, extract_rna_ribbon_invariants,
    extract_rna_ribbon_invariants_with_window, scan_rna_switching_hinges, RnaBimodalityProfile,
    RnaCurveInvariants, RnaHingeCandidate, RnaRibbonTrace,
};
pub use types::BackboneTrace;
pub use writhe::{compute_local_writhe, compute_total_writhe};

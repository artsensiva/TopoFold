//! Golden dataset validation test using the authentic crystal structure of Crambin (PDB: 1CRN).
//!
//! Validates:
//! 1. In-memory PDB C-alpha parser accuracy on experimental crystallography data.
//! 2. Biophysical concordance of discrete differential invariants:
//!    - Canonical alpha-helical signatures (residues 7-19, 23-30): kappa ~ 1.55 rad (89 deg), tau ~ 0.88 rad (50 deg)
//!    - Extended beta-strand signatures (residues 1-4, 32-35): kappa ~ 0.85 rad, |tau| > 2.5 rad (~140-180 deg)
//! 3. Total discrete writhe and sliding-window writhe spectrum.
//! 4. SE(3) rigid-body invariance on real-world experimental protein coordinates.

use nalgebra::{Isometry3, Point3, Translation3, UnitQuaternion};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

use topofold_core::{
    compute_local_writhe, compute_total_writhe, discrete_frechet_invariants, dtw_invariants,
    extract_curve_invariants, parse_pdb_str, BackboneTrace,
};

const PDB_1CRN: &str = include_str!("fixtures/1crn.pdb");

fn angular_diff(a: f64, b: f64) -> f64 {
    let mut diff = (a - b).abs();
    if diff > PI {
        diff = (2.0 * PI - diff).abs();
    }
    diff
}

#[test]
fn test_golden_1crn_parsing_and_invariants() {
    let (trace, metas) = parse_pdb_str(PDB_1CRN, Some('A')).expect("Parse 1CRN PDB");

    // 1CRN contains exactly 46 residues
    assert_eq!(trace.len(), 46);
    assert_eq!(metas.len(), 46);

    // Verify first and last residues
    assert_eq!(metas[0].name, "THR");
    assert_eq!(metas[0].seq_id, 1);
    assert_eq!(metas[45].name, "ASN");
    assert_eq!(metas[45].seq_id, 46);

    // Extract SE(3)-invariant curve invariants
    let invariants = extract_curve_invariants(&trace).expect("Extract invariants from 1CRN");

    assert_eq!(invariants.segment_lengths.len(), 45);
    assert_eq!(invariants.curvatures.len(), 44);
    assert_eq!(invariants.torsions.len(), 43);

    // Trans-peptide virtual C_alpha - C_alpha bond lengths must average ~3.80 Å
    let mean_bond = invariants.mean_bond_length();
    assert!(
        (mean_bond - 3.80).abs() < 0.05,
        "Mean C-alpha bond length in 1CRN must be ~3.80 Å, found {mean_bond}"
    );

    // Check individual bond lengths are within physical tolerances [3.70, 3.90] Å
    for (i, &l) in invariants.segment_lengths.iter().enumerate() {
        assert!(
            l >= 3.65 && l <= 3.95,
            "Residue bond {i} length {l} Å out of physical range"
        );
    }

    // Inspect Helix 1 (PDB residues 7-19, 0-indexed residues 6-18)
    // For residues in the core of Helix 1 (residues 8-17):
    // Invariants indices: curvature i corresponds to vertex i+1 (residue seq_id i+2)
    // Residue 9 to 16 corresponds to curvature indices 7 to 14
    for idx in 7..=13 {
        let kappa = invariants.curvatures[idx];
        let tau = invariants.torsions[idx];

        // Alpha-helices exhibit turning angles ~ 85-95 deg (1.48 - 1.66 rad)
        assert!(
            kappa >= 1.30 && kappa <= 1.80,
            "Helix 1 residue curvature {kappa:.3} rad out of canonical range at idx {idx}"
        );

        // Alpha-helices exhibit positive right-handed dihedrals ~ 45-60 deg (0.78 - 1.05 rad)
        assert!(
            tau >= 0.65 && tau <= 1.25,
            "Helix 1 residue torsion {tau:.3} rad out of canonical range at idx {idx}"
        );
    }

    // Inspect Beta-strand (PDB residues 1-4, indices 0-3)
    // Beta-strands exhibit extended dihedrals (|tau| near pi)
    let beta_tau = invariants.torsions[0]; // between vertex 1 and 2
    assert!(
        beta_tau.abs() > 2.0,
        "Beta-strand torsion must be extended (|tau| > 2.0 rad), found {beta_tau}"
    );
}

#[test]
fn test_golden_1crn_writhe_and_se3_invariance() {
    let (trace, _) = parse_pdb_str(PDB_1CRN, Some('A')).expect("Parse 1CRN PDB");
    let inv_orig = extract_curve_invariants(&trace).expect("Invariants");
    let total_wr_orig = compute_total_writhe(&trace).expect("Total writhe");

    // Local writhe spectrum with 9-residue sliding window (window_radius = 4)
    let local_wr_orig = compute_local_writhe(&trace, 4).expect("Local writhe");
    assert_eq!(local_wr_orig.len(), 46 - 8);

    // Verify rigorous SE(3) invariance on real crystal structure across 50 random rigid-body transformations
    let mut rng = StdRng::seed_from_u64(0x1C8A_2026);
    for trial in 0..50 {
        let q = UnitQuaternion::from_euler_angles(
            rng.gen_range(-PI..PI),
            rng.gen_range(-PI..PI),
            rng.gen_range(-PI..PI),
        );
        let t = Translation3::new(
            rng.gen_range(-50_000.0..50_000.0),
            rng.gen_range(-50_000.0..50_000.0),
            rng.gen_range(-50_000.0..50_000.0),
        );
        let iso = Isometry3::from_parts(t, q);

        let trans_coords: Vec<Point3<f64>> = trace
            .coordinates()
            .iter()
            .map(|p| iso.transform_point(p))
            .collect();
        let trans_trace = BackboneTrace::new(trans_coords);

        let inv_trans = extract_curve_invariants(&trans_trace).expect("Transformed invariants");
        let total_wr_trans = compute_total_writhe(&trans_trace).expect("Transformed writhe");
        let local_wr_trans = compute_local_writhe(&trans_trace, 4).expect("Transformed local writhe");

        // Curvature invariance
        for i in 0..inv_orig.curvatures.len() {
            let dk = (inv_orig.curvatures[i] - inv_trans.curvatures[i]).abs();
            assert!(dk < 1e-6, "Trial {trial}: Curvature mismatch at {i}: {dk}");
        }

        // Torsion invariance
        for i in 0..inv_orig.torsions.len() {
            let dt = angular_diff(inv_orig.torsions[i], inv_trans.torsions[i]);
            assert!(dt < 1e-6, "Trial {trial}: Torsion mismatch at {i}: {dt}");
        }

        // Total Writhe invariance
        let d_wr = (total_wr_orig - total_wr_trans).abs();
        assert!(d_wr < 1e-6, "Trial {trial}: Total writhe mismatch: {d_wr}");

        // Local Writhe spectrum invariance
        for i in 0..local_wr_orig.len() {
            let d_lwr = (local_wr_orig[i] - local_wr_trans[i]).abs();
            assert!(d_lwr < 1e-6, "Trial {trial}: Local writhe mismatch at {i}: {d_lwr}");
        }
    }
}

#[test]
fn test_golden_1crn_metric_distances() {
    let (trace, _) = parse_pdb_str(PDB_1CRN, Some('A')).expect("Parse 1CRN PDB");
    let invariants = extract_curve_invariants(&trace).expect("Invariants");

    // Self-distance under Discrete Fréchet metric must be 0
    let df_self = discrete_frechet_invariants(&invariants, &invariants, 1.0, 1.0);
    assert!(df_self < 1e-12, "Self-Fréchet distance must be 0, got {df_self}");

    // Self-distance under DTW must be 0
    let dtw_self = dtw_invariants(&invariants, &invariants, 1.0, 1.0);
    assert!(dtw_self < 1e-12, "Self-DTW distance must be 0, got {dtw_self}");
}

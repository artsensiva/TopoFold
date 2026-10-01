//! Integration tests verifying SE(3) invariance and biophysical curve signatures.

use nalgebra::{Isometry3, Point3, Translation3, UnitQuaternion};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;
use topofold_core::{extract_curve_invariants, BackboneTrace};

fn angular_diff(a: f64, b: f64) -> f64 {
    let mut diff = (a - b).abs();
    if diff > PI {
        diff = (2.0 * PI - diff).abs();
    }
    diff
}

#[test]
fn test_integration_se3_random_walk_invariance() {
    let mut rng = StdRng::seed_from_u64(0xCAFE_BABE);
    let n_points = 100;

    // Generate a pseudo-random C-alpha chain with realistic 3.8 Å bond lengths
    let mut coords = Vec::with_capacity(n_points);
    coords.push(Point3::new(0.0, 0.0, 0.0));

    for i in 1..n_points {
        let prev = coords[i - 1];
        // Sample direction on S^2
        let phi: f64 = rng.gen_range(0.0..2.0 * PI);
        let cos_theta: f64 = rng.gen_range(-0.8..0.8); // avoid extreme folding
        let sin_theta: f64 = (1.0 - cos_theta * cos_theta).sqrt();

        let step = nalgebra::Vector3::new(
            sin_theta * phi.cos(),
            sin_theta * phi.sin(),
            cos_theta,
        ) * 3.80; // Trans-peptide C_alpha - C_alpha distance

        coords.push(prev + step);
    }

    let trace_orig = BackboneTrace::new(coords);
    let inv_orig = extract_curve_invariants(&trace_orig).expect("Random walk invariants");

    for trial in 0..50 {
        let q = UnitQuaternion::from_euler_angles(
            rng.gen_range(-PI..PI),
            rng.gen_range(-PI..PI),
            rng.gen_range(-PI..PI),
        );
        let t = Translation3::new(
            rng.gen_range(-500.0..500.0),
            rng.gen_range(-500.0..500.0),
            rng.gen_range(-500.0..500.0),
        );
        let iso = Isometry3::from_parts(t, q);

        let trans_coords: Vec<Point3<f64>> = trace_orig
            .coordinates()
            .iter()
            .map(|p| iso.transform_point(p))
            .collect();

        let trace_trans = BackboneTrace::new(trans_coords);
        let inv_trans = extract_curve_invariants(&trace_trans).expect("Transformed invariants");

        for i in 0..inv_orig.curvatures.len() {
            let dk = (inv_orig.curvatures[i] - inv_trans.curvatures[i]).abs();
            assert!(
                dk < 1e-6,
                "Trial {trial}: Curvature mismatch at {i}: {dk}"
            );
        }

        for i in 0..inv_orig.torsions.len() {
            let dt = angular_diff(inv_orig.torsions[i], inv_trans.torsions[i]);
            assert!(
                dt < 1e-6,
                "Trial {trial}: Torsion mismatch at {i}: {dt}"
            );
        }
    }

    // Verify subcurve extraction exactness:
    // Extracting invariants of trace[20..=60] must match inv_orig.subcurve(20, 60)
    let start = 20;
    let end = 60;
    let sub_trace = BackboneTrace::new(trace_orig.coordinates()[start..=end].to_vec());
    let sub_inv_from_trace = extract_curve_invariants(&sub_trace).expect("Sub-trace invariants");
    let sub_inv_from_slice = inv_orig.subcurve(start, end).expect("Subcurve from inv");

    assert_eq!(
        sub_inv_from_slice.segment_lengths.len(),
        sub_inv_from_trace.segment_lengths.len()
    );
    assert_eq!(
        sub_inv_from_slice.curvatures.len(),
        sub_inv_from_trace.curvatures.len()
    );
    assert_eq!(
        sub_inv_from_slice.torsions.len(),
        sub_inv_from_trace.torsions.len()
    );

    for i in 0..sub_inv_from_trace.curvatures.len() {
        let diff = (sub_inv_from_trace.curvatures[i] - sub_inv_from_slice.curvatures[i]).abs();
        assert!(diff < 1e-12, "Subcurve curvature mismatch at {i}: {diff}");
    }

    for i in 0..sub_inv_from_trace.torsions.len() {
        let diff = (sub_inv_from_trace.torsions[i] - sub_inv_from_slice.torsions[i]).abs();
        assert!(diff < 1e-12, "Subcurve torsion mismatch at {i}: {diff}");
    }
}

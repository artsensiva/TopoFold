//! M1 Writhe Normalization and Sign Convention Regression Tests
//!
//! This test suite validates the production writhe implementation against independent
//! mathematical references. Tests are designed to fail before the M1 fix and pass after.
//!
//! Semantic convention: `segment_pair_writhe()` returns the complete mutual contribution
//! of one unordered non-adjacent segment pair to self-writhe. Enumeration via i < j covers
//! each unordered pair exactly once.

use nalgebra::{Point3, Vector3};
use std::f64::consts::PI;
use topofold_core::{compute_local_writhe, compute_total_writhe, BackboneTrace};

// ============================================================================
// INDEPENDENT REFERENCE IMPLEMENTATIONS (TEST-ONLY)
// ============================================================================

/// Independent numerical Gauss integral using midpoint quadrature.
///
/// Computes Wr(C) = (1/4π) ∫∫ ((r(s)-r(t))·(dr/ds × dr/dt)) / |r(s)-r(t)|³ ds dt
/// over the full ordered domain, avoiding adjacent and self segments.
///
/// Returns (writhe, convergence_error_estimate) where convergence is estimated
/// from the difference between m/2 and m resolution.
fn numerical_gauss_writhe(coords: &[Point3<f64>], m: usize) -> (f64, f64) {
    let n = coords.len();
    let mut total = 0.0;
    let t_vals: Vec<f64> = (0..m).map(|k| ((k as f64) + 0.5) / (m as f64)).collect();

    for i in 0..n - 1 {
        let a0 = coords[i];
        let a1 = coords[i + 1];
        let da = a1 - a0;

        for j in 0..n - 1 {
            if (i as i32 - j as i32).abs() <= 1 {
                continue; // Skip adjacent and self
            }

            let b0 = coords[j];
            let b1 = coords[j + 1];
            let db = b1 - b0;
            let cross_product = da.cross(&db);

            for &s in &t_vals {
                let p = a0 + da * s;
                for &t in &t_vals {
                    let q = b0 + db * t;
                    let diff = p - q;
                    let r_cubed = diff.norm().powi(3);
                    if r_cubed > 1e-12 {
                        total += diff.dot(&cross_product) / r_cubed;
                    }
                }
            }
        }
    }

    let wr = total / (4.0 * PI * (m * m) as f64);

    // Convergence estimate: run at m/2 and compare
    let m_half = m / 2;
    if m_half >= 2 {
        let (wr_half, _) = numerical_gauss_writhe_internal(coords, m_half);
        let error = (wr - wr_half).abs();
        (wr, error)
    } else {
        (wr, f64::NAN)
    }
}

fn numerical_gauss_writhe_internal(coords: &[Point3<f64>], m: usize) -> (f64, f64) {
    let n = coords.len();
    let mut total = 0.0;
    let t_vals: Vec<f64> = (0..m).map(|k| ((k as f64) + 0.5) / (m as f64)).collect();

    for i in 0..n - 1 {
        let a0 = coords[i];
        let a1 = coords[i + 1];
        let da = a1 - a0;

        for j in 0..n - 1 {
            if (i as i32 - j as i32).abs() <= 1 {
                continue;
            }

            let b0 = coords[j];
            let b1 = coords[j + 1];
            let db = b1 - b0;
            let cross_product = da.cross(&db);

            for &s in &t_vals {
                let p = a0 + da * s;
                for &t in &t_vals {
                    let q = b0 + db * t;
                    let diff = p - q;
                    let r_cubed = diff.norm().powi(3);
                    if r_cubed > 1e-12 {
                        total += diff.dot(&cross_product) / r_cubed;
                    }
                }
            }
        }
    }

    (total / (4.0 * PI * (m * m) as f64), 0.0)
}

/// Independent polygonal writhe reference using Klenin-Langowski spherical quadrilateral method.
///
/// This implementation is mathematically independent of the production segment_pair_writhe().
/// It computes the signed solid angle using the spherical quadrilateral normal method.
///
/// Returns the complete writhe by summing over all unordered segment pairs i < j.
fn klenin_langowski_writhe(coords: &[Point3<f64>]) -> f64 {
    let n = coords.len();
    let mut total = 0.0;

    for i in 0..n - 2 {
        let r1 = coords[i];
        let r2 = coords[i + 1];

        for j in (i + 2)..n - 1 {
            let r3 = coords[j];
            let r4 = coords[j + 1];

            let r13 = r3 - r1;
            let r14 = r4 - r1;
            let r23 = r3 - r2;
            let r24 = r4 - r2;

            // Four normals of the spherical quadrilateral
            let n1 = r13.cross(&r14);
            let n2 = r14.cross(&r24);
            let n3 = r24.cross(&r23);
            let n4 = r23.cross(&r13);

            let normals = [n1, n2, n3, n4];
            let unit_normals: Vec<Vector3<f64>> = normals
                .iter()
                .map(|n| {
                    let norm = n.norm();
                    if norm > 1e-12 {
                        n / norm
                    } else {
                        Vector3::zeros()
                    }
                })
                .collect();

            // Sum of spherical angles
            let mut omega = 0.0;
            for k in 0..4 {
                let k_next = (k + 1) % 4;
                let dot = unit_normals[k].dot(&unit_normals[k_next]);
                let clamped = dot.clamp(-1.0, 1.0);
                omega += clamped.asin();
            }

            // Sign from segment orientation
            let cross_seg = (r4 - r3).cross(&(r2 - r1));
            let sign = cross_seg.dot(&r13).signum();

            // Klenin-Langowski formula gives contribution per unordered pair
            // Multiply by 2 because we enumerate each unordered pair once with i < j
            total += 2.0 * omega * sign / (4.0 * PI);
        }
    }

    total
}

// ============================================================================
// TEST FIXTURES
// ============================================================================

fn ideal_helix(n: usize, right_handed: bool) -> Vec<Point3<f64>> {
    let r = 2.3;
    let dz = 1.5;
    let theta_deg = if right_handed { 100.0 } else { -100.0 };
    let theta_rad = theta_deg * PI / 180.0;

    (0..n)
        .map(|i| {
            let t = (i as f64) * theta_rad;
            let z = (i as f64) * dz;
            Point3::new(r * t.cos(), r * t.sin(), z)
        })
        .collect()
}

fn random_walk_3d(n: usize, seed: u64) -> Vec<Point3<f64>> {
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    let mut rng = StdRng::seed_from_u64(seed);
    let mut coords = Vec::with_capacity(n);
    let mut pos = Point3::origin();

    for _ in 0..n {
        coords.push(pos);
        let step = Vector3::new(
            rng.gen_range(-2.0..2.0),
            rng.gen_range(-2.0..2.0),
            rng.gen_range(-2.0..2.0),
        );
        pos += step;
    }

    coords
}

fn straight_chain(n: usize) -> Vec<Point3<f64>> {
    (0..n)
        .map(|i| Point3::new(0.0, 0.0, i as f64 * 3.8))
        .collect()
}

fn planar_zigzag(n: usize) -> Vec<Point3<f64>> {
    (0..n)
        .map(|i| {
            let x = (i as f64) * 1.5;
            let y = if i % 2 == 0 { 1.0 } else { -1.0 };
            Point3::new(x, y, 0.0)
        })
        .collect()
}

// ============================================================================
// TEST 1: PRODUCTION VS NUMERICAL GAUSS REFERENCE
// ============================================================================

#[test]
fn test_total_writhe_vs_numerical_gauss() {
    let helix = ideal_helix(25, true);
    let trace = BackboneTrace::new(helix.clone());
    let wr_prod = compute_total_writhe(&trace).expect("compute_total_writhe failed");

    // Test convergence at multiple resolutions
    let resolutions = [20, 40, 80, 160];
    let mut wr_numerical = Vec::new();
    let mut errors = Vec::new();

    println!("\n=== Test 1: Production vs Numerical Gauss ===");
    println!("Resolution (m) | Wr(numerical) | Convergence Error");
    println!("---------------|---------------|------------------");

    for &m in &resolutions {
        let (wr, err) = numerical_gauss_writhe(&helix, m);
        println!("{:14} | {:+13.6} | {:17.6e}", m, wr, err);
        wr_numerical.push(wr);
        errors.push(err);
    }

    // Use highest resolution as reference
    let wr_ref = wr_numerical[wr_numerical.len() - 1];
    let tol = errors[errors.len() - 1] * 2.0; // Tolerance based on convergence error

    println!("\nProduction writhe: {:+.6}", wr_prod);
    println!("Reference writhe (m=160): {:+.6}", wr_ref);
    println!("Tolerance (2×convergence error): {:.6e}", tol);
    println!("Difference: {:+.6e}", (wr_prod - wr_ref).abs());

    assert!(
        (wr_prod - wr_ref).abs() < tol,
        "Production writhe {:.6} disagrees with numerical reference {:.6} (tolerance {:.6e})",
        wr_prod,
        wr_ref,
        tol
    );
}

// ============================================================================
// TEST 2: PRODUCTION VS INDEPENDENT POLYGONAL REFERENCE
// ============================================================================

#[test]
fn test_total_writhe_vs_klenin_langowski() {
    println!("\n=== Test 2: Production vs Klenin-Langowski Polygonal Reference ===");
    println!("Curve                          | Production | K-L Reference | Difference");
    println!("-------------------------------|------------|---------------|------------");

    let test_cases = vec![
        ("Right-handed helix (25)", ideal_helix(25, true)),
        ("Left-handed helix (25)", ideal_helix(25, false)),
        ("Random walk (30, seed=0)", random_walk_3d(30, 0)),
        ("Random walk (30, seed=1)", random_walk_3d(30, 1)),
        ("Random walk (30, seed=42)", random_walk_3d(30, 42)),
    ];

    for (name, coords) in test_cases {
        let trace = BackboneTrace::new(coords.clone());
        let wr_prod = compute_total_writhe(&trace).expect("compute_total_writhe failed");
        let wr_kl = klenin_langowski_writhe(&coords);
        let diff = (wr_prod - wr_kl).abs();

        println!(
            "{:<30} | {:+10.4} | {:+13.4} | {:+10.4e}",
            name, wr_prod, wr_kl, diff
        );

        assert!(
            diff < 1e-7,
            "Production writhe {:.6} disagrees with K-L reference {:.6} for {}",
            wr_prod,
            wr_kl,
            name
        );
    }
}

// ============================================================================
// TEST 3: STRAIGHT CHAIN (ZERO WRITHE)
// ============================================================================

#[test]
fn test_straight_chain_zero_writhe() {
    println!("\n=== Test 3: Straight Chain (Zero Writhe) ===");

    let straight = straight_chain(20);
    let trace = BackboneTrace::new(straight);
    let wr = compute_total_writhe(&trace).expect("compute_total_writhe failed");

    println!("Writhe of straight chain: {:+.6e}", wr);
    println!("Tolerance: 1e-10");

    assert!(
        wr.abs() < 1e-10,
        "Straight chain must have zero writhe, got {:.6e}",
        wr
    );
}

// ============================================================================
// TEST 4: PLANAR CURVE (ZERO WRITHE)
// ============================================================================

#[test]
fn test_planar_curve_zero_writhe() {
    println!("\n=== Test 4: Planar Curve (Zero Writhe) ===");

    let planar = planar_zigzag(20);
    let trace = BackboneTrace::new(planar);
    let wr = compute_total_writhe(&trace).expect("compute_total_writhe failed");

    println!("Writhe of planar zigzag: {:+.6e}", wr);
    println!("Tolerance: 1e-10");

    assert!(
        wr.abs() < 1e-10,
        "Planar curve must have zero writhe, got {:.6e}",
        wr
    );
}

// ============================================================================
// TEST 5: SE(3) INVARIANCE
// ============================================================================

#[test]
fn test_se3_invariance() {
    use nalgebra::{Rotation3, Translation3};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    println!("\n=== Test 5: SE(3) Invariance ===");

    let helix = ideal_helix(25, true);
    let trace = BackboneTrace::new(helix.clone());
    let wr_original = compute_total_writhe(&trace).expect("compute_total_writhe failed");

    let mut rng = StdRng::seed_from_u64(42);
    let mut max_diff: f64 = 0.0;

    for _ in 0..20 {
        // Random rotation
        let axis = Vector3::new(rng.gen(), rng.gen(), rng.gen()).normalize();
        let angle = rng.gen_range(-PI..PI);
        let rotation = Rotation3::from_axis_angle(&nalgebra::Unit::new_normalize(axis), angle);

        // Random translation
        let translation = Translation3::new(
            rng.gen_range(-1000.0..1000.0),
            rng.gen_range(-1000.0..1000.0),
            rng.gen_range(-1000.0..1000.0),
        );

        // Apply transformation
        let transformed: Vec<Point3<f64>> =
            helix.iter().map(|p| translation * (rotation * p)).collect();

        let trace_trans = BackboneTrace::new(transformed);
        let wr_trans = compute_total_writhe(&trace_trans).expect("compute_total_writhe failed");

        let diff = (wr_original - wr_trans).abs();
        max_diff = max_diff.max(diff);
    }

    println!("Original writhe: {:+.6}", wr_original);
    println!(
        "Max difference over 20 random SE(3) transforms: {:.2e}",
        max_diff
    );
    println!("Tolerance: 1e-12");

    assert!(
        max_diff < 1e-12,
        "SE(3) invariance violated: max difference {:.2e}",
        max_diff
    );
}

// ============================================================================
// TEST 6: SPATIAL REFLECTION
// ============================================================================

#[test]
fn test_spatial_reflection() {
    println!("\n=== Test 6: Spatial Reflection ===");

    let helix = ideal_helix(25, true);
    let trace = BackboneTrace::new(helix.clone());
    let wr_original = compute_total_writhe(&trace).expect("compute_total_writhe failed");

    // Reflection via negation (determinant = -1)
    let reflected: Vec<Point3<f64>> = helix
        .iter()
        .map(|p| Point3::new(-p.x, -p.y, -p.z))
        .collect();
    let trace_refl = BackboneTrace::new(reflected);
    let wr_reflected = compute_total_writhe(&trace_refl).expect("compute_total_writhe failed");

    let sum = wr_original + wr_reflected;

    println!("Wr(C) = {:+.6}", wr_original);
    println!("Wr(-C) = {:+.6}", wr_reflected);
    println!("Sum (should be ≈0) = {:+.6e}", sum);
    println!("Tolerance: 1e-12");

    assert!(
        sum.abs() < 1e-12,
        "Reflection law violated: Wr(C) + Wr(-C) = {:.6e}, expected ≈0",
        sum
    );
}

// ============================================================================
// TEST 7: WHOLE-CHAIN REVERSAL
// ============================================================================

#[test]
fn test_chain_reversal_invariance() {
    println!("\n=== Test 7: Whole-Chain Reversal (INVARIANCE) ===");

    let helix = ideal_helix(25, true);
    let trace = BackboneTrace::new(helix.clone());
    let wr_original = compute_total_writhe(&trace).expect("compute_total_writhe failed");

    // Reverse the chain
    let mut reversed = helix.clone();
    reversed.reverse();
    let trace_rev = BackboneTrace::new(reversed);
    let wr_reversed = compute_total_writhe(&trace_rev).expect("compute_total_writhe failed");

    let diff = wr_original - wr_reversed;

    println!("Wr(C) = {:+.6}", wr_original);
    println!("Wr(reversed) = {:+.6}", wr_reversed);
    println!("Difference (should be ≈0 for self-writhe) = {:+.6e}", diff);
    println!("Tolerance: 1e-12");

    assert!(
        diff.abs() < 1e-12,
        "Chain reversal invariance violated: Wr(C) - Wr(reversed) = {:.6e}, expected ≈0",
        diff
    );
}

// ============================================================================
// TEST 8: LOCAL WRITHE CONSISTENCY
// ============================================================================

#[test]
fn test_local_writhe_consistency() {
    println!("\n=== Test 8: Local Writhe Consistency ===");

    let helix = ideal_helix(50, true);
    let trace = BackboneTrace::new(helix.clone());

    let window_radius = 10;
    let local_spectrum =
        compute_local_writhe(&trace, window_radius).expect("compute_local_writhe failed");

    // Pick a specific window (e.g., center at index 25)
    let center_idx = 25;
    let start = center_idx - window_radius;
    let end = center_idx + window_radius + 1; // exclusive

    // Extract the same subchain
    let subchain: Vec<Point3<f64>> = helix[start..end].to_vec();
    let subtrace = BackboneTrace::new(subchain.clone());
    let wr_total_on_sub = compute_total_writhe(&subtrace).expect("compute_total_writhe failed");

    // Local spectrum index
    let local_idx = center_idx - window_radius;
    let wr_local = local_spectrum[local_idx];

    println!("Window: center={}, radius={}", center_idx, window_radius);
    println!("Local writhe from spectrum: {:+.6}", wr_local);
    println!("Total writhe on subchain: {:+.6}", wr_total_on_sub);
    println!("Difference: {:+.6e}", (wr_local - wr_total_on_sub).abs());
    println!("Tolerance: 1e-12");

    assert!(
        (wr_local - wr_total_on_sub).abs() < 1e-12,
        "Local writhe {:.6} disagrees with total writhe on subchain {:.6}",
        wr_local,
        wr_total_on_sub
    );

    // Also compare subchain against numerical reference with convergence-based tolerance
    let (wr_num_40, _) = numerical_gauss_writhe(&subchain, 40);
    let (wr_num_80, _) = numerical_gauss_writhe(&subchain, 80);
    let convergence_error = (wr_num_80 - wr_num_40).abs();
    let tol_local = convergence_error * 2.0;

    println!("\nSubchain vs numerical Gauss:");
    println!("  m=40: {:+.6}", wr_num_40);
    println!("  m=80: {:+.6}", wr_num_80);
    println!("  Convergence error: {:.6e}", convergence_error);
    println!("  Tolerance (2×error): {:.6e}", tol_local);
    println!(
        "  Production vs m=80: {:+.6e}",
        (wr_total_on_sub - wr_num_80).abs()
    );

    assert!(
        (wr_total_on_sub - wr_num_80).abs() < tol_local,
        "Subchain writhe {:.6} disagrees with numerical reference {:.6} (tolerance {:.6e})",
        wr_total_on_sub,
        wr_num_80,
        tol_local
    );
}

//! Integration test for the two-tier conformational trajectory search engine.

use nalgebra::Point3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

use topofold_core::{compute_local_writhe, extract_curve_invariants, BackboneTrace};
use topofold_index::{ConformationalIndex, ConformerFrame};

fn generate_perturbed_helix(n_atoms: usize, noise_scale: f64, rng: &mut StdRng) -> BackboneTrace {
    let r = 2.3;
    let delta_z = 1.5;
    let delta_theta = 100.0 * PI / 180.0;

    let coords: Vec<Point3<f64>> = (0..n_atoms)
        .map(|i| {
            let theta = (i as f64) * delta_theta;
            let z = (i as f64) * delta_z;
            let dx = rng.gen_range(-noise_scale..noise_scale);
            let dy = rng.gen_range(-noise_scale..noise_scale);
            let dz = rng.gen_range(-noise_scale..noise_scale);
            Point3::new(r * theta.cos() + dx, r * theta.sin() + dy, z + dz)
        })
        .collect();

    BackboneTrace::new(coords)
}

#[test]
fn test_ensemble_indexing_and_sub_linear_search() {
    let mut rng = StdRng::seed_from_u64(0x7090_19DE_8000);
    let n_atoms = 30;
    let n_frames = 100;
    let mut frames = Vec::with_capacity(n_frames);

    // Generate an ensemble of 100 frames:
    // Frames 0..80: Thermal fluctuations around a canonical alpha-helix (noise = 0.05 Å)
    // Frames 80..100: Decoy conformers (e.g. high noise or altered pitch, representing folded vs unfolded states)
    for frame_id in 0..n_frames {
        let noise = if frame_id < 80 { 0.05 } else { 0.80 };
        let trace = generate_perturbed_helix(n_atoms, noise, &mut rng);
        let invariants = extract_curve_invariants(&trace).expect("Invariants");
        let writhe_spectrum = compute_local_writhe(&trace, 3).expect("Local writhe");

        frames.push(ConformerFrame {
            frame_id: frame_id as u64,
            time_ps: frame_id as f64 * 10.0,
            invariants,
            writhe_spectrum,
        });
    }

    let target_frame = frames[10].clone();

    // Build the metric index
    let index = ConformationalIndex::build(frames);
    assert_eq!(index.len(), n_frames);

    // Query for conformations structurally similar to Frame 10:
    // Tight tolerances: should retrieve nearby thermal fluctuations, strictly excluding decoys (frames 80..100)
    let hits = index.query_conformations(
        &target_frame.invariants,
        &target_frame.writhe_spectrum,
        0.15, // Fréchet tolerance on (kappa, tau)
        0.05, // Coarse writhe spectrum tolerance
    );

    assert!(!hits.is_empty(), "Query must find at least the exact target");
    assert_eq!(hits[0].frame_id, 10, "Top hit must be Frame 10 itself");
    assert!(hits[0].frechet_distance < 1e-12, "Self distance must be 0");

    // Verify all returned hits are within thermal ensemble (< 80) and decoys are excluded
    for hit in &hits {
        assert!(
            hit.frame_id < 80,
            "Decoy frame {} should have been pruned by metric index",
            hit.frame_id
        );
    }
}

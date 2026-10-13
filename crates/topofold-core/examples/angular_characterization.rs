//! WP-A angular characterization harness: Rust paths R1 and R2 (SYNTHETIC fixtures).
//!
//! Harness code, not production. It calls only public `topofold_core` functions:
//! - R1: `circular_mean`, `mean_resultant_length`, `centered_deviations`, `StreamingMoments`,
//!   `sarles_bimodality_coefficient` on fixture angle arrays (f64).
//! - R2: `compute_bimodality_profile` on traces built from fixture coordinates (f64), plus an
//!   exact replication of its per-window operations via `extract_curve_invariants`. The
//!   replicated BC must equal the profile's BC bit for bit (protocol §3).
//!
//! Usage: cargo run --release -p topofold-core --example angular_characterization --
//!        <fixtures.tsv> <fixture_dir> <out.json> <out_arrays.f64>
//!
//! Floats are written with `{:e}` (shortest round-trip representation); non-finite values are
//! written as `null`.

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::fs;
use std::io::Write as _;

use nalgebra::Point3;
use topofold_core::{
    centered_deviations, circular_mean, compute_bimodality_profile, extract_curve_invariants,
    mean_resultant_length, sarles_bimodality_coefficient, BackboneTrace, StreamingMoments,
};

const PERTURBATIONS: [f64; 4] = [1e-9, -1e-9, 1e-6, -1e-6];

fn num(x: f64) -> String {
    if x.is_finite() {
        format!("{x:e}")
    } else {
        "null".to_string()
    }
}

fn wrap(a: f64) -> f64 {
    let w = (a + PI).rem_euclid(2.0 * PI) - PI;
    if w <= -PI {
        w + 2.0 * PI
    } else {
        w
    }
}

fn read_f64(path: &str) -> Vec<f64> {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    assert!(
        bytes.len().is_multiple_of(8),
        "{path}: size not a multiple of 8"
    );
    bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| f64::from_le_bytes(*c))
        .collect()
}

struct Arrays {
    data: Vec<f64>,
}

impl Arrays {
    /// Appends values and returns their (offset, length) in the arrays file.
    fn push(&mut self, v: &[f64]) -> String {
        let off = self.data.len();
        self.data.extend_from_slice(v);
        format!("[{off},{}]", v.len())
    }
}

fn moments_json(m: &StreamingMoments, sarles: f64) -> String {
    let bc = m.bimodality_coefficient();
    format!(
        "{{\"n\":{},\"variance\":{},\"skewness\":{},\"excess_kurtosis\":{},\"bc\":{},\
         \"sarles_fn_bc\":{},\"sarles_fn_bits_equal\":{}}}",
        m.count(),
        num(m.variance()),
        num(m.skewness()),
        num(m.excess_kurtosis()),
        num(bc),
        num(sarles),
        sarles.to_bits() == bc.to_bits()
    )
}

fn r1_angular(rows: &[Vec<f64>], perturb: bool, arrays: &mut Arrays) -> String {
    let means: Vec<f64> = rows.iter().map(|r| circular_mean(r)).collect();
    let rr: Vec<f64> = rows.iter().map(|r| mean_resultant_length(r)).collect();
    let across_mean = circular_mean(&means);
    let across_r = mean_resultant_length(&means);
    let devs = centered_deviations(&means);
    let mut m = StreamingMoments::new();
    for &d in &devs {
        m.update(d);
    }
    let sarles = sarles_bimodality_coefficient(&devs);
    let mut s = format!(
        "{{\"across_mean\":{},\"across_R\":{},\"moments\":{},\"means\":{},\"R\":{},\"deviations\":{}",
        num(across_mean),
        num(across_r),
        moments_json(&m, sarles),
        arrays.push(&means),
        arrays.push(&rr),
        arrays.push(&devs)
    );
    if perturb {
        let mut parts = Vec::new();
        for d in PERTURBATIONS {
            let mut row = rows[0].clone();
            row[1] += d;
            parts.push(format!(
                "{{\"delta\":{},\"frame0_mean\":{}}}",
                num(d),
                num(circular_mean(&row))
            ));
        }
        let _ = write!(s, ",\"perturbation_frame0\":[{}]", parts.join(","));
    }
    s.push('}');
    s
}

fn r2(rows: &[Vec<f64>], coords: &[f64], w: usize, kappa: f64, arrays: &mut Arrays) -> String {
    let f = rows.len();
    assert_eq!(coords.len(), f * w * 3, "coordinate array size");
    let traces: Vec<BackboneTrace> = (0..f)
        .map(|i| {
            let pts = (0..w)
                .map(|j| {
                    let o = (i * w + j) * 3;
                    Point3::new(coords[o], coords[o + 1], coords[o + 2])
                })
                .collect();
            BackboneTrace::new(pts)
        })
        .collect();

    // Production call.
    let profile = compute_bimodality_profile(&traces, w).expect("compute_bimodality_profile");
    assert_eq!(profile.len(), 1, "fixtures are single-window traces");
    let p = &profile[0];

    // Exact replication of the per-window operations (window i = 0), same order as production.
    let n_kappas = w - 2;
    let n_torsions = w - 3;
    let mut tau_means = Vec::with_capacity(f);
    let mut kappa_means = Vec::with_capacity(f);
    let mut frame_r = Vec::with_capacity(f);
    let mut max_tau_err = 0.0_f64;
    let mut max_kappa_err = 0.0_f64;
    for (i, t) in traces.iter().enumerate() {
        let inv = extract_curve_invariants(t).expect("extract_curve_invariants");
        let k_sum: f64 = inv.curvatures[0..n_kappas].iter().sum();
        kappa_means.push(k_sum / (n_kappas as f64));
        let t_slice = &inv.torsions[0..n_torsions];
        tau_means.push(circular_mean(t_slice));
        frame_r.push(mean_resultant_length(t_slice)); // harness diagnostic, not used by production
        for (j, &tau) in t_slice.iter().enumerate() {
            max_tau_err = max_tau_err.max(wrap(tau - rows[i][j]).abs());
        }
        for &k in &inv.curvatures {
            max_kappa_err = max_kappa_err.max((k - kappa).abs());
        }
    }
    let tau_devs = centered_deviations(&tau_means);
    let mut m_tau = StreamingMoments::new();
    let mut m_kappa = StreamingMoments::new();
    for &k in &kappa_means {
        m_kappa.update(k);
    }
    for &d in &tau_devs {
        m_tau.update(d);
    }
    let rep_tau = m_tau.bimodality_coefficient();
    let rep_kappa = m_kappa.bimodality_coefficient();
    format!(
        "{{\"profile\":{{\"bc_tau\":{},\"bc_kappa\":{},\"bc_theta\":{},\"score\":{}}},\
         \"replication\":{{\"bc_tau\":{},\"bc_kappa\":{},\"bc_tau_bits_equal\":{},\
         \"bc_kappa_bits_equal\":{}}},\"recovery\":{{\"max_tau_circ_err\":{},\"max_kappa_err\":{}}},\
         \"across_mean\":{},\"across_R\":{},\"tau_moments\":{},\"kappa_variance\":{},\
         \"means\":{},\"R\":{},\"deviations\":{}}}",
        num(p.bc_tau),
        num(p.bc_kappa),
        num(p.bc_theta),
        num(p.score),
        num(rep_tau),
        num(rep_kappa),
        rep_tau.to_bits() == p.bc_tau.to_bits(),
        rep_kappa.to_bits() == p.bc_kappa.to_bits(),
        num(max_tau_err),
        num(max_kappa_err),
        num(circular_mean(&tau_means)),
        num(mean_resultant_length(&tau_means)),
        moments_json(&m_tau, sarles_bimodality_coefficient(&tau_devs)),
        num(m_kappa.variance()),
        arrays.push(&tau_means),
        arrays.push(&frame_r),
        arrays.push(&tau_devs)
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        5,
        "usage: angular_characterization <fixtures.tsv> <fixture_dir> <out.json> <out_arrays.f64>"
    );
    let (tsv, dir, out_json, out_arrays) = (&args[1], &args[2], &args[3], &args[4]);
    let manifest = fs::read_to_string(tsv).expect("read manifest");
    let mut arrays = Arrays { data: Vec::new() };
    let mut cases = Vec::new();

    for line in manifest.lines().filter(|l| !l.trim().is_empty()) {
        let c: Vec<&str> = line.split('\t').collect();
        assert_eq!(c.len(), 9, "manifest line: {line}");
        let (id, kind) = (c[0], c[1]);
        let f: usize = c[2].parse().expect("F");
        let ncols: usize = c[3].parse().expect("ncols");
        let w: usize = c[4].parse().expect("W");
        let kappa: f64 = c[5].parse().expect("kappa");
        let flat = read_f64(&format!("{dir}/{}", c[6]));
        assert_eq!(flat.len(), f * ncols, "{id}: torsion array size");
        let rows: Vec<Vec<f64>> = flat.chunks_exact(ncols).map(<[f64]>::to_vec).collect();

        let body = if kind == "scalar" {
            let vals: Vec<f64> = rows.iter().map(|r| r[0]).collect();
            let mut m = StreamingMoments::new();
            for &v in &vals {
                m.update(v);
            }
            format!(
                "\"r1\":{{\"moments\":{}}}",
                moments_json(&m, sarles_bimodality_coefficient(&vals))
            )
        } else {
            let mut b = format!("\"r1\":{}", r1_angular(&rows, c[8] == "1", &mut arrays));
            if c[7] != "-" {
                let coords = read_f64(&format!("{dir}/{}", c[7]));
                let _ = write!(b, ",\"r2\":{}", r2(&rows, &coords, w, kappa, &mut arrays));
            }
            b
        };
        cases.push(format!("\"{id}\":{{{body}}}"));
    }

    let json = format!(
        "{{\"_meta\":{{\"type\":\"SYNTHETIC\",\"producer\":\"crates/topofold-core/examples/angular_characterization.rs\",\
         \"arrays_file\":\"{out_arrays}\",\"array_refs\":\"[offset,length] into arrays_file (little-endian f64)\",\
         \"n_cases\":{}}},\"cases\":{{{}}}}}\n",
        cases.len(),
        cases.join(",")
    );
    fs::write(out_json, json).expect("write json");
    let mut fh = fs::File::create(out_arrays).expect("create arrays");
    for v in &arrays.data {
        fh.write_all(&v.to_le_bytes()).expect("write arrays");
    }
    println!("{} cases -> {out_json}", cases.len());
}

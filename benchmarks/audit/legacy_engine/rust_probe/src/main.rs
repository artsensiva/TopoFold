//! Legacy-engine edge-case probe (research harness, SYNTHETIC inputs). Prints one JSON object.
use topofold_core::distance::{discrete_frechet_invariants, dtw_invariants, writhe_spectrum_distance};
use topofold_core::invariants::CurveInvariants;

fn inv(k: &[f64], t: &[f64]) -> CurveInvariants {
    CurveInvariants {
        segment_lengths: vec![3.8; k.len() + 1],
        curvatures: k.to_vec(),
        torsions: t.to_vec(),
        sidechain_dihedrals: vec![],
    }
}

fn f(a: &CurveInvariants, b: &CurveInvariants) -> f64 {
    discrete_frechet_invariants(a, b, 1.0, 1.0)
}

fn js(x: f64) -> String {
    if x.is_finite() { format!("{x:.17e}") } else { format!("\"{x}\"") }
}

fn main() {
    let k = [1.0, 1.1, 1.2, 1.3, 1.4];
    let t = [0.5, -0.4, 2.9, -3.0];
    let base = inv(&k, &t);
    let mut out: Vec<(String, String)> = Vec::new();

    // 1. empty input vs non-empty
    out.push(("empty_vs_nonempty".into(), js(f(&inv(&[], &[]), &base))));
    // 2. last kappa (vertex W-2) changed by 1 rad: legacy pairing truncates to min(len) = 4 pairs
    let mut k2 = k;
    k2[4] += 1.0;
    out.push(("last_kappa_changed_by_1rad".into(), js(f(&base, &inv(&k2, &t)))));
    // 3. first kappa changed by 1 rad (control)
    let mut k3 = k;
    k3[0] += 1.0;
    out.push(("first_kappa_changed_by_1rad".into(), js(f(&base, &inv(&k3, &t)))));
    // 4. NaN in one torsion of b; compare with the same pair set to a large difference
    let mut t_nan = t;
    t_nan[1] = f64::NAN;
    out.push(("nan_torsion".into(), js(f(&base, &inv(&k, &t_nan)))));
    let mut t_far = t;
    t_far[1] = t[1] + 3.0;
    out.push(("far_torsion_same_position".into(), js(f(&base, &inv(&k, &t_far)))));
    // 5. all-NaN torsions
    out.push(("all_nan_torsions".into(), js(f(&base, &inv(&k, &[f64::NAN; 4])))));
    // 6. torsion wrap: pi-0.01 vs -pi+0.01 at one position
    let mut ta = t;
    let mut tb = t;
    ta[2] = std::f64::consts::PI - 0.01;
    tb[2] = -std::f64::consts::PI + 0.01;
    out.push(("wrap_pair_expected_0p02".into(), js(f(&inv(&k, &ta), &inv(&k, &tb)))));
    // 7. misregistration: sequences [x,y,y,y] vs [x,x,y,y] (equal kappa)
    let kk = [1.0; 5];
    out.push((
        "misregistration_x_y_vs_x_x".into(),
        js(f(&inv(&kk, &[0.0, 2.0, 2.0, 2.0]), &inv(&kk, &[0.0, 0.0, 2.0, 2.0]))),
    ));
    // 8. writhe spectrum distance: unequal lengths truncate; empty -> 0
    out.push(("writhe_unequal_len".into(), js(writhe_spectrum_distance(&[1.0, 2.0], &[1.0, 2.0, 50.0]))));
    out.push(("writhe_empty".into(), js(writhe_spectrum_distance(&[], &[1.0]))));
    // 9. DTW empty
    out.push(("dtw_empty".into(), js(dtw_invariants(&inv(&[], &[]), &base, 1.0, 1.0))));
    // 10. curvature-only array longer than torsion array by more than one (unequal lengths)
    out.push(("kappa_len5_tau_len2".into(), js(f(&inv(&k, &t[..2]), &inv(&k2, &t[..2])))));

    let body: Vec<String> = out.iter().map(|(a, b)| format!("  \"{a}\": {b}")).collect();
    println!("{{\n{}\n}}", body.join(",\n"));
}

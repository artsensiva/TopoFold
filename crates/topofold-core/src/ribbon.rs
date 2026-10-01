//! Side-Chain Orientation via C-beta Ribbon Geometry.
//!
//! Extends 1D polygonal space curves into 2D oriented ribbons by tracking the
//! relative orientation of side-chain C-beta atoms ($\text{C}_\beta$) with respect
//! to the discrete Frenet-Serret frame.
//!
//! # Biophysical Motivation
//!
//! Pure C-alpha traces capture backbone curvature and torsion but can be blind to
//! "rotameric cryptic pocket gating," where bulky aromatic side chains (Trp, Tyr, Phe)
//! swing into or out of a binding pocket without significant C-alpha displacement.
//!
//! By evaluating the side-chain orientation dihedral angle $\theta_\beta$ between the
//! local osculating plane ($\mathbf{T} \times \mathbf{N}$) and the $\text{C}_\alpha \to \text{C}_\beta$
//! unit vector $\mathbf{v}_\beta$, TopoFold constructs an SE(3)-invariant ribbon representation.
//!
//! # Mathematical Foundations
//!
//! 1. **C-beta Unit Vector**:
//!    $$\mathbf{v}_\beta = \frac{\mathbf{r}_{\text{C}\beta} - \mathbf{r}_{\text{C}\alpha}}{\|\mathbf{r}_{\text{C}\beta} - \mathbf{r}_{\text{C}\alpha}\|}$$
//!
//! 2. **Glycine Pseudo-Cbeta**:
//!    Glycine lacks a C-beta atom. A deterministic pseudo-Cbeta vector is constructed from
//!    the backbone peptide bisector:
//!    $$\mathbf{v}_{\text{bisect}} = \frac{(\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_N) + (\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_C)}{\|(\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_N) + (\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_C)\|}$$
//!
//! 3. **Ribbon Dihedral Angle $\theta_\beta$**:
//!    At internal vertex $i$ with Frenet frame $(\mathbf{T}_{\text{vertex}, i}, \mathbf{N}_i, \mathbf{B}_i)$:
//!    - $\mathbf{B}_i = \frac{\mathbf{T}_{i-1} \times \mathbf{T}_i}{\|\mathbf{T}_{i-1} \times \mathbf{T}_i\|}$ (unit normal to the osculating plane)
//!    - $\mathbf{T}_{\text{vertex}, i} = \frac{\mathbf{T}_{i-1} + \mathbf{T}_i}{\|\mathbf{T}_{i-1} + \mathbf{T}_i\|}$ (tangent along the curve)
//!    - $\mathbf{N}_i = \mathbf{B}_i \times \mathbf{T}_{\text{vertex}, i}$ (principal normal in the osculating plane)
//!
//!    The ribbon orientation angle $\theta_{\beta, i}$ is the signed dihedral angle around $\mathbf{T}_{\text{vertex}, i}$:
//!    $$\theta_{\beta, i} = \operatorname{atan2}(\mathbf{v}_{\beta, i} \cdot \mathbf{B}_i, \; \mathbf{v}_{\beta, i} \cdot \mathbf{N}_i) \in (-\pi, \pi]$$
//!
//! 4. **Strict SE(3) Invariance**:
//!    Under any rigid rotation $\mathbf{R} \in SO(3)$ and translation $\mathbf{t} \in \mathbb{R}^3$,
//!    all dot products are preserved:
//!    $$(\mathbf{R} \mathbf{v}_\beta) \cdot (\mathbf{R} \mathbf{B}_i) = \mathbf{v}_\beta \cdot \mathbf{B}_i, \quad (\mathbf{R} \mathbf{v}_\beta) \cdot (\mathbf{R} \mathbf{N}_i) = \mathbf{v}_\beta \cdot \mathbf{N}_i$$
//!    guaranteeing strict gauge invariance down to $< 10^{-12}$.

#![forbid(unsafe_code)]

use nalgebra::{Point3, UnitVector3, Vector3};

use crate::error::GeometryError;
use crate::geometry::GEOMETRY_EPSILON;
use crate::invariants::CurveInvariants;
use crate::types::BackboneTrace;

/// Standard tetrahedral C-alpha to C-beta bond length (in Ångströms).
pub const STANDARD_CA_CB_BOND_LENGTH: f64 = 1.52;

/// A protein backbone trace augmented with side-chain C-beta (or pseudo-C-beta) coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct RibbonTrace {
    /// 3D coordinates of C-alpha atoms (in Ångströms).
    pub(crate) ca_coords: Vec<Point3<f64>>,
    /// 3D coordinates of C-beta (or pseudo-C-beta) atoms (in Ångströms).
    pub(crate) cb_coords: Vec<Point3<f64>>,
}

impl RibbonTrace {
    /// Creates a new `RibbonTrace` from matching C-alpha and C-beta coordinates.
    ///
    /// # Errors
    /// Returns [`GeometryError::InsufficientPoints`] if fewer than 2 atoms are provided,
    /// or if `ca_coords.len() != cb_coords.len()`.
    pub fn new(
        ca_coords: Vec<Point3<f64>>,
        cb_coords: Vec<Point3<f64>>,
    ) -> Result<Self, GeometryError> {
        if ca_coords.len() != cb_coords.len() {
            return Err(GeometryError::InsufficientPoints {
                required: ca_coords.len(),
                actual: cb_coords.len(),
            });
        }
        if ca_coords.len() < 2 {
            return Err(GeometryError::InsufficientPoints {
                required: 2,
                actual: ca_coords.len(),
            });
        }
        Ok(Self { ca_coords, cb_coords })
    }

    /// Creates a `RibbonTrace` from slices of `[f64; 3]` arrays.
    ///
    /// # Errors
    /// Returns [`GeometryError`] if coordinate lengths do not match or are insufficient.
    pub fn from_arrays(ca: &[[f64; 3]], cb: &[[f64; 3]]) -> Result<Self, GeometryError> {
        let ca_pts: Vec<Point3<f64>> = ca.iter().map(|&[x, y, z]| Point3::new(x, y, z)).collect();
        let cb_pts: Vec<Point3<f64>> = cb.iter().map(|&[x, y, z]| Point3::new(x, y, z)).collect();
        Self::new(ca_pts, cb_pts)
    }

    /// Returns the number of residues in the ribbon trace.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.ca_coords.len()
    }

    /// Returns `true` if the ribbon trace is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ca_coords.is_empty()
    }

    /// Provides slice access to C-alpha coordinates.
    #[inline]
    #[must_use]
    pub fn ca_coords(&self) -> &[Point3<f64>] {
        &self.ca_coords
    }

    /// Provides slice access to C-beta coordinates.
    #[inline]
    #[must_use]
    pub fn cb_coords(&self) -> &[Point3<f64>] {
        &self.cb_coords
    }

    /// Converts the ribbon trace into a standard [`BackboneTrace`] with attached C-beta coordinates.
    #[must_use]
    pub fn to_backbone_trace(&self) -> BackboneTrace {
        BackboneTrace::with_cbeta(self.ca_coords.clone(), self.cb_coords.clone())
    }

    /// Computes the unit vectors pointing from C-alpha to C-beta for each residue.
    ///
    /// # Errors
    /// Returns [`GeometryError::DegenerateSegment`] if any C-alpha to C-beta distance
    /// is less than [`GEOMETRY_EPSILON`].
    pub fn cbeta_unit_vectors(&self) -> Result<Vec<UnitVector3<f64>>, GeometryError> {
        let mut vectors = Vec::with_capacity(self.len());
        for i in 0..self.len() {
            let diff: Vector3<f64> = self.cb_coords[i] - self.ca_coords[i];
            let len = diff.norm();
            if len < GEOMETRY_EPSILON {
                return Err(GeometryError::DegenerateSegment {
                    index: i,
                    next_index: i,
                    length: len,
                });
            }
            vectors.push(UnitVector3::new_unchecked(diff / len));
        }
        Ok(vectors)
    }
}

/// Constructs a deterministic pseudo-Cbeta coordinate for Glycine residues using the
/// backbone peptide bisector:
///
/// $$\mathbf{v}_{\text{bisect}} = \frac{(\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_N) + (\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_C)}{\|(\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_N) + (\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_C)\|}$$
/// $$\mathbf{r}_{\text{pseudo-CB}} = \mathbf{r}_{\text{C}\alpha} + 1.52 \cdot \mathbf{v}_{\text{bisect}}$$
///
/// If the bisector vector has norm below [`GEOMETRY_EPSILON`], a deterministic perpendicular
/// vector is chosen to guarantee smooth, non-NaN numerical behavior.
#[must_use]
pub fn compute_glycine_pseudo_cbeta(
    r_ca: Point3<f64>,
    r_n: Point3<f64>,
    r_c: Point3<f64>,
) -> Point3<f64> {
    let d_n: Vector3<f64> = r_ca - r_n;
    let d_c: Vector3<f64> = r_ca - r_c;
    let sum = d_n + d_c;
    let norm = sum.norm();

    let v_bisect = if norm >= GEOMETRY_EPSILON {
        sum / norm
    } else {
        // Fallback: cross product of d_n and d_c
        let cross = d_n.cross(&d_c);
        let cross_norm = cross.norm();
        if cross_norm >= GEOMETRY_EPSILON {
            cross / cross_norm
        } else {
            // Perpendicular to d_n
            let perp = if d_n.x.abs() > 0.1 || d_n.y.abs() > 0.1 {
                Vector3::new(-d_n.y, d_n.x, 0.0)
            } else {
                Vector3::new(0.0, -d_n.z, d_n.y)
            };
            let perp_norm = perp.norm();
            if perp_norm >= GEOMETRY_EPSILON {
                perp / perp_norm
            } else {
                Vector3::new(1.0, 0.0, 0.0)
            }
        }
    };

    r_ca + STANDARD_CA_CB_BOND_LENGTH * v_bisect
}

/// Constructs a fallback pseudo-Cbeta coordinate for internal residues when N and C atoms
/// are not available, using the adjacent C-alpha backbone bisector:
///
/// $$\mathbf{d}_1 = \mathbf{r}_{\text{C}\alpha} - \mathbf{r}_{\text{prev}}$$
/// $$\mathbf{d}_2 = \mathbf{r}_{\text{C}\alpha} - \mathbf{r}_{\text{next}}$$
/// $$\mathbf{v}_{\text{bisect}} = \frac{\mathbf{d}_1 + \mathbf{d}_2}{\|\mathbf{d}_1 + \mathbf{d}_2\|}$$
///
/// This provides a smooth, singularity-free ribbon vector even when working with C-alpha-only datasets.
#[must_use]
pub fn compute_pseudo_cbeta_from_ca(
    prev_ca: Point3<f64>,
    curr_ca: Point3<f64>,
    next_ca: Point3<f64>,
) -> Point3<f64> {
    let d1: Vector3<f64> = curr_ca - prev_ca;
    let d2: Vector3<f64> = curr_ca - next_ca;
    let sum = d1 + d2;
    let norm = sum.norm();

    let v_bisect = if norm >= GEOMETRY_EPSILON {
        sum / norm
    } else {
        let diff: Vector3<f64> = next_ca - prev_ca;
        let perp = if diff.x.abs() > 0.1 || diff.y.abs() > 0.1 {
            Vector3::new(-diff.y, diff.x, 0.0)
        } else {
            Vector3::new(0.0, -diff.z, diff.y)
        };
        let perp_norm = perp.norm();
        if perp_norm >= GEOMETRY_EPSILON {
            perp / perp_norm
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        }
    };

    curr_ca + STANDARD_CA_CB_BOND_LENGTH * v_bisect
}

/// Computes the side-chain orientation dihedral angles $\theta_\beta$ for all internal residues.
///
/// For internal vertex $i \in \{1, \dots, N-2\}$:
/// - Evaluates discrete Frenet-Serret framing $(\mathbf{T}_{\text{vertex}, i}, \mathbf{N}_i, \mathbf{B}_i)$
/// - Evaluates C-beta unit direction $\mathbf{v}_{\beta, i} = \frac{\mathbf{r}_{\text{C}\beta, i} - \mathbf{r}_{\text{C}\alpha, i}}{\|\mathbf{r}_{\text{C}\beta, i} - \mathbf{r}_{\text{C}\alpha, i}\|}$
/// - Returns $\theta_{\beta, i} = \operatorname{atan2}(\mathbf{v}_{\beta, i} \cdot \mathbf{B}_i, \; \mathbf{v}_{\beta, i} \cdot \mathbf{N}_i) \in (-\pi, \pi]$
///
/// Output length: `N - 2` (matching discrete curvature count).
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `ca_coords.len() < 3` or if lengths do not match.
/// Returns [`GeometryError::DegenerateSegment`] if any bond length is $< \text{GEOMETRY\_EPSILON}$.
/// Returns [`GeometryError::CollinearSegments`] if any consecutive tangents are collinear.
pub fn compute_sidechain_dihedrals(
    ca_coords: &[Point3<f64>],
    cb_coords: &[Point3<f64>],
) -> Result<Vec<f64>, GeometryError> {
    let n = ca_coords.len();
    if n != cb_coords.len() {
        return Err(GeometryError::InsufficientPoints {
            required: n,
            actual: cb_coords.len(),
        });
    }
    if n < 3 {
        return Err(GeometryError::InsufficientPoints {
            required: 3,
            actual: n,
        });
    }

    // 1. Compute unit tangent vectors along C-alpha segments (length N - 1)
    let mut tangents = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let diff: Vector3<f64> = ca_coords[i + 1] - ca_coords[i];
        let len = diff.norm();
        if len < GEOMETRY_EPSILON {
            return Err(GeometryError::DegenerateSegment {
                index: i,
                next_index: i + 1,
                length: len,
            });
        }
        tangents.push(UnitVector3::new_unchecked(diff / len));
    }

    // 2. Compute C-beta unit vectors (length N)
    let mut cb_vectors = Vec::with_capacity(n);
    for i in 0..n {
        let diff: Vector3<f64> = cb_coords[i] - ca_coords[i];
        let len = diff.norm();
        if len < GEOMETRY_EPSILON {
            return Err(GeometryError::DegenerateSegment {
                index: i,
                next_index: i,
                length: len,
            });
        }
        cb_vectors.push(UnitVector3::new_unchecked(diff / len));
    }

    // 3. Compute theta_beta for each internal vertex i = 1 .. N - 2
    let mut dihedrals = Vec::with_capacity(n - 2);

    for i in 1..n - 1 {
        let t_prev = &tangents[i - 1];
        let t_curr = &tangents[i];

        // Binormal B_i = normalized(T_{i-1} x T_i)
        let cross = t_prev.cross(t_curr);
        let cross_norm = cross.norm();
        if cross_norm < GEOMETRY_EPSILON {
            let dot = t_prev.dot(t_curr).clamp(-1.0, 1.0);
            return Err(GeometryError::CollinearSegments {
                vertex_index: i,
                turning_angle: dot.acos(),
            });
        }
        let b_i = UnitVector3::new_unchecked(cross / cross_norm);

        // Vertex tangent T_vertex = normalized(T_{i-1} + T_i)
        let t_sum = t_prev.as_ref() + t_curr.as_ref();
        let t_sum_norm = t_sum.norm();
        let t_vertex = if t_sum_norm >= GEOMETRY_EPSILON {
            UnitVector3::new_unchecked(t_sum / t_sum_norm)
        } else {
            *t_curr
        };

        // In-plane principal normal N_i = B_i x T_vertex
        let n_i = UnitVector3::new_unchecked(b_i.cross(&t_vertex));

        // Sidechain C-beta unit vector at vertex i
        let v_beta = &cb_vectors[i];

        let proj_normal = v_beta.dot(&n_i);
        let proj_binormal = v_beta.dot(&b_i);

        let theta = proj_binormal.atan2(proj_normal);
        dihedrals.push(theta);
    }

    Ok(dihedrals)
}

/// Extracts complete SE(3)-invariant curve and ribbon invariants for a [`RibbonTrace`].
///
/// Returns [`CurveInvariants`] populated with:
/// - `segment_lengths`: length $N - 1$
/// - `curvatures`: length $N - 2$
/// - `torsions`: length $N - 3$
/// - `sidechain_dihedrals`: length $N - 2$
///
/// # Errors
/// Returns [`GeometryError`] if trace has fewer than 4 points or degenerate geometry.
pub fn extract_ribbon_invariants(trace: &RibbonTrace) -> Result<CurveInvariants, GeometryError> {
    if trace.len() < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: trace.len(),
        });
    }

    let bb_trace = trace.to_backbone_trace();
    let (segment_lengths, tangents) = crate::geometry::compute_tangents(&bb_trace)?;
    let curvatures = crate::geometry::compute_curvatures(&tangents)?;
    let torsions = crate::geometry::compute_torsions(&tangents)?;
    let sidechain_dihedrals = compute_sidechain_dihedrals(trace.ca_coords(), trace.cb_coords())?;

    Ok(CurveInvariants {
        segment_lengths,
        curvatures,
        torsions,
        sidechain_dihedrals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Isometry3, Translation3, UnitQuaternion};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::f64::consts::PI;

    #[test]
    fn test_se3_strict_invariance_theta_beta() {
        let mut rng = StdRng::seed_from_u64(42);

        // Generate synthetic backbone (10 residues)
        let n = 10;
        let mut ca_coords = Vec::with_capacity(n);
        let mut cb_coords = Vec::with_capacity(n);

        for i in 0..n {
            let t = i as f64 * 0.8;
            let ca = Point3::new(
                2.3 * t.cos(),
                2.3 * t.sin(),
                1.5 * t,
            );
            // Sidechain pointing outward with some vertical tilt
            let cb = ca + Vector3::new(
                1.2 * (t + 0.3).cos(),
                1.2 * (t + 0.3).sin(),
                0.8,
            );
            ca_coords.push(ca);
            cb_coords.push(cb);
        }

        let theta_orig = compute_sidechain_dihedrals(&ca_coords, &cb_coords).unwrap();
        assert_eq!(theta_orig.len(), n - 2);

        // Test across 100 random SO(3) rotations and R^3 translations
        for _ in 0..100 {
            let roll: f64 = rng.gen_range(0.0..2.0 * PI);
            let pitch: f64 = rng.gen_range(0.0..2.0 * PI);
            let yaw: f64 = rng.gen_range(0.0..2.0 * PI);
            let rot = UnitQuaternion::from_euler_angles(roll, pitch, yaw);

            let tx: f64 = rng.gen_range(-100.0..100.0);
            let ty: f64 = rng.gen_range(-100.0..100.0);
            let tz: f64 = rng.gen_range(-100.0..100.0);
            let trans = Translation3::new(tx, ty, tz);

            let isometry = Isometry3::from_parts(trans, rot);

            let ca_rot: Vec<Point3<f64>> = ca_coords.iter().map(|p| isometry.transform_point(p)).collect();
            let cb_rot: Vec<Point3<f64>> = cb_coords.iter().map(|p| isometry.transform_point(p)).collect();

            let theta_rot = compute_sidechain_dihedrals(&ca_rot, &cb_rot).unwrap();

            for (orig, rot_val) in theta_orig.iter().zip(theta_rot.iter()) {
                let diff = (orig - rot_val).abs();
                let diff_s1 = diff.min((2.0 * PI - diff).abs());
                assert!(
                    diff_s1 < 1e-12,
                    "SE(3) invariance violated: orig={orig}, rot={rot_val}, diff={diff_s1}"
                );
            }
        }
    }

    #[test]
    fn test_glycine_pseudo_cbeta_smoothness() {
        let r_ca = Point3::new(0.0, 0.0, 0.0);
        let r_n = Point3::new(-1.0, 1.0, 0.0);
        let r_c = Point3::new(1.0, 1.0, 0.0);

        let cb = compute_glycine_pseudo_cbeta(r_ca, r_n, r_c);
        let diff: Vector3<f64> = cb - r_ca;
        let len = diff.norm();

        assert!(!cb.x.is_nan() && !cb.y.is_nan() && !cb.z.is_nan());
        assert!((len - STANDARD_CA_CB_BOND_LENGTH).abs() < 1e-9);

        // Backbone bisector of (-1, 1, 0) and (1, 1, 0) from (0,0,0) is in -Y direction:
        // (r_ca - r_n) + (r_ca - r_c) = (1, -1, 0) + (-1, -1, 0) = (0, -2, 0) -> direction (0, -1, 0)
        assert!(cb.y < -1.0);
        assert!(cb.x.abs() < 1e-9);
        assert!(cb.z.abs() < 1e-9);
    }

    #[test]
    fn test_pseudo_cbeta_from_ca_smoothness() {
        let prev = Point3::new(-1.0, 1.0, 0.0);
        let curr = Point3::new(0.0, 0.0, 0.0);
        let next = Point3::new(1.0, 1.0, 0.0);

        let cb = compute_pseudo_cbeta_from_ca(prev, curr, next);
        let diff: Vector3<f64> = cb - curr;
        let len = diff.norm();

        assert!(!cb.x.is_nan() && !cb.y.is_nan() && !cb.z.is_nan());
        assert!((len - STANDARD_CA_CB_BOND_LENGTH).abs() < 1e-9);
        assert!(cb.y < -1.0);
    }
}

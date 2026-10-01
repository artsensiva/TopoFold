//! Core data structures representing protein backbone geometry.

use nalgebra::Point3;

/// An ordered sequence of C-alpha coordinates in 3D space representing a protein backbone trace.
///
/// In standard structural biology, consecutive C-alpha atoms in a *trans*-peptide backbone
/// are separated by approximately 3.80 Å.
#[derive(Debug, Clone, PartialEq)]
pub struct BackboneTrace {
    /// 3D coordinates of C-alpha atoms (in Ångströms).
    pub(crate) coordinates: Vec<Point3<f64>>,
    /// Optional 3D coordinates of C-beta (or pseudo-C-beta) atoms (in Ångströms).
    pub(crate) cb_coordinates: Option<Vec<Point3<f64>>>,
}

impl BackboneTrace {
    /// Creates a new `BackboneTrace` from a vector of 3D C-alpha points.
    #[must_use]
    pub fn new(coordinates: Vec<Point3<f64>>) -> Self {
        Self {
            coordinates,
            cb_coordinates: None,
        }
    }

    /// Creates a `BackboneTrace` with matching C-alpha and C-beta coordinates.
    #[must_use]
    pub fn with_cbeta(coordinates: Vec<Point3<f64>>, cb_coordinates: Vec<Point3<f64>>) -> Self {
        Self {
            coordinates,
            cb_coordinates: Some(cb_coordinates),
        }
    }

    /// Creates a `BackboneTrace` from a slice of `[f64; 3]` arrays.
    #[must_use]
    pub fn from_arrays(coords: &[[f64; 3]]) -> Self {
        let points = coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect();
        Self {
            coordinates: points,
            cb_coordinates: None,
        }
    }

    /// Creates a `BackboneTrace` with C-beta coordinates from slices of `[f64; 3]` arrays.
    #[must_use]
    pub fn from_arrays_with_cbeta(ca_coords: &[[f64; 3]], cb_coords: &[[f64; 3]]) -> Self {
        let ca_points = ca_coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect();
        let cb_points = cb_coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect();
        Self {
            coordinates: ca_points,
            cb_coordinates: Some(cb_points),
        }
    }

    /// Returns the number of C-alpha atoms in the backbone trace.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.coordinates.len()
    }

    /// Returns `true` if the backbone trace contains no atoms.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.coordinates.is_empty()
    }

    /// Provides slice access to the Cartesian coordinates.
    #[inline]
    #[must_use]
    pub fn coordinates(&self) -> &[Point3<f64>] {
        &self.coordinates
    }

    /// Returns a coordinate at index `i`, if in bounds.
    #[inline]
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&Point3<f64>> {
        self.coordinates.get(index)
    }

    /// Provides slice access to C-beta coordinates if present.
    #[inline]
    #[must_use]
    pub fn cb_coordinates(&self) -> Option<&[Point3<f64>]> {
        self.cb_coordinates.as_deref()
    }

    /// Returns `true` if the backbone trace contains C-beta coordinates.
    #[inline]
    #[must_use]
    pub fn has_cbeta(&self) -> bool {
        self.cb_coordinates.is_some()
    }

    /// Attaches or updates C-beta coordinates for the backbone trace.
    pub fn set_cb_coordinates(&mut self, cb_coords: Vec<Point3<f64>>) {
        self.cb_coordinates = Some(cb_coords);
    }
}

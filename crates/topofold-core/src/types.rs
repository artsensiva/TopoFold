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
}

impl BackboneTrace {
    /// Creates a new `BackboneTrace` from a vector of 3D points.
    #[must_use]
    pub fn new(coordinates: Vec<Point3<f64>>) -> Self {
        Self { coordinates }
    }

    /// Creates a `BackboneTrace` from a slice of `[f64; 3]` arrays.
    #[must_use]
    pub fn from_arrays(coords: &[[f64; 3]]) -> Self {
        let points = coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect();
        Self { coordinates: points }
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
}

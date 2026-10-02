//! Python bindings (PyO3) for the TopoFold conformational trajectory engine.

use nalgebra::Point3;
use numpy::{
    PyArray1, PyArray2, PyArray3, PyReadonlyArray2, PyReadonlyArray3, PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufReader;

use topofold_core::error::GeometryError;
use topofold_core::invariants::CurveInvariants;
use topofold_core::pdb::{parse_pdb_ca, parse_pdb_rna, PdbError};
use topofold_core::rna::{
    compute_rna_bimodality_profile as core_rna_bimodality_profile,
    extract_rna_ribbon_invariants_with_window,
    scan_rna_switching_hinges as core_scan_rna_switching_hinges, RnaRibbonTrace,
};
use topofold_core::{
    compute_bimodality_profile as core_bimodality_profile, compute_idp_density_profile,
    compute_intermolecular_allosteric_network as core_intermolecular_allosteric_network,
    compute_intrinsic_allosteric_network, compute_local_writhe,
    compute_ternary_cooperativity_index as core_ternary_cooperativity_index,
    detect_bistable_segments, detect_transient_motifs_with_params,
    discrete_frechet_distance_coords, extract_curve_invariants, BackboneTrace,
};
use topofold_index::{ConformationalIndex as RustConformationalIndex, ConformerFrame};
use topofold_io::{DcdReader, TrajectoryError};

/// Helper to convert a GeometryError to a PyValueError.
fn to_py_err(err: GeometryError) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// Helper to convert a PdbError to a PyValueError.
fn to_pdb_err(err: PdbError) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// Helper to convert a TrajectoryError to a PyValueError.
fn to_traj_err(err: TrajectoryError) -> PyErr {
    PyValueError::new_err(err.to_string())
}

fn trace_to_array2(trace: &BackboneTrace) -> ndarray::Array2<f32> {
    let n = trace.len();
    let mut arr = ndarray::Array2::<f32>::zeros((n, 3));
    for (i, pt) in trace.coordinates().iter().enumerate() {
        arr[[i, 0]] = pt.x as f32;
        arr[[i, 1]] = pt.y as f32;
        arr[[i, 2]] = pt.z as f32;
    }
    arr
}

fn traces_to_array3(traces: &[BackboneTrace]) -> Result<ndarray::Array3<f32>, PyErr> {
    if traces.is_empty() {
        return Err(PyValueError::new_err("Trajectory contains 0 frames"));
    }
    let f_count = traces.len();
    let n_residues = traces[0].len();
    let mut arr = ndarray::Array3::<f32>::zeros((f_count, n_residues, 3));
    for (f_idx, trace) in traces.iter().enumerate() {
        if trace.len() != n_residues {
            return Err(PyValueError::new_err(format!(
                "Frame {} has {} residues, expected {}",
                f_idx,
                trace.len(),
                n_residues
            )));
        }
        for (r_idx, pt) in trace.coordinates().iter().enumerate() {
            arr[[f_idx, r_idx, 0]] = pt.x as f32;
            arr[[f_idx, r_idx, 1]] = pt.y as f32;
            arr[[f_idx, r_idx, 2]] = pt.z as f32;
        }
    }
    Ok(arr)
}

fn cb_trace_to_array2(trace: &BackboneTrace) -> Result<ndarray::Array2<f32>, PyErr> {
    let cb_coords = trace
        .cb_coordinates()
        .ok_or_else(|| PyValueError::new_err("BackboneTrace contains no C-beta coordinates"))?;
    let n = cb_coords.len();
    let mut arr = ndarray::Array2::<f32>::zeros((n, 3));
    for (i, pt) in cb_coords.iter().enumerate() {
        arr[[i, 0]] = pt.x as f32;
        arr[[i, 1]] = pt.y as f32;
        arr[[i, 2]] = pt.z as f32;
    }
    Ok(arr)
}

fn cb_traces_to_array3(traces: &[BackboneTrace]) -> Result<ndarray::Array3<f32>, PyErr> {
    if traces.is_empty() {
        return Err(PyValueError::new_err("Trajectory contains 0 frames"));
    }
    let f_count = traces.len();
    let n_residues = traces[0].len();
    let mut arr = ndarray::Array3::<f32>::zeros((f_count, n_residues, 3));
    for (f_idx, trace) in traces.iter().enumerate() {
        let cb_coords = trace.cb_coordinates().ok_or_else(|| {
            PyValueError::new_err(format!("Frame {} has no C-beta coordinates", f_idx))
        })?;
        if cb_coords.len() != n_residues {
            return Err(PyValueError::new_err(format!(
                "Frame {} has {} C-beta atoms, expected {}",
                f_idx,
                cb_coords.len(),
                n_residues
            )));
        }
        for (r_idx, pt) in cb_coords.iter().enumerate() {
            arr[[f_idx, r_idx, 0]] = pt.x as f32;
            arr[[f_idx, r_idx, 1]] = pt.y as f32;
            arr[[f_idx, r_idx, 2]] = pt.z as f32;
        }
    }
    Ok(arr)
}

fn traces_to_index(
    traces: Vec<BackboneTrace>,
    window_size: usize,
    coarse_radius: f64,
) -> Result<ConformationalIndex, PyErr> {
    if traces.is_empty() {
        return Err(PyValueError::new_err("Trajectory contains 0 frames"));
    }
    let n_residues = traces[0].len();
    let window_radius = (window_size / 2).max(1);
    if n_residues < 2 * window_radius + 2 {
        return Err(PyValueError::new_err(format!(
            "Backbone length ({}) is too short for window_size {} (requires at least {} residues)",
            n_residues,
            window_size,
            2 * window_radius + 2
        )));
    }

    let frames: Result<Vec<ConformerFrame>, PyErr> = traces
        .into_par_iter()
        .enumerate()
        .map(|(f_idx, trace)| {
            if trace.len() != n_residues {
                return Err(PyValueError::new_err(format!(
                    "Frame {} has {} residues, expected {}",
                    f_idx,
                    trace.len(),
                    n_residues
                )));
            }
            let invariants = extract_curve_invariants(&trace).map_err(to_py_err)?;
            let writhe_spectrum = compute_local_writhe(&trace, window_radius).map_err(to_py_err)?;
            Ok(ConformerFrame {
                frame_id: f_idx as u64,
                time_ps: f_idx as f64,
                invariants,
                writhe_spectrum,
            })
        })
        .collect();

    let frames = frames?;
    let inner = RustConformationalIndex::build(frames);
    Ok(ConformationalIndex {
        window_size,
        coarse_radius,
        expected_residues: Some(n_residues),
        inner: Some(inner),
    })
}

/// Helper to convert a 2D ndarray view of shape (N, 3) to BackboneTrace.
fn view_to_trace(view: ndarray::ArrayView2<'_, f32>) -> Result<BackboneTrace, PyErr> {
    let n = view.shape()[0];
    if view.shape()[1] != 3 {
        return Err(PyValueError::new_err(format!(
            "Expected second dimension to be 3 for (x, y, z), found {}",
            view.shape()[1]
        )));
    }

    let mut points = Vec::with_capacity(n);
    for row in view.rows() {
        points.push(Point3::new(row[0] as f64, row[1] as f64, row[2] as f64));
    }

    Ok(BackboneTrace::new(points))
}

/// Computes the SE(3)-invariant discrete differential geometry primitives of a C-alpha backbone.
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Cartesian coordinates of C-alpha atoms in Ångströms.
/// Computes discrete differential geometry invariants (kappa, tau, writhe, and optional theta_beta).
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Coordinates of N C-alpha atoms in Ångströms.
/// cb_coords : numpy.ndarray of shape (N, 3), dtype=float32, optional
///     Optional coordinates of N C-beta atoms in Ångströms.
///
/// Returns
/// -------
/// tuple
///     - If cb_coords is None: `(kappa, tau, writhe_spectrum)`
///     - If cb_coords is provided: `(kappa, tau, writhe_spectrum, theta_beta)`
#[pyfunction]
#[pyo3(signature = (coords, cb_coords = None))]
pub fn compute_invariants<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f32>,
    cb_coords: Option<PyReadonlyArray2<'py, f32>>,
) -> PyResult<Py<PyAny>> {
    let shape = coords.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err(format!(
            "Expected 2D array of shape (N, 3), found {:?}",
            shape
        )));
    }

    let n = shape[0];
    if n < 4 {
        return Err(PyValueError::new_err(format!(
            "Protein backbone requires at least 4 C-alpha atoms, found {}",
            n
        )));
    }

    // Adaptive window radius: half-width of writhe window
    let window_radius = if n >= 12 {
        4
    } else if n >= 8 {
        2
    } else {
        1
    };

    let array_view = coords.as_array();
    let opt_cb_pts: Option<Vec<Point3<f64>>> = if let Some(ref cb) = cb_coords {
        let cb_view = cb.as_array();
        if cb_view.shape() != shape {
            return Err(PyValueError::new_err(format!(
                "cb_coords shape {:?} does not match coords shape {:?}",
                cb_view.shape(),
                shape
            )));
        }
        Some(
            cb_view
                .outer_iter()
                .map(|r| Point3::new(r[0] as f64, r[1] as f64, r[2] as f64))
                .collect(),
        )
    } else {
        None
    };

    let (invariants, writhe_spectrum) =
        py.detach(|| -> Result<(CurveInvariants, Vec<f64>), PyErr> {
            let mut trace = view_to_trace(array_view)?;
            if let Some(cb_pts) = opt_cb_pts {
                trace.set_cb_coordinates(cb_pts);
            }
            let inv = extract_curve_invariants(&trace).map_err(to_py_err)?;
            let wr = compute_local_writhe(&trace, window_radius).map_err(to_py_err)?;
            Ok((inv, wr))
        })?;

    let py_kappa = PyArray1::from_vec(py, invariants.curvatures).unbind();
    let py_tau = PyArray1::from_vec(py, invariants.torsions).unbind();
    let py_writhe = PyArray1::from_vec(py, writhe_spectrum).unbind();

    if cb_coords.is_some() {
        let py_theta = PyArray1::from_vec(py, invariants.sidechain_dihedrals).unbind();
        let tuple = (py_kappa, py_tau, py_writhe, py_theta);
        tuple.into_py_any(py)
    } else {
        let tuple = (py_kappa, py_tau, py_writhe);
        tuple.into_py_any(py)
    }
}

/// Computes the side-chain orientation dihedral angles (theta_beta) between the osculating plane and C-beta.
///
/// Parameters
/// ----------
/// ca_coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Coordinates of N C-alpha atoms in Ångströms.
/// cb_coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Coordinates of N C-beta atoms in Ångströms.
///
/// Returns
/// -------
/// numpy.ndarray of shape (N - 2,), dtype=float64
///     Orientation dihedral angles in radians in (-pi, pi].
#[pyfunction]
#[pyo3(signature = (ca_coords, cb_coords))]
pub fn compute_theta_beta<'py>(
    py: Python<'py>,
    ca_coords: PyReadonlyArray2<'py, f32>,
    cb_coords: PyReadonlyArray2<'py, f32>,
) -> PyResult<Py<PyArray1<f64>>> {
    let ca_view = ca_coords.as_array();
    let cb_view = cb_coords.as_array();

    if ca_view.shape() != cb_view.shape() {
        return Err(PyValueError::new_err(format!(
            "ca_coords shape {:?} does not match cb_coords shape {:?}",
            ca_view.shape(),
            cb_view.shape()
        )));
    }
    if ca_view.shape()[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 coordinates"));
    }
    let n = ca_view.shape()[0];
    if n < 3 {
        return Err(PyValueError::new_err(
            "Theta_beta requires at least 3 residues",
        ));
    }

    let ca_pts: Vec<Point3<f64>> = ca_view
        .outer_iter()
        .map(|r| Point3::new(r[0] as f64, r[1] as f64, r[2] as f64))
        .collect();
    let cb_pts: Vec<Point3<f64>> = cb_view
        .outer_iter()
        .map(|r| Point3::new(r[0] as f64, r[1] as f64, r[2] as f64))
        .collect();

    let thetas = py.detach(|| {
        topofold_core::ribbon::compute_sidechain_dihedrals(&ca_pts, &cb_pts).map_err(to_py_err)
    })?;

    Ok(PyArray1::from_vec(py, thetas).unbind())
}

/// Two-tier metric search index for protein conformational trajectories.
///
/// Indexes conformational ensembles using a Vantage-Point Tree over localized
/// Writhe spectra ($L_2$ metric) for sub-linear coarse pruning, followed by fine
/// verification using the exact Discrete Fréchet distance on $(\kappa, \tau)$.
#[pyclass]
pub struct ConformationalIndex {
    window_size: usize,
    coarse_radius: f64,
    expected_residues: Option<usize>,
    inner: Option<RustConformationalIndex>,
}

#[pymethods]
impl ConformationalIndex {
    /// Creates a new, unpopulated ConformationalIndex.
    ///
    /// Parameters
    /// ----------
    /// window_size : int, default=10
    ///     Window length in residues for the localized Writhe spectrum.
    /// coarse_radius : float, default=0.5
    ///     Threshold radius for Tier 1 VP-Tree candidate filtering.
    #[new]
    #[pyo3(signature = (window_size = 10, coarse_radius = 0.5))]
    pub fn new(window_size: usize, coarse_radius: f64) -> Self {
        Self {
            window_size,
            coarse_radius,
            expected_residues: None,
            inner: None,
        }
    }

    /// Fits the index on an ensemble of conformational frames.
    ///
    /// Parameters
    /// ----------
    /// coords : numpy.ndarray of shape (F, N, 3), dtype=float32
    ///     Trajectory coordinates of F frames, each containing N residues.
    pub fn fit(&mut self, py: Python<'_>, coords: PyReadonlyArray3<f32>) -> PyResult<()> {
        let shape = coords.shape();
        if shape.len() != 3 || shape[2] != 3 {
            return Err(PyValueError::new_err(format!(
                "Expected 3D coordinates array of shape (F, N, 3), found {:?}",
                shape
            )));
        }

        let f_count = shape[0];
        let n_residues = shape[1];

        if f_count == 0 {
            return Err(PyValueError::new_err(
                "Trajectory ensemble must contain at least 1 frame",
            ));
        }
        if n_residues < 4 {
            return Err(PyValueError::new_err(format!(
                "Protein backbone requires at least 4 residues, found {}",
                n_residues
            )));
        }

        let window_radius = (self.window_size / 2).max(1);
        if n_residues < 2 * window_radius + 2 {
            return Err(PyValueError::new_err(format!(
                "Backbone length ({}) is too short for window_size {} (requires at least {} residues)",
                n_residues,
                self.window_size,
                2 * window_radius + 2
            )));
        }

        let array_view = coords.as_array();

        // Process all frames in parallel across CPU cores using Rayon without Python GIL
        let frames: Result<Vec<ConformerFrame>, PyErr> = py.detach(|| {
            (0..f_count)
                .into_par_iter()
                .map(|f_idx| {
                    let frame_slice = array_view.slice(ndarray::s![f_idx, .., ..]);
                    let mut points = Vec::with_capacity(n_residues);
                    for row in frame_slice.rows() {
                        points.push(Point3::new(row[0] as f64, row[1] as f64, row[2] as f64));
                    }
                    let trace = BackboneTrace::new(points);
                    let invariants = extract_curve_invariants(&trace).map_err(to_py_err)?;
                    let writhe_spectrum =
                        compute_local_writhe(&trace, window_radius).map_err(to_py_err)?;

                    Ok(ConformerFrame {
                        frame_id: f_idx as u64,
                        time_ps: f_idx as f64,
                        invariants,
                        writhe_spectrum,
                    })
                })
                .collect()
        });

        let frames = frames?;
        self.inner = Some(RustConformationalIndex::build(frames));
        self.expected_residues = Some(n_residues);

        Ok(())
    }

    /// Queries the index for the k most structurally similar conformers.
    ///
    /// Parameters
    /// ----------
    /// query_coords : numpy.ndarray of shape (N, 3), dtype=float32
    ///     Cartesian coordinates of query C-alpha conformation.
    /// k : int, default=10
    ///     Number of nearest neighbors to retrieve.
    ///
    /// Returns
    /// -------
    /// hits : list of tuple (frame_idx: int, frechet_distance: float)
    ///     Sorted in ascending order of Discrete Fréchet distance.
    #[pyo3(signature = (query_coords, k = 10))]
    pub fn query(
        &self,
        py: Python<'_>,
        query_coords: PyReadonlyArray2<f32>,
        k: usize,
    ) -> PyResult<Vec<(u64, f64)>> {
        let index = self.inner.as_ref().ok_or_else(|| {
            PyValueError::new_err("ConformationalIndex has not been fitted yet. Call fit() first.")
        })?;

        let shape = query_coords.shape();
        if shape.len() != 2 || shape[1] != 3 {
            return Err(PyValueError::new_err(format!(
                "Expected 2D query array of shape (N, 3), found {:?}",
                shape
            )));
        }

        let n = shape[0];
        if let Some(expected) = self.expected_residues {
            if n != expected {
                return Err(PyValueError::new_err(format!(
                    "Residue count mismatch: index was fitted on {} residues, but query has {}",
                    expected, n
                )));
            }
        }

        let window_radius = (self.window_size / 2).max(1);
        let array_view = query_coords.as_array();

        let (target_inv, target_writhe) =
            py.detach(|| -> Result<(CurveInvariants, Vec<f64>), PyErr> {
                let trace = view_to_trace(array_view)?;
                let inv = extract_curve_invariants(&trace).map_err(to_py_err)?;
                let wr = compute_local_writhe(&trace, window_radius).map_err(to_py_err)?;
                Ok((inv, wr))
            })?;

        let coarse_radius = self.coarse_radius;
        let hits =
            py.detach(|| index.query_k_nearest(&target_inv, &target_writhe, k, coarse_radius));

        let results: Vec<(u64, f64)> = hits
            .into_iter()
            .map(|h| (h.frame_id, h.frechet_distance))
            .collect();

        Ok(results)
    }

    /// Queries the index for the k most structurally similar conformers along a localized subcurve / residue range.
    ///
    /// Parameters
    /// ----------
    /// start_res : int
    ///     Starting residue index (0-based, inclusive).
    /// end_res : int
    ///     Ending residue index (0-based, inclusive). Must satisfy `end_res >= start_res + 3`.
    /// query_coords : numpy.ndarray of shape (M, 3), dtype=float32
    ///     Cartesian coordinates of query conformation. May be either the sub-segment of length
    ///     `(end_res - start_res + 1)` or the entire backbone.
    /// k : int, default=10
    ///     Number of nearest neighbors to retrieve.
    ///
    /// Returns
    /// -------
    /// hits : list of tuple (frame_idx: int, frechet_distance: float)
    ///     Sorted in ascending order of Discrete Fréchet distance.
    #[pyo3(signature = (start_res, end_res, query_coords, k = 10))]
    pub fn query_subcurve(
        &self,
        py: Python<'_>,
        start_res: usize,
        end_res: usize,
        query_coords: PyReadonlyArray2<f32>,
        k: usize,
    ) -> PyResult<Vec<(u64, f64)>> {
        let index = self.inner.as_ref().ok_or_else(|| {
            PyValueError::new_err("ConformationalIndex has not been fitted yet. Call fit() first.")
        })?;

        if start_res >= end_res || end_res - start_res < 3 {
            return Err(PyValueError::new_err(format!(
                "Invalid residue range [{}, {}]: subcurve requires at least 4 residues (end_res >= start_res + 3)",
                start_res, end_res
            )));
        }

        if let Some(expected) = self.expected_residues {
            if end_res >= expected {
                return Err(PyValueError::new_err(format!(
                    "end_res ({}) exceeds indexed chain length ({})",
                    end_res, expected
                )));
            }
        }

        let shape = query_coords.shape();
        if shape.len() != 2 || shape[1] != 3 {
            return Err(PyValueError::new_err(format!(
                "Expected 2D query array of shape (N, 3), found {:?}",
                shape
            )));
        }

        let n = shape[0];
        let sub_len = end_res - start_res + 1;
        let array_view = query_coords.as_array();

        let sub_inv = py.detach(|| -> Result<CurveInvariants, PyErr> {
            let trace = view_to_trace(array_view)?;
            if n == sub_len {
                extract_curve_invariants(&trace).map_err(to_py_err)
            } else if n > end_res {
                let full_inv = extract_curve_invariants(&trace).map_err(to_py_err)?;
                full_inv.subcurve(start_res, end_res).ok_or_else(|| {
                    PyValueError::new_err("Failed to extract subcurve invariants")
                })
            } else {
                Err(PyValueError::new_err(format!(
                    "query_coords length ({}) does not match subcurve length ({}) nor full chain length",
                    n, sub_len
                )))
            }
        })?;

        let hits = py.detach(|| index.query_subcurve_k_nearest(start_res, end_res, &sub_inv, k));

        let results: Vec<(u64, f64)> = hits
            .into_iter()
            .map(|h| (h.frame_id, h.frechet_distance))
            .collect();

        Ok(results)
    }

    /// Directly builds a ConformationalIndex from a multi-model PDB trajectory file.
    #[staticmethod]
    #[pyo3(signature = (path, chain = None, window_size = 10, coarse_radius = 0.5))]
    pub fn from_pdb_trajectory(
        py: Python<'_>,
        path: &str,
        chain: Option<char>,
        window_size: usize,
        coarse_radius: f64,
    ) -> PyResult<Self> {
        let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let traces = py.detach(|| {
            topofold_io::read_pdb_trajectory(BufReader::new(file), chain).map_err(to_traj_err)
        })?;
        py.detach(|| traces_to_index(traces, window_size, coarse_radius))
    }

    /// Directly builds a ConformationalIndex from a DCD binary trajectory file.
    #[staticmethod]
    #[pyo3(signature = (path, ca_indices, window_size = 10, coarse_radius = 0.5))]
    pub fn from_dcd(
        py: Python<'_>,
        path: &str,
        ca_indices: Vec<usize>,
        window_size: usize,
        coarse_radius: f64,
    ) -> PyResult<Self> {
        let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let traces = py.detach(|| {
            topofold_io::read_dcd_trajectory(BufReader::new(file), &ca_indices).map_err(to_traj_err)
        })?;
        py.detach(|| traces_to_index(traces, window_size, coarse_radius))
    }

    /// Returns the total number of indexed conformational frames.
    pub fn __len__(&self) -> usize {
        self.inner.as_ref().map_or(0, |idx| idx.len())
    }

    /// String representation.
    pub fn __repr__(&self) -> String {
        format!(
            "ConformationalIndex(frames={}, window_size={}, coarse_radius={})",
            self.__len__(),
            self.window_size,
            self.coarse_radius
        )
    }
}

/// Reads C-alpha coordinates (and optional C-beta coordinates) from a PDB file.
///
/// Parameters
/// ----------
/// path : str
///     Path to PDB file.
/// chain : str, optional
///     Single-character chain ID filter.
/// extract_cbeta : bool, default=False
///     If True, returns a tuple `(ca_coords, cb_coords)`.
///     For Glycine residues lacking C-beta, deterministic pseudo-Cbeta coordinates are computed.
///
/// Returns
/// -------
/// numpy.ndarray of shape (N, 3), dtype=float32, or tuple of two such arrays if `extract_cbeta=True`.
#[pyfunction]
#[pyo3(signature = (path, chain = None, extract_cbeta = false))]
pub fn read_pdb<'py>(
    py: Python<'py>,
    path: &str,
    chain: Option<char>,
    extract_cbeta: bool,
) -> PyResult<Py<PyAny>> {
    let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (trace, _) = py.detach(|| parse_pdb_ca(BufReader::new(file), chain).map_err(to_pdb_err))?;
    let ca_arr = trace_to_array2(&trace);
    let py_ca = PyArray2::from_owned_array(py, ca_arr);
    if extract_cbeta {
        let cb_arr = cb_trace_to_array2(&trace)?;
        let py_cb = PyArray2::from_owned_array(py, cb_arr);
        (py_ca, py_cb).into_py_any(py)
    } else {
        Ok(py_ca.into_any().unbind())
    }
}

/// Reads a multi-model PDB trajectory into a 3D numpy array of shape (F, N, 3).
/// If `extract_cbeta=True`, returns a tuple `(ca_coords, cb_coords)`.
#[pyfunction]
#[pyo3(signature = (path, chain = None, extract_cbeta = false))]
pub fn read_pdb_trajectory<'py>(
    py: Python<'py>,
    path: &str,
    chain: Option<char>,
    extract_cbeta: bool,
) -> PyResult<Py<PyAny>> {
    let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let traces = py.detach(|| {
        topofold_io::read_pdb_trajectory(BufReader::new(file), chain).map_err(to_traj_err)
    })?;
    let ca_arr = traces_to_array3(&traces)?;
    let py_ca = PyArray3::from_owned_array(py, ca_arr);
    if extract_cbeta {
        let cb_arr = cb_traces_to_array3(&traces)?;
        let py_cb = PyArray3::from_owned_array(py, cb_arr);
        (py_ca, py_cb).into_py_any(py)
    } else {
        Ok(py_ca.into_any().unbind())
    }
}

/// Reads a binary DCD trajectory file into a 3D numpy array of shape (F, N, 3).
/// If `cb_indices` is provided, returns `(ca_coords, cb_coords)`.
#[pyfunction]
#[pyo3(signature = (path, ca_indices = None, cb_indices = None))]
pub fn read_dcd<'py>(
    py: Python<'py>,
    path: &str,
    ca_indices: Option<Vec<usize>>,
    cb_indices: Option<Vec<Option<usize>>>,
) -> PyResult<Py<PyAny>> {
    let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let has_cb = cb_indices.is_some();
    let traces = py.detach(|| -> Result<Vec<BackboneTrace>, PyErr> {
        let mut reader = DcdReader::new(BufReader::new(file)).map_err(to_traj_err)?;
        let mut result = Vec::new();
        match (ca_indices.as_ref(), cb_indices.as_ref()) {
            (Some(ca_idxs), Some(cb_idxs)) => {
                while let Some(frame) = reader.next_frame().map_err(to_traj_err)? {
                    let trace = frame
                        .to_backbone_trace_with_cbeta(ca_idxs, cb_idxs)
                        .map_err(to_traj_err)?;
                    result.push(trace);
                }
            }
            (Some(ca_idxs), None) => {
                while let Some(frame) = reader.next_frame().map_err(to_traj_err)? {
                    let trace = frame.to_backbone_trace(ca_idxs).map_err(to_traj_err)?;
                    result.push(trace);
                }
            }
            (None, Some(_)) => {
                return Err(PyValueError::new_err(
                    "cb_indices requires ca_indices to be specified",
                ));
            }
            (None, None) => {
                while let Some(frame) = reader.next_frame().map_err(to_traj_err)? {
                    let trace = frame.to_all_atoms_trace();
                    result.push(trace);
                }
            }
        }
        Ok(result)
    })?;

    let ca_arr = traces_to_array3(&traces)?;
    let py_ca = PyArray3::from_owned_array(py, ca_arr);
    if has_cb {
        let cb_arr = cb_traces_to_array3(&traces)?;
        let py_cb = PyArray3::from_owned_array(py, cb_arr);
        (py_ca, py_cb).into_py_any(py)
    } else {
        Ok(py_ca.into_any().unbind())
    }
}

/// Helper to convert 3D ndarray views of shape (F, N, 3) to Vec<BackboneTrace>.
fn array3_to_traces(
    coords: &ndarray::ArrayView3<'_, f32>,
    cb_coords: Option<&ndarray::ArrayView3<'_, f32>>,
) -> Result<Vec<BackboneTrace>, PyErr> {
    let shape = coords.shape();
    if shape.len() != 3 || shape[2] != 3 {
        return Err(PyValueError::new_err(format!(
            "Expected 3D array of shape (F, N, 3), found {:?}",
            shape
        )));
    }
    let f_count = shape[0];
    let n_residues = shape[1];

    if let Some(cb) = cb_coords {
        if cb.shape() != shape {
            return Err(PyValueError::new_err(format!(
                "cb_coords shape {:?} does not match coords shape {:?}",
                cb.shape(),
                shape
            )));
        }
    }

    let mut traces = Vec::with_capacity(f_count);
    for f in 0..f_count {
        let frame_view = coords.index_axis(ndarray::Axis(0), f);
        let mut ca_points = Vec::with_capacity(n_residues);
        for row in frame_view.rows() {
            ca_points.push(Point3::new(row[0] as f64, row[1] as f64, row[2] as f64));
        }

        if let Some(cb) = cb_coords {
            let cb_frame_view = cb.index_axis(ndarray::Axis(0), f);
            let mut cb_points = Vec::with_capacity(n_residues);
            for row in cb_frame_view.rows() {
                cb_points.push(Point3::new(row[0] as f64, row[1] as f64, row[2] as f64));
            }
            traces.push(BackboneTrace::with_cbeta(ca_points, cb_points));
        } else {
            traces.push(BackboneTrace::new(ca_points));
        }
    }
    Ok(traces)
}

/// Scans a conformational trajectory for bistable cryptic pockets and mobile functional loops without prior hints.
///
/// Slides a window of length `window_size` (default: 8 residues) along the backbone C-alpha trace,
/// evaluates Sarle's Bimodality Coefficient (BC) on intrinsic discrete curve invariants (curvature, torsion,
/// and optional side-chain orientation dihedral theta_beta),
/// identifies continuous residue clusters where BC >= `bc_threshold` (default: 0.60),
/// and ranks candidate segments by peak transition score.
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Trajectory coordinates of F frames across N residues.
/// window_size : int, default=8
///     Length of the sliding subcurve window in residues (must be >= 4).
/// bc_threshold : float, default=0.60
///     Threshold for Sarle's Bimodality Coefficient (values > 0.555 indicate bimodal transitions).
/// cb_coords : numpy.ndarray of shape (F, N, 3), dtype=float32, optional
///     Optional C-beta coordinates of F frames across N residues for rotameric gating detection.
///
/// Returns
/// -------
/// List[Tuple[int, int, float]]
///     List of detected candidate segments as `(start_res, end_res, peak_score)` tuples,
///     sorted by transition score in descending order. Residue indices are 0-based inclusive.
#[pyfunction]
#[pyo3(signature = (coords, window_size = 8, bc_threshold = 0.6, cb_coords = None))]
pub fn scan_cryptic_pockets<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    window_size: usize,
    bc_threshold: f64,
    cb_coords: Option<PyReadonlyArray3<'py, f32>>,
) -> PyResult<Vec<(usize, usize, f64)>> {
    let array_view = coords.as_array();
    let cb_view = cb_coords.as_ref().map(|cb| cb.as_array());
    let traces = array3_to_traces(&array_view, cb_view.as_ref())?;
    let candidates = py.detach(|| detect_bistable_segments(&traces, window_size, bc_threshold));
    let result = candidates
        .into_iter()
        .map(|c| (c.start_res, c.end_res, c.score))
        .collect();
    Ok(result)
}

/// Computes the sliding-window bimodality coefficient profile along the protein sequence.
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Trajectory coordinates of F frames across N residues.
/// window_size : int, default=8
///     Length of the sliding subcurve window in residues.
/// cb_coords : numpy.ndarray of shape (F, N, 3), dtype=float32, optional
///     Optional C-beta coordinates of F frames across N residues.
///
/// Returns
/// -------
/// tuple
///     - If cb_coords is None: `(scores, bc_tau, bc_kappa)` of dtype=float64
///     - If cb_coords is provided: `(scores, bc_tau, bc_kappa, bc_theta)` of dtype=float64
#[pyfunction]
#[pyo3(signature = (coords, window_size = 8, cb_coords = None))]
pub fn compute_bimodality_profile<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    window_size: usize,
    cb_coords: Option<PyReadonlyArray3<'py, f32>>,
) -> PyResult<Py<PyAny>> {
    let array_view = coords.as_array();
    let cb_view = cb_coords.as_ref().map(|cb| cb.as_array());
    let has_cb = cb_coords.is_some();
    let traces = array3_to_traces(&array_view, cb_view.as_ref())?;
    let profile = py.detach(|| -> Result<Vec<topofold_core::WindowBimodality>, PyErr> {
        core_bimodality_profile(&traces, window_size).map_err(to_py_err)
    })?;

    let n = profile.len();
    let mut scores = Vec::with_capacity(n);
    let mut bc_tau = Vec::with_capacity(n);
    let mut bc_kappa = Vec::with_capacity(n);
    let mut bc_theta = Vec::with_capacity(n);

    for w in profile {
        scores.push(w.score);
        bc_tau.push(w.bc_tau);
        bc_kappa.push(w.bc_kappa);
        bc_theta.push(w.bc_theta);
    }

    let py_scores = PyArray1::from_vec(py, scores).unbind();
    let py_tau = PyArray1::from_vec(py, bc_tau).unbind();
    let py_kappa = PyArray1::from_vec(py, bc_kappa).unbind();

    if has_cb {
        let py_theta = PyArray1::from_vec(py, bc_theta).unbind();
        (py_scores, py_tau, py_kappa, py_theta).into_py_any(py)
    } else {
        (py_scores, py_tau, py_kappa).into_py_any(py)
    }
}

/// Computes the Intrinsic Allosteric Communication Network matrix across a trajectory.
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Trajectory coordinates of F frames across N residues.
/// cb_coords : numpy.ndarray of shape (F, N, 3), dtype=float32, optional
///     Optional C-beta coordinates for side-chain rotamer orientation.
/// regularizer_eps : float, default=1e-6
///     Diagonal regularizer jitter added to covariance matrices.
///
/// Returns
/// -------
/// numpy.ndarray of shape (N, N), dtype=float64
///     Symmetric matrix of generalized correlation coefficients r_MI in [0, 1].
#[pyfunction]
#[pyo3(signature = (coords, cb_coords = None, regularizer_eps = 1e-6))]
pub fn compute_allosteric_network<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    cb_coords: Option<PyReadonlyArray3<'py, f32>>,
    regularizer_eps: f64,
) -> PyResult<Py<PyArray2<f64>>> {
    let array_view = coords.as_array();
    let cb_view = cb_coords.as_ref().map(|cb| cb.as_array());
    let traces = array3_to_traces(&array_view, cb_view.as_ref())?;
    let matrix = py.detach(|| {
        compute_intrinsic_allosteric_network(&traces, regularizer_eps).map_err(to_py_err)
    })?;
    let py_arr = PyArray2::from_owned_array(py, matrix);
    Ok(py_arr.unbind())
}

/// Computes the N_A x N_B Inter-Molecular Allosteric Communication Network matrix between two chains.
///
/// Parameters
/// ----------
/// coords_a : numpy.ndarray of shape (F, N_A, 3), dtype=float32
///     Trajectory coordinates of Chain A across F frames.
/// coords_b : numpy.ndarray of shape (F, N_B, 3), dtype=float32
///     Trajectory coordinates of Chain B across F frames.
/// cb_coords_a : numpy.ndarray of shape (F, N_A, 3), dtype=float32, optional
///     Optional C-beta coordinates for Chain A.
/// cb_coords_b : numpy.ndarray of shape (F, N_B, 3), dtype=float32, optional
///     Optional C-beta coordinates for Chain B.
/// regularizer_eps : float, default=1e-6
///     Diagonal regularizer jitter added to covariance matrices.
///
/// Returns
/// -------
/// numpy.ndarray of shape (N_A, N_B), dtype=float64
///     Cross-chain matrix of generalized correlation coefficients r_MI in [0, 1].
#[pyfunction]
#[pyo3(signature = (coords_a, coords_b, cb_coords_a = None, cb_coords_b = None, regularizer_eps = 1e-6))]
pub fn compute_intermolecular_allosteric_network<'py>(
    py: Python<'py>,
    coords_a: PyReadonlyArray3<'py, f32>,
    coords_b: PyReadonlyArray3<'py, f32>,
    cb_coords_a: Option<PyReadonlyArray3<'py, f32>>,
    cb_coords_b: Option<PyReadonlyArray3<'py, f32>>,
    regularizer_eps: f64,
) -> PyResult<Py<PyArray2<f64>>> {
    let a_view = coords_a.as_array();
    let b_view = coords_b.as_array();
    let cb_a_view = cb_coords_a.as_ref().map(|cb| cb.as_array());
    let cb_b_view = cb_coords_b.as_ref().map(|cb| cb.as_array());

    let traces_a = array3_to_traces(&a_view, cb_a_view.as_ref())?;
    let traces_b = array3_to_traces(&b_view, cb_b_view.as_ref())?;

    let matrix = py.detach(|| {
        core_intermolecular_allosteric_network(&traces_a, &traces_b, regularizer_eps)
            .map_err(to_py_err)
    })?;
    let py_arr = PyArray2::from_owned_array(py, matrix);
    Ok(py_arr.unbind())
}

/// Computes the average dynamic allosteric cooperativity index across an inter-protein contact interface.
///
/// Parameters
/// ----------
/// inter_matrix : numpy.ndarray of shape (N_A, N_B), dtype=float64
///     The cross-chain generalized correlation matrix r_MI.
/// interface_pairs : list of tuple of int (i_A, j_B)
///     List of 0-based residue pairs located at the physical contact interface.
///
/// Returns
/// -------
/// float
///     The mean dynamic cooperativity score across the interface contacts in [0, 1].
#[pyfunction]
pub fn compute_ternary_cooperativity_index(
    inter_matrix: PyReadonlyArray2<'_, f64>,
    interface_pairs: Vec<(usize, usize)>,
) -> PyResult<f64> {
    let arr = inter_matrix.as_array();
    let score = core_ternary_cooperativity_index(&arr.to_owned(), &interface_pairs);
    Ok(score)
}

/// Computes the Discrete Fréchet Distance between two 3D polygonal curves P and Q.
///
/// Parameters
/// ----------
/// coords_a : numpy.ndarray of shape (N, 3), dtype=float32
///     First curve coordinates.
/// coords_b : numpy.ndarray of shape (M, 3), dtype=float32
///     Second curve coordinates.
///
/// Returns
/// -------
/// float
///     Minimax coupling discrete Fréchet distance in Ångströms.
#[pyfunction]
pub fn discrete_frechet_distance(
    coords_a: PyReadonlyArray2<'_, f32>,
    coords_b: PyReadonlyArray2<'_, f32>,
) -> PyResult<f64> {
    let a_view = coords_a.as_array();
    let b_view = coords_b.as_array();
    let trace_a = view_to_trace(a_view)?;
    let trace_b = view_to_trace(b_view)?;
    let dist = discrete_frechet_distance_coords(trace_a.coordinates(), trace_b.coordinates());
    Ok(dist)
}

/// Computes the Spectral Topological Density S_topo across an ensemble of IDP conformations.
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Trajectory coordinates of F frames across N residues.
/// window_radius : int, default=4
///     Half-width of the sliding window in residues (e.g. 4 for 9-residue subchains).
///
/// Returns
/// -------
/// numpy.ndarray of shape (N,), dtype=float64
///     Mean topological compactness S_topo(i) across all frames.
#[pyfunction]
#[pyo3(signature = (coords, window_radius = 4))]
pub fn compute_idp_topological_density<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    window_radius: usize,
) -> PyResult<Py<PyArray1<f64>>> {
    let array_view = coords.as_array();
    let traces = array3_to_traces(&array_view, None)?;
    let profile =
        py.detach(|| compute_idp_density_profile(&traces, window_radius).map_err(to_py_err))?;
    let py_arr = PyArray1::from_vec(py, profile.mean_density);
    Ok(py_arr.unbind())
}

/// Return type for the full Spectral Topological Density profile: (mean, variance, std, z_scores).
type IdpProfileArrays = (
    Py<PyArray1<f64>>,
    Py<PyArray1<f64>>,
    Py<PyArray1<f64>>,
    Py<PyArray1<f64>>,
);

/// Transient motif tuple: (start_res, end_res, peak_res, mean_compactness, peak_compactness, z_score).
type IdpMotifTuple = (usize, usize, usize, f64, f64, f64);

/// Computes the complete Spectral Topological Density profile (mean, variance, std, z-scores).
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Trajectory coordinates of F frames across N residues.
/// window_radius : int, default=4
///     Half-width of the sliding window in residues (e.g. 4 for 9-residue subchains).
///
/// Returns
/// -------
/// tuple of 4 numpy.ndarrays of shape (N,), dtype=float64
///     `(mean_density, variance_density, std_density, z_scores)`
#[pyfunction(name = "compute_idp_density_profile")]
#[pyo3(signature = (coords, window_radius = 4))]
pub fn compute_idp_density_profile_py<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    window_radius: usize,
) -> PyResult<IdpProfileArrays> {
    let array_view = coords.as_array();
    let traces = array3_to_traces(&array_view, None)?;
    let profile =
        py.detach(|| compute_idp_density_profile(&traces, window_radius).map_err(to_py_err))?;
    let py_mean = PyArray1::from_vec(py, profile.mean_density).unbind();
    let py_var = PyArray1::from_vec(py, profile.variance_density).unbind();
    let py_std = PyArray1::from_vec(py, profile.std_density).unbind();
    let py_z = PyArray1::from_vec(py, profile.z_scores).unbind();
    Ok((py_mean, py_var, py_std, py_z))
}

/// Autonomously detects transiently structured pre-nucleation motifs within an IDP ensemble.
///
/// Parameters
/// ----------
/// coords : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Trajectory coordinates of F frames across N residues.
/// window_size : int, default=8
///     Width of the sliding window in residues.
/// z_threshold : float, default=1.0
///     Minimum sequence Z-score to include a residue in a motif.
///
/// Returns
/// -------
/// list of tuples: `[(start_res, end_res, peak_res, mean_compactness, peak_compactness, z_score), ...]`
#[pyfunction(name = "detect_transient_motifs")]
#[pyo3(signature = (coords, window_size = 8, z_threshold = 1.0))]
pub fn detect_transient_motifs_py<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    window_size: usize,
    z_threshold: f64,
) -> PyResult<Vec<IdpMotifTuple>> {
    let array_view = coords.as_array();
    let traces = array3_to_traces(&array_view, None)?;
    let window_radius = (window_size / 2).max(1);
    let motifs = py.detach(|| {
        detect_transient_motifs_with_params(&traces, window_radius, z_threshold, 2)
            .map_err(to_py_err)
    })?;
    let result = motifs
        .into_iter()
        .map(|m| {
            (
                m.start_res,
                m.end_res,
                m.peak_res,
                m.mean_compactness,
                m.peak_compactness,
                m.z_score,
            )
        })
        .collect();
    Ok(result)
}

/// Python wrapper for an RNA ribonucleic ribbon trace.
#[pyclass(name = "RnaRibbonTrace", from_py_object)]
#[derive(Clone)]
pub struct PyRnaRibbonTrace {
    pub(crate) inner: topofold_core::rna::RnaRibbonTrace,
}

#[pymethods]
impl PyRnaRibbonTrace {
    /// 3D coordinates of Phosphorus atoms as (N, 3) float32 numpy array.
    #[getter]
    pub fn p_coords<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f32>>> {
        let n = self.inner.len();
        let mut arr = ndarray::Array2::<f32>::zeros((n, 3));
        for (i, pt) in self.inner.p_coords().iter().enumerate() {
            arr[[i, 0]] = pt.x as f32;
            arr[[i, 1]] = pt.y as f32;
            arr[[i, 2]] = pt.z as f32;
        }
        Ok(PyArray2::from_owned_array(py, arr).unbind())
    }

    /// 3D coordinates of ribose C4' atoms as (N, 3) float32 numpy array.
    #[getter]
    pub fn c4_coords<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f32>>> {
        let n = self.inner.len();
        let mut arr = ndarray::Array2::<f32>::zeros((n, 3));
        for (i, pt) in self.inner.c4_coords().iter().enumerate() {
            arr[[i, 0]] = pt.x as f32;
            arr[[i, 1]] = pt.y as f32;
            arr[[i, 2]] = pt.z as f32;
        }
        Ok(PyArray2::from_owned_array(py, arr).unbind())
    }

    /// 3D coordinates of ribose C1' atoms as (N, 3) float32 numpy array.
    #[getter]
    pub fn c1_coords<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f32>>> {
        let n = self.inner.len();
        let mut arr = ndarray::Array2::<f32>::zeros((n, 3));
        for (i, pt) in self.inner.c1_coords().iter().enumerate() {
            arr[[i, 0]] = pt.x as f32;
            arr[[i, 1]] = pt.y as f32;
            arr[[i, 2]] = pt.z as f32;
        }
        Ok(PyArray2::from_owned_array(py, arr).unbind())
    }

    /// 3D coordinates of glycosidic nitrogen atoms as (N, 3) float32 numpy array.
    #[getter]
    pub fn n_coords<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f32>>> {
        let n = self.inner.len();
        let mut arr = ndarray::Array2::<f32>::zeros((n, 3));
        for (i, pt) in self.inner.n_coords().iter().enumerate() {
            arr[[i, 0]] = pt.x as f32;
            arr[[i, 1]] = pt.y as f32;
            arr[[i, 2]] = pt.z as f32;
        }
        Ok(PyArray2::from_owned_array(py, arr).unbind())
    }

    /// Base orientation unit vectors (pointing from C1' to N) as (N, 3) float32 numpy array.
    #[getter]
    pub fn base_vectors<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f32>>> {
        let vecs = self.inner.base_unit_vectors().map_err(to_py_err)?;
        let n = vecs.len();
        let mut arr = ndarray::Array2::<f32>::zeros((n, 3));
        for (i, v) in vecs.iter().enumerate() {
            arr[[i, 0]] = v.x as f32;
            arr[[i, 1]] = v.y as f32;
            arr[[i, 2]] = v.z as f32;
        }
        Ok(PyArray2::from_owned_array(py, arr).unbind())
    }

    /// Nucleotide residue names (e.g. ['C', 'G', ...]).
    #[getter]
    pub fn names(&self) -> Vec<String> {
        self.inner.names().to_vec()
    }

    /// Residue sequence numbers.
    #[getter]
    pub fn seq_ids(&self) -> Vec<i32> {
        self.inner.seq_ids().to_vec()
    }

    /// Chain identifiers.
    #[getter]
    pub fn chain_ids(&self) -> Vec<String> {
        self.inner
            .chain_ids()
            .iter()
            .map(|c| c.to_string())
            .collect()
    }

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn __repr__(&self) -> String {
        format!("<RnaRibbonTrace: {} nucleotides>", self.inner.len())
    }
}

/// Reads an RNA structure from a PDB file into an SE(3)-invariant [`RnaRibbonTrace`].
///
/// Parameters
/// ----------
/// path : str
///     Path to the input PDB file.
/// chain : str or None, optional
///     Single-character chain ID filter (e.g. 'X'). If None, parses all chains.
///
/// Returns
/// -------
/// RnaRibbonTrace
///     Oriented ribonucleic ribbon trace containing Phosphorus and base coordinates.
#[pyfunction]
#[pyo3(signature = (path, chain = None))]
pub fn read_pdb_rna(path: &str, chain: Option<char>) -> PyResult<PyRnaRibbonTrace> {
    let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (trace, _) = parse_pdb_rna(BufReader::new(file), chain).map_err(to_pdb_err)?;
    Ok(PyRnaRibbonTrace { inner: trace })
}

/// Computes SE(3)-invariant ribonucleic curve invariants (kappa_P, tau_P, Writhe_P, theta_base).
///
/// Parameters
/// ----------
/// p_coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Coordinates of N Phosphorus atoms in Ångströms.
/// base_coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Coordinates of N base orientation vectors or glycosidic nitrogen atoms.
/// window_radius : int, default=2
///     Sliding window half-width for localized Gauss linking writhe spectrum.
///
/// Returns
type RnaInvariantsTuple = (
    Py<PyArray1<f64>>,
    Py<PyArray1<f64>>,
    Py<PyArray1<f64>>,
    Py<PyArray1<f64>>,
);

/// Computes SE(3)-invariant ribonucleic curve invariants (kappa_P, tau_P, Writhe_P, theta_base).
///
/// Parameters
/// ----------
/// p_coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Coordinates of N Phosphorus atoms in Ångströms.
/// base_coords : numpy.ndarray of shape (N, 3), dtype=float32
///     Coordinates of N base orientation vectors or glycosidic nitrogen atoms.
/// window_radius : int, default=2
///     Sliding window half-width for localized Gauss linking writhe spectrum.
///
/// Returns
/// -------
/// tuple of 4 numpy.ndarrays:
///     - kappa : (N - 2,) float64 (Phosphorus curvature turning angles in [0, pi])
///     - tau : (N - 3,) float64 (Phosphorus torsion dihedrals in (-pi, pi])
///     - writhe_spectrum : (N - 2 * window_radius,) float64 (Local Phosphorus writhe)
///     - theta_base : (N - 2,) float64 (Glycosidic base ribbon orientation dihedrals in (-pi, pi])
#[pyfunction]
#[pyo3(signature = (p_coords, base_coords, window_radius = 2))]
pub fn compute_rna_invariants<'py>(
    py: Python<'py>,
    p_coords: PyReadonlyArray2<'py, f32>,
    base_coords: PyReadonlyArray2<'py, f32>,
    window_radius: usize,
) -> PyResult<RnaInvariantsTuple> {
    let p_view = p_coords.as_array();
    let b_view = base_coords.as_array();

    if p_view.shape() != b_view.shape() {
        return Err(PyValueError::new_err(format!(
            "p_coords shape {:?} does not match base_coords shape {:?}",
            p_view.shape(),
            b_view.shape()
        )));
    }
    if p_view.shape()[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 coordinates"));
    }
    let n = p_view.shape()[0];
    if n < 4 {
        return Err(PyValueError::new_err(format!(
            "RNA ribbon requires at least 4 nucleotides, found {}",
            n
        )));
    }

    let p_pts: Vec<Point3<f64>> = p_view
        .outer_iter()
        .map(|r| Point3::new(r[0] as f64, r[1] as f64, r[2] as f64))
        .collect();
    let c1_pts: Vec<Point3<f64>> = vec![Point3::new(0.0, 0.0, 0.0); n];
    let n_pts: Vec<Point3<f64>> = b_view
        .outer_iter()
        .map(|r| Point3::new(r[0] as f64, r[1] as f64, r[2] as f64))
        .collect();
    let c4_pts = p_pts.clone();

    let trace = RnaRibbonTrace::new(p_pts, c4_pts, c1_pts, n_pts).map_err(to_py_err)?;
    let invariants = py.detach(|| {
        extract_rna_ribbon_invariants_with_window(&trace, window_radius).map_err(to_py_err)
    })?;

    let py_kappa = PyArray1::from_vec(py, invariants.curvatures).unbind();
    let py_tau = PyArray1::from_vec(py, invariants.torsions).unbind();
    let py_writhe = PyArray1::from_vec(py, invariants.local_writhes).unbind();
    let py_theta = PyArray1::from_vec(py, invariants.base_dihedrals).unbind();

    Ok((py_kappa, py_tau, py_writhe, py_theta))
}

/// Scans an RNA trajectory ensemble for bistable conformational switching hinges.
///
/// Parameters
/// ----------
/// coords_ensemble : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Trajectory of Phosphorus coordinates across F frames.
/// base_ensemble : numpy.ndarray of shape (F, N, 3), dtype=float32, optional
///     Optional base orientation vectors across F frames.
/// window_size : int, default=3
///     Sliding window width in nucleotides.
/// threshold : float, default=0.70
///     Minimum Sarle's Bimodality Coefficient cutoff to report a candidate hinge.
///
/// Returns
/// -------
/// list of dict
///     Ranked candidate switching hinges: `[{'start': int, 'end': int, 'score': float, 'bc_tau': float, 'bc_kappa': float, 'bc_theta': float, 'peak_res': int}, ...]`
#[pyfunction]
#[pyo3(signature = (coords_ensemble, base_ensemble = None, window_size = 3, threshold = 0.70))]
pub fn scan_rna_switching_hinges<'py>(
    py: Python<'py>,
    coords_ensemble: PyReadonlyArray3<'py, f32>,
    base_ensemble: Option<PyReadonlyArray3<'py, f32>>,
    window_size: usize,
    threshold: f64,
) -> PyResult<Vec<Bound<'py, pyo3::types::PyDict>>> {
    let p_arr = coords_ensemble.as_array();
    let shape = p_arr.shape();
    let f_count = shape[0];
    let n_nt = shape[1];

    let b_view = base_ensemble.as_ref().map(|b| b.as_array());

    let mut traces = Vec::with_capacity(f_count);
    for f in 0..f_count {
        let mut p_pts = Vec::with_capacity(n_nt);
        let mut c1_pts = Vec::with_capacity(n_nt);
        let mut n_pts = Vec::with_capacity(n_nt);
        let mut c4_pts = Vec::with_capacity(n_nt);

        for i in 0..n_nt {
            let px = p_arr[[f, i, 0]] as f64;
            let py_val = p_arr[[f, i, 1]] as f64;
            let pz = p_arr[[f, i, 2]] as f64;
            p_pts.push(Point3::new(px, py_val, pz));
            c4_pts.push(Point3::new(px, py_val, pz));
            c1_pts.push(Point3::new(0.0, 0.0, 0.0));

            if let Some(ref b_arr) = b_view {
                let bx = b_arr[[f, i, 0]] as f64;
                let by = b_arr[[f, i, 1]] as f64;
                let bz = b_arr[[f, i, 2]] as f64;
                n_pts.push(Point3::new(bx, by, bz));
            } else {
                n_pts.push(Point3::new(px + 1.0, py_val, pz));
            }
        }

        let trace = RnaRibbonTrace::new(p_pts, c4_pts, c1_pts, n_pts).map_err(to_py_err)?;
        traces.push(trace);
    }

    let candidates = py.detach(|| {
        core_scan_rna_switching_hinges(&traces, window_size, threshold).map_err(to_py_err)
    })?;

    let mut result = Vec::with_capacity(candidates.len());
    for c in candidates {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("start", c.start_res)?;
        dict.set_item("end", c.end_res)?;
        dict.set_item("score", c.score)?;
        dict.set_item("bc_tau", c.bc_tau)?;
        dict.set_item("bc_kappa", c.bc_kappa)?;
        dict.set_item("bc_theta", c.bc_theta)?;
        dict.set_item("peak_res", c.peak_res)?;
        result.push(dict);
    }

    Ok(result)
}

/// Computes sequence-wide Sarle's bimodality profile for an RNA trajectory.
///
/// Parameters
/// ----------
/// coords_ensemble : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Phosphorus coordinates across F frames.
/// base_ensemble : numpy.ndarray of shape (F, N, 3), dtype=float32, optional
///     Optional base orientation unit vectors across F frames.
/// window_size : int, default=3
///     Sliding window width in nucleotides.
///
/// Returns
/// -------
/// tuple of numpy.ndarray
type PyRnaBimodalityProfile<'py> = (
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
);

/// Computes sequence-wide Sarle's bimodality profile for an RNA trajectory.
///
/// Parameters
/// ----------
/// coords_ensemble : numpy.ndarray of shape (F, N, 3), dtype=float32
///     Phosphorus coordinates across F frames.
/// base_ensemble : numpy.ndarray of shape (F, N, 3), dtype=float32, optional
///     Optional base orientation unit vectors across F frames.
/// window_size : int, default=3
///     Sliding window width in nucleotides.
///
/// Returns
/// -------
/// tuple of numpy.ndarray
///     `(scores, bc_theta, bc_kappa, bc_tau)`:
///     - `scores`: composite max(BC_theta, BC_kappa, BC_tau) (1D float64)
///     - `bc_theta`: glycosidic ribbon dihedral bimodality (1D float64)
///     - `bc_kappa`: backbone curvature bimodality (1D float64)
///     - `bc_tau`: backbone torsion bimodality (1D float64)
#[pyfunction]
#[pyo3(signature = (coords_ensemble, base_ensemble = None, window_size = 3))]
pub fn compute_rna_bimodality_profile<'py>(
    py: Python<'py>,
    coords_ensemble: PyReadonlyArray3<'py, f32>,
    base_ensemble: Option<PyReadonlyArray3<'py, f32>>,
    window_size: usize,
) -> PyResult<PyRnaBimodalityProfile<'py>> {
    let p_arr = coords_ensemble.as_array();
    let shape = p_arr.shape();
    let f_count = shape[0];
    let n_nt = shape[1];

    let b_view = base_ensemble.as_ref().map(|b| b.as_array());

    let mut traces = Vec::with_capacity(f_count);
    for f in 0..f_count {
        let mut p_pts = Vec::with_capacity(n_nt);
        let mut c1_pts = Vec::with_capacity(n_nt);
        let mut n_pts = Vec::with_capacity(n_nt);
        let mut c4_pts = Vec::with_capacity(n_nt);

        for i in 0..n_nt {
            let px = p_arr[[f, i, 0]] as f64;
            let py_val = p_arr[[f, i, 1]] as f64;
            let pz = p_arr[[f, i, 2]] as f64;
            p_pts.push(Point3::new(px, py_val, pz));
            c4_pts.push(Point3::new(px, py_val, pz));
            c1_pts.push(Point3::new(0.0, 0.0, 0.0));

            if let Some(ref b_arr) = b_view {
                let bx = b_arr[[f, i, 0]] as f64;
                let by = b_arr[[f, i, 1]] as f64;
                let bz = b_arr[[f, i, 2]] as f64;
                n_pts.push(Point3::new(bx, by, bz));
            } else {
                n_pts.push(Point3::new(px + 1.0, py_val, pz));
            }
        }

        let trace = RnaRibbonTrace::new(p_pts, c4_pts, c1_pts, n_pts).map_err(to_py_err)?;
        traces.push(trace);
    }

    let (scores, bc_th, bc_k, bc_t) =
        py.detach(|| core_rna_bimodality_profile(&traces, window_size).map_err(to_py_err))?;

    Ok((
        PyArray1::from_vec(py, scores),
        PyArray1::from_vec(py, bc_th),
        PyArray1::from_vec(py, bc_k),
        PyArray1::from_vec(py, bc_t),
    ))
}

/// TopoFold Python module definition.
#[pymodule]
fn topofold(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", "High-performance SE(3)-invariant discrete differential geometry engine for protein trajectories")?;
    m.add_function(wrap_pyfunction!(compute_invariants, m)?)?;
    m.add_function(wrap_pyfunction!(compute_theta_beta, m)?)?;
    m.add_function(wrap_pyfunction!(scan_cryptic_pockets, m)?)?;
    m.add_function(wrap_pyfunction!(compute_bimodality_profile, m)?)?;
    m.add_function(wrap_pyfunction!(compute_allosteric_network, m)?)?;
    m.add_function(wrap_pyfunction!(
        compute_intermolecular_allosteric_network,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(compute_ternary_cooperativity_index, m)?)?;
    m.add_function(wrap_pyfunction!(compute_idp_topological_density, m)?)?;
    m.add_function(wrap_pyfunction!(compute_idp_density_profile_py, m)?)?;
    m.add_function(wrap_pyfunction!(detect_transient_motifs_py, m)?)?;
    m.add_function(wrap_pyfunction!(discrete_frechet_distance, m)?)?;
    m.add_function(wrap_pyfunction!(read_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(read_pdb_trajectory, m)?)?;
    m.add_function(wrap_pyfunction!(read_dcd, m)?)?;
    m.add_function(wrap_pyfunction!(read_pdb_rna, m)?)?;
    m.add_function(wrap_pyfunction!(compute_rna_invariants, m)?)?;
    m.add_function(wrap_pyfunction!(scan_rna_switching_hinges, m)?)?;
    m.add_function(wrap_pyfunction!(compute_rna_bimodality_profile, m)?)?;
    m.add_class::<ConformationalIndex>()?;
    m.add_class::<PyRnaRibbonTrace>()?;
    Ok(())
}

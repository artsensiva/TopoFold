//! Python bindings (PyO3) for the TopoFold conformational trajectory engine.

use std::fs::File;
use std::io::BufReader;
use nalgebra::Point3;
use numpy::{PyArray1, PyArray2, PyArray3, PyReadonlyArray2, PyReadonlyArray3, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

use topofold_core::error::GeometryError;
use topofold_core::invariants::CurveInvariants;
use topofold_core::pdb::{parse_pdb_ca, PdbError};
use topofold_core::{
    compute_bimodality_profile as core_bimodality_profile, compute_local_writhe,
    detect_bistable_segments, extract_curve_invariants, BackboneTrace,
};
use topofold_index::{ConformerFrame, ConformationalIndex as RustConformationalIndex};
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
///
/// Returns
/// -------
/// (kappa, tau, writhe_spectrum) : tuple of numpy.ndarray of dtype=float64
///     - kappa: Discrete turning angles in [0, pi] (length N - 2)
///     - tau: Discrete signed dihedrals in (-pi, pi] (length N - 3)
///     - writhe_spectrum: Localized sliding-window writhe (length N - 2 * window_radius)
/// Type alias for the tuple returned by compute_invariants.
pub type InvariantsTuple = (Py<PyArray1<f64>>, Py<PyArray1<f64>>, Py<PyArray1<f64>>);

#[pyfunction]
pub fn compute_invariants<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f32>,
) -> PyResult<InvariantsTuple> {
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

    let (invariants, writhe_spectrum) = py.detach(|| -> Result<(CurveInvariants, Vec<f64>), PyErr> {
        let trace = view_to_trace(array_view)?;
        let inv = extract_curve_invariants(&trace).map_err(to_py_err)?;
        let wr = compute_local_writhe(&trace, window_radius).map_err(to_py_err)?;
        Ok((inv, wr))
    })?;

    let py_kappa = PyArray1::from_vec(py, invariants.curvatures).unbind();
    let py_tau = PyArray1::from_vec(py, invariants.torsions).unbind();
    let py_writhe = PyArray1::from_vec(py, writhe_spectrum).unbind();

    Ok((py_kappa, py_tau, py_writhe))
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
            return Err(PyValueError::new_err("Trajectory ensemble must contain at least 1 frame"));
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
                    let writhe_spectrum = compute_local_writhe(&trace, window_radius).map_err(to_py_err)?;

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

        let (target_inv, target_writhe) = py.detach(|| -> Result<(CurveInvariants, Vec<f64>), PyErr> {
            let trace = view_to_trace(array_view)?;
            let inv = extract_curve_invariants(&trace).map_err(to_py_err)?;
            let wr = compute_local_writhe(&trace, window_radius).map_err(to_py_err)?;
            Ok((inv, wr))
        })?;

        let coarse_radius = self.coarse_radius;
        let hits = py.detach(|| {
            index.query_k_nearest(&target_inv, &target_writhe, k, coarse_radius)
        });

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

        let hits = py.detach(|| {
            index.query_subcurve_k_nearest(start_res, end_res, &sub_inv, k)
        });

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
        let traces = py.detach(|| topofold_io::read_pdb_trajectory(BufReader::new(file), chain).map_err(to_traj_err))?;
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
        let traces = py.detach(|| topofold_io::read_dcd_trajectory(BufReader::new(file), &ca_indices).map_err(to_traj_err))?;
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

/// Reads C-alpha coordinates from a PDB file into a 2D numpy array of shape (N, 3).
#[pyfunction]
#[pyo3(signature = (path, chain = None))]
pub fn read_pdb<'py>(
    py: Python<'py>,
    path: &str,
    chain: Option<char>,
) -> PyResult<Py<PyArray2<f32>>> {
    let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (trace, _) = py.detach(|| parse_pdb_ca(BufReader::new(file), chain).map_err(to_pdb_err))?;
    let arr = trace_to_array2(&trace);
    Ok(PyArray2::from_owned_array(py, arr).unbind())
}

/// Reads a multi-model PDB trajectory into a 3D numpy array of shape (F, N, 3).
#[pyfunction]
#[pyo3(signature = (path, chain = None))]
pub fn read_pdb_trajectory<'py>(
    py: Python<'py>,
    path: &str,
    chain: Option<char>,
) -> PyResult<Py<PyArray3<f32>>> {
    let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let traces = py.detach(|| topofold_io::read_pdb_trajectory(BufReader::new(file), chain).map_err(to_traj_err))?;
    let arr = traces_to_array3(&traces)?;
    Ok(PyArray3::from_owned_array(py, arr).unbind())
}

/// Reads a binary DCD trajectory file into a 3D numpy array of shape (F, N, 3).
#[pyfunction]
#[pyo3(signature = (path, ca_indices = None))]
pub fn read_dcd<'py>(
    py: Python<'py>,
    path: &str,
    ca_indices: Option<Vec<usize>>,
) -> PyResult<Py<PyArray3<f32>>> {
    let file = File::open(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let traces = py.detach(|| -> Result<Vec<BackboneTrace>, PyErr> {
        let mut reader = DcdReader::new(BufReader::new(file)).map_err(to_traj_err)?;
        let mut result = Vec::new();
        while let Some(frame) = reader.next_frame().map_err(to_traj_err)? {
            let trace = match &ca_indices {
                Some(indices) => frame.to_backbone_trace(indices).map_err(to_traj_err)?,
                None => frame.to_all_atoms_trace(),
            };
            result.push(trace);
        }
        Ok(result)
    })?;
    let arr = traces_to_array3(&traces)?;
    Ok(PyArray3::from_owned_array(py, arr).unbind())
}

/// Helper to convert a 3D ndarray view of shape (F, N, 3) to Vec<BackboneTrace>.
fn array3_to_traces(coords: &ndarray::ArrayView3<'_, f32>) -> Result<Vec<BackboneTrace>, PyErr> {
    let shape = coords.shape();
    if shape.len() != 3 || shape[2] != 3 {
        return Err(PyValueError::new_err(format!(
            "Expected 3D array of shape (F, N, 3), found {:?}",
            shape
        )));
    }
    let f_count = shape[0];
    let mut traces = Vec::with_capacity(f_count);
    for f in 0..f_count {
        let frame_view = coords.index_axis(ndarray::Axis(0), f);
        traces.push(view_to_trace(frame_view)?);
    }
    Ok(traces)
}

/// Scans a conformational trajectory for bistable cryptic pockets and mobile functional loops without prior hints.
///
/// Slides a window of length `window_size` (default: 8 residues) along the backbone C-alpha trace,
/// evaluates Sarle's Bimodality Coefficient (BC) on intrinsic discrete curve invariants (curvature, torsion),
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
///
/// Returns
/// -------
/// List[Tuple[int, int, float]]
///     List of detected candidate segments as `(start_res, end_res, peak_score)` tuples,
///     sorted by transition score in descending order. Residue indices are 0-based inclusive.
#[pyfunction]
#[pyo3(signature = (coords, window_size = 8, bc_threshold = 0.6))]
pub fn scan_cryptic_pockets<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    window_size: usize,
    bc_threshold: f64,
) -> PyResult<Vec<(usize, usize, f64)>> {
    let array_view = coords.as_array();
    let traces = array3_to_traces(&array_view)?;
    let candidates = py.detach(|| {
        detect_bistable_segments(&traces, window_size, bc_threshold)
    });
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
///
/// Returns
/// -------
/// (scores, bc_tau, bc_kappa) : tuple of numpy.ndarray of dtype=float64
///     Arrays of length N - window_size + 1 containing:
///     - scores: Composite bimodality score max(bc_tau, bc_kappa)
///     - bc_tau: Sarle's bimodality coefficient for discrete torsion tau
///     - bc_kappa: Sarle's bimodality coefficient for discrete curvature kappa
#[pyfunction]
#[pyo3(signature = (coords, window_size = 8))]
pub fn compute_bimodality_profile<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray3<'py, f32>,
    window_size: usize,
) -> PyResult<(Py<PyArray1<f64>>, Py<PyArray1<f64>>, Py<PyArray1<f64>>)> {
    let array_view = coords.as_array();
    let traces = array3_to_traces(&array_view)?;
    let profile = py.detach(|| -> Result<Vec<topofold_core::WindowBimodality>, PyErr> {
        core_bimodality_profile(&traces, window_size).map_err(to_py_err)
    })?;

    let n = profile.len();
    let mut scores = Vec::with_capacity(n);
    let mut bc_tau = Vec::with_capacity(n);
    let mut bc_kappa = Vec::with_capacity(n);

    for w in profile {
        scores.push(w.score);
        bc_tau.push(w.bc_tau);
        bc_kappa.push(w.bc_kappa);
    }

    let py_scores = PyArray1::from_vec(py, scores).unbind();
    let py_tau = PyArray1::from_vec(py, bc_tau).unbind();
    let py_kappa = PyArray1::from_vec(py, bc_kappa).unbind();

    Ok((py_scores, py_tau, py_kappa))
}

/// TopoFold Python module definition.
#[pymodule]
fn topofold(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__doc__", "High-performance SE(3)-invariant discrete differential geometry engine for protein trajectories")?;
    m.add_function(wrap_pyfunction!(compute_invariants, m)?)?;
    m.add_function(wrap_pyfunction!(scan_cryptic_pockets, m)?)?;
    m.add_function(wrap_pyfunction!(compute_bimodality_profile, m)?)?;
    m.add_function(wrap_pyfunction!(read_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(read_pdb_trajectory, m)?)?;
    m.add_function(wrap_pyfunction!(read_dcd, m)?)?;
    m.add_class::<ConformationalIndex>()?;
    Ok(())
}

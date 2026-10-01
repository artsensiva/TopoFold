//! High-throughput molecular dynamics trajectory streaming I/O (DCD, PDB) for TopoFold.

#![forbid(unsafe_code)]

pub mod dcd;
pub mod error;
pub mod pdb_traj;

pub use dcd::{read_dcd_trajectory, DcdFrame, DcdHeader, DcdReader, DcdWriter, UnitCell};
pub use error::TrajectoryError;
pub use pdb_traj::read_pdb_trajectory;

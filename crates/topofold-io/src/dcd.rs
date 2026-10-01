//! High-performance streaming parser and writer for DCD trajectory files.
//!
//! DCD is the standard binary trajectory format used by CHARMM, NAMD, and OpenMM.
//! This implementation is pure Rust, zero-copy, and supports streaming without
//! loading the full trajectory into memory.

#![forbid(unsafe_code)]

use std::io::{Read, Write};
use nalgebra::Point3;
use topofold_core::types::BackboneTrace;

use crate::error::TrajectoryError;

/// Header metadata for a DCD trajectory file.
#[derive(Debug, Clone, PartialEq)]
pub struct DcdHeader {
    /// Format signature (typically "CORD" for coordinates or "VELD" for velocities).
    pub signature: String,
    /// Total number of frames in the trajectory (0 if undetermined).
    pub n_frames: usize,
    /// Starting timestep index.
    pub start_step: usize,
    /// Timestep interval between frames.
    pub step_interval: usize,
    /// Total number of integration steps.
    pub n_steps: usize,
    /// Integration timestep (in AKMA units or picoseconds).
    pub timestep: f32,
    /// Whether each frame includes a 48-byte unit cell dimension block.
    pub has_unit_cell: bool,
    /// CHARMM format version (e.g. 24 or 31).
    pub charmm_version: i32,
    /// Title lines stored in the trajectory.
    pub titles: Vec<String>,
    /// Total number of atoms per frame.
    pub n_atoms: usize,
}

/// Unit cell dimensions for a periodic simulation box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitCell {
    /// Box edge length A (Å).
    pub a: f64,
    /// Angle gamma (degrees).
    pub gamma: f64,
    /// Box edge length B (Å).
    pub b: f64,
    /// Angle beta (degrees).
    pub beta: f64,
    /// Angle alpha (degrees).
    pub alpha: f64,
    /// Box edge length C (Å).
    pub c: f64,
}

/// Single coordinate frame parsed from a DCD stream.
#[derive(Debug, Clone)]
pub struct DcdFrame {
    /// 0-based frame index within the trajectory.
    pub frame_idx: usize,
    /// Optional periodic box unit cell.
    pub unit_cell: Option<UnitCell>,
    /// Cartesian X coordinates (in Ångströms).
    pub x: Vec<f32>,
    /// Cartesian Y coordinates (in Ångströms).
    pub y: Vec<f32>,
    /// Cartesian Z coordinates (in Ångströms).
    pub z: Vec<f32>,
}

impl DcdFrame {
    /// Converts full-atom coordinates into a `BackboneTrace` using specific atom indices (e.g. C-alpha).
    ///
    /// # Errors
    /// Returns [`TrajectoryError::IndexOutOfBounds`] if any index is `>= n_atoms`.
    pub fn to_backbone_trace(&self, atom_indices: &[usize]) -> Result<BackboneTrace, TrajectoryError> {
        let n_atoms = self.x.len();
        let mut points = Vec::with_capacity(atom_indices.len());

        for &idx in atom_indices {
            if idx >= n_atoms {
                return Err(TrajectoryError::IndexOutOfBounds {
                    index: idx,
                    n_atoms,
                });
            }
            points.push(Point3::new(
                self.x[idx] as f64,
                self.y[idx] as f64,
                self.z[idx] as f64,
            ));
        }

        Ok(BackboneTrace::new(points))
    }

    /// Converts all atom coordinates directly into a `BackboneTrace`.
    #[must_use]
    pub fn to_all_atoms_trace(&self) -> BackboneTrace {
        let n_atoms = self.x.len();
        let mut points = Vec::with_capacity(n_atoms);
        for i in 0..n_atoms {
            points.push(Point3::new(
                self.x[i] as f64,
                self.y[i] as f64,
                self.z[i] as f64,
            ));
        }
        BackboneTrace::new(points)
    }
}

/// Reads a single Fortran unformatted binary record.
///
/// Returns `Ok(None)` if clean EOF is reached at the start of a record.
fn read_fortran_record<R: Read>(reader: &mut R, buf: &mut Vec<u8>) -> Result<Option<usize>, TrajectoryError> {
    let mut len_bytes = [0u8; 4];
    match reader.read_exact(&mut len_bytes) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(TrajectoryError::Io(e)),
    }

    let record_len = u32::from_le_bytes(len_bytes) as usize;
    buf.resize(record_len, 0);
    reader.read_exact(buf)?;

    let mut trailing_len_bytes = [0u8; 4];
    reader.read_exact(&mut trailing_len_bytes)?;
    let trailing_len = u32::from_le_bytes(trailing_len_bytes) as usize;

    if record_len != trailing_len {
        return Err(TrajectoryError::RecordSizeMismatch {
            expected: record_len,
            found: trailing_len,
        });
    }

    Ok(Some(record_len))
}

/// Writes a single Fortran unformatted binary record.
fn write_fortran_record<W: Write>(writer: &mut W, data: &[u8]) -> Result<(), TrajectoryError> {
    let len = data.len() as u32;
    writer.write_all(&len.to_le_bytes())?;
    writer.write_all(data)?;
    writer.write_all(&len.to_le_bytes())?;
    Ok(())
}

/// Streaming reader for DCD trajectory files.
pub struct DcdReader<R> {
    reader: R,
    header: DcdHeader,
    current_frame: usize,
    buf: Vec<u8>,
}

impl<R: Read> DcdReader<R> {
    /// Initializes a DCD reader by parsing the header records from the stream.
    ///
    /// # Errors
    /// Returns [`TrajectoryError`] if the header is malformed.
    pub fn new(mut reader: R) -> Result<Self, TrajectoryError> {
        let mut buf = Vec::new();

        // 1. Read Header Record (84 bytes)
        let rec_len = read_fortran_record(&mut reader, &mut buf)?
            .ok_or(TrajectoryError::UnexpectedEof)?;

        if rec_len != 84 {
            return Err(TrajectoryError::InvalidHeader(format!(
                "Expected 84 bytes in DCD header record, found {rec_len}"
            )));
        }

        let signature = String::from_utf8_lossy(&buf[0..4]).to_string();
        if signature != "CORD" && signature != "VELD" {
            return Err(TrajectoryError::InvalidHeader(format!(
                "Invalid DCD signature: '{signature}', expected 'CORD' or 'VELD'"
            )));
        }

        let n_frames = i32::from_le_bytes(buf[4..8].try_into().unwrap()) as usize;
        let start_step = i32::from_le_bytes(buf[8..12].try_into().unwrap()) as usize;
        let step_interval = i32::from_le_bytes(buf[12..16].try_into().unwrap()) as usize;
        let n_steps = i32::from_le_bytes(buf[16..20].try_into().unwrap()) as usize;
        let unit_cell_flag = i32::from_le_bytes(buf[36..40].try_into().unwrap());
        let has_unit_cell = unit_cell_flag != 0;
        let timestep = f32::from_le_bytes(buf[40..44].try_into().unwrap());
        let charmm_version = i32::from_le_bytes(buf[80..84].try_into().unwrap());

        // 2. Read Title Record
        let title_len = read_fortran_record(&mut reader, &mut buf)?
            .ok_or(TrajectoryError::UnexpectedEof)?;

        if title_len < 4 {
            return Err(TrajectoryError::InvalidHeader(
                "Title block is truncated".to_string(),
            ));
        }

        let n_titles = i32::from_le_bytes(buf[0..4].try_into().unwrap()) as usize;
        let mut titles = Vec::with_capacity(n_titles);
        let mut offset = 4;
        for _ in 0..n_titles {
            if offset + 80 <= title_len {
                let title_line = String::from_utf8_lossy(&buf[offset..offset + 80])
                    .trim()
                    .to_string();
                titles.push(title_line);
                offset += 80;
            }
        }

        // 3. Read Atom Count Record (4 bytes)
        let atoms_len = read_fortran_record(&mut reader, &mut buf)?
            .ok_or(TrajectoryError::UnexpectedEof)?;

        if atoms_len != 4 {
            return Err(TrajectoryError::InvalidHeader(format!(
                "Expected 4 bytes in atom count record, found {atoms_len}"
            )));
        }

        let n_atoms = i32::from_le_bytes(buf[0..4].try_into().unwrap()) as usize;

        let header = DcdHeader {
            signature,
            n_frames,
            start_step,
            step_interval,
            n_steps,
            timestep,
            has_unit_cell,
            charmm_version,
            titles,
            n_atoms,
        };

        Ok(Self {
            reader,
            header,
            current_frame: 0,
            buf,
        })
    }

    /// Returns a reference to the parsed DCD header metadata.
    #[must_use]
    pub fn header(&self) -> &DcdHeader {
        &self.header
    }

    /// Reads the next frame from the trajectory stream.
    ///
    /// Returns `Ok(Some(frame))` if a frame was read, or `Ok(None)` if the end of the
    /// trajectory was reached cleanly.
    ///
    /// # Errors
    /// Returns [`TrajectoryError`] if the stream is truncated or record sizes do not match.
    pub fn next_frame(&mut self) -> Result<Option<DcdFrame>, TrajectoryError> {
        let n_atoms = self.header.n_atoms;
        let expected_coord_bytes = n_atoms * 4;

        // 1. Read optional unit cell block
        let mut unit_cell = None;
        if self.header.has_unit_cell {
            match read_fortran_record(&mut self.reader, &mut self.buf)? {
                Some(len) => {
                    if len != 48 {
                        return Err(TrajectoryError::RecordSizeMismatch {
                            expected: 48,
                            found: len,
                        });
                    }
                    unit_cell = Some(UnitCell {
                        a: f64::from_le_bytes(self.buf[0..8].try_into().unwrap()),
                        gamma: f64::from_le_bytes(self.buf[8..16].try_into().unwrap()),
                        b: f64::from_le_bytes(self.buf[16..24].try_into().unwrap()),
                        beta: f64::from_le_bytes(self.buf[24..32].try_into().unwrap()),
                        alpha: f64::from_le_bytes(self.buf[32..40].try_into().unwrap()),
                        c: f64::from_le_bytes(self.buf[40..48].try_into().unwrap()),
                    });
                }
                None => return Ok(None), // Clean EOF before frame
            }
        }

        // 2. Read X coordinates
        let x_len = match read_fortran_record(&mut self.reader, &mut self.buf)? {
            Some(len) => len,
            None => {
                if unit_cell.is_some() {
                    return Err(TrajectoryError::UnexpectedEof);
                }
                return Ok(None); // Clean EOF
            }
        };

        if x_len != expected_coord_bytes {
            return Err(TrajectoryError::RecordSizeMismatch {
                expected: expected_coord_bytes,
                found: x_len,
            });
        }

        let mut x = Vec::with_capacity(n_atoms);
        for &chunk in self.buf.as_chunks::<4>().0 {
            x.push(f32::from_le_bytes(chunk));
        }

        // 3. Read Y coordinates
        let y_len = read_fortran_record(&mut self.reader, &mut self.buf)?
            .ok_or(TrajectoryError::UnexpectedEof)?;

        if y_len != expected_coord_bytes {
            return Err(TrajectoryError::RecordSizeMismatch {
                expected: expected_coord_bytes,
                found: y_len,
            });
        }

        let mut y = Vec::with_capacity(n_atoms);
        for &chunk in self.buf.as_chunks::<4>().0 {
            y.push(f32::from_le_bytes(chunk));
        }

        // 4. Read Z coordinates
        let z_len = read_fortran_record(&mut self.reader, &mut self.buf)?
            .ok_or(TrajectoryError::UnexpectedEof)?;

        if z_len != expected_coord_bytes {
            return Err(TrajectoryError::RecordSizeMismatch {
                expected: expected_coord_bytes,
                found: z_len,
            });
        }

        let mut z = Vec::with_capacity(n_atoms);
        for &chunk in self.buf.as_chunks::<4>().0 {
            z.push(f32::from_le_bytes(chunk));
        }

        let frame = DcdFrame {
            frame_idx: self.current_frame,
            unit_cell,
            x,
            y,
            z,
        };

        self.current_frame += 1;
        Ok(Some(frame))
    }
}

/// Streaming writer for creating DCD trajectory files.
pub struct DcdWriter<W> {
    writer: W,
    header: DcdHeader,
    frames_written: usize,
}

impl<W: Write> DcdWriter<W> {
    /// Creates a new `DcdWriter` and writes the header records to the destination stream.
    ///
    /// # Errors
    /// Returns [`TrajectoryError`] if writing fails.
    pub fn new(mut writer: W, header: DcdHeader) -> Result<Self, TrajectoryError> {
        // 1. Write Header Record (84 bytes)
        let mut rec1 = vec![0u8; 84];
        let sig = header.signature.as_bytes();
        let sig_len = sig.len().min(4);
        rec1[0..sig_len].copy_from_slice(&sig[0..sig_len]);

        rec1[4..8].copy_from_slice(&(header.n_frames as i32).to_le_bytes());
        rec1[8..12].copy_from_slice(&(header.start_step as i32).to_le_bytes());
        rec1[12..16].copy_from_slice(&(header.step_interval as i32).to_le_bytes());
        rec1[16..20].copy_from_slice(&(header.n_steps as i32).to_le_bytes());
        let unit_cell_flag: i32 = if header.has_unit_cell { 1 } else { 0 };
        rec1[36..40].copy_from_slice(&unit_cell_flag.to_le_bytes());
        rec1[40..44].copy_from_slice(&header.timestep.to_le_bytes());
        rec1[80..84].copy_from_slice(&header.charmm_version.to_le_bytes());
        write_fortran_record(&mut writer, &rec1)?;

        // 2. Write Title Record
        let n_titles = header.titles.len();
        let mut rec2 = Vec::with_capacity(4 + n_titles * 80);
        rec2.extend_from_slice(&(n_titles as i32).to_le_bytes());
        for title in &header.titles {
            let mut line = [b' '; 80];
            let bytes = title.as_bytes();
            let copy_len = bytes.len().min(80);
            line[0..copy_len].copy_from_slice(&bytes[0..copy_len]);
            rec2.extend_from_slice(&line);
        }
        write_fortran_record(&mut writer, &rec2)?;

        // 3. Write Atom Count Record
        let rec3 = (header.n_atoms as i32).to_le_bytes();
        write_fortran_record(&mut writer, &rec3)?;

        Ok(Self {
            writer,
            header,
            frames_written: 0,
        })
    }

    /// Writes a single frame coordinates into the DCD trajectory stream.
    ///
    /// # Errors
    /// Returns [`TrajectoryError::AtomCountMismatch`] if coordinates length does not match `n_atoms`.
    pub fn write_frame(&mut self, frame: &DcdFrame) -> Result<(), TrajectoryError> {
        let n_atoms = self.header.n_atoms;
        if frame.x.len() != n_atoms || frame.y.len() != n_atoms || frame.z.len() != n_atoms {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: n_atoms,
                found: frame.x.len(),
            });
        }

        // 1. Write optional unit cell
        if self.header.has_unit_cell {
            let mut cell_buf = [0u8; 48];
            if let Some(cell) = frame.unit_cell {
                cell_buf[0..8].copy_from_slice(&cell.a.to_le_bytes());
                cell_buf[8..16].copy_from_slice(&cell.gamma.to_le_bytes());
                cell_buf[16..24].copy_from_slice(&cell.b.to_le_bytes());
                cell_buf[24..32].copy_from_slice(&cell.beta.to_le_bytes());
                cell_buf[32..40].copy_from_slice(&cell.alpha.to_le_bytes());
                cell_buf[40..48].copy_from_slice(&cell.c.to_le_bytes());
            }
            write_fortran_record(&mut self.writer, &cell_buf)?;
        }

        // 2. Write X coordinates
        let mut x_bytes = Vec::with_capacity(n_atoms * 4);
        for &val in &frame.x {
            x_bytes.extend_from_slice(&val.to_le_bytes());
        }
        write_fortran_record(&mut self.writer, &x_bytes)?;

        // 3. Write Y coordinates
        let mut y_bytes = Vec::with_capacity(n_atoms * 4);
        for &val in &frame.y {
            y_bytes.extend_from_slice(&val.to_le_bytes());
        }
        write_fortran_record(&mut self.writer, &y_bytes)?;

        // 4. Write Z coordinates
        let mut z_bytes = Vec::with_capacity(n_atoms * 4);
        for &val in &frame.z {
            z_bytes.extend_from_slice(&val.to_le_bytes());
        }
        write_fortran_record(&mut self.writer, &z_bytes)?;

        self.frames_written += 1;
        Ok(())
    }

    /// Number of frames written so far.
    #[must_use]
    pub fn frames_written(&self) -> usize {
        self.frames_written
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_dcd_round_trip() {
        let n_atoms = 10;
        let n_frames = 5;

        let header = DcdHeader {
            signature: "CORD".to_string(),
            n_frames,
            start_step: 0,
            step_interval: 100,
            n_steps: n_frames * 100,
            timestep: 0.0416,
            has_unit_cell: true,
            charmm_version: 31,
            titles: vec!["TopoFold Unit Test Trajectory".to_string()],
            n_atoms,
        };

        let mut buffer = Vec::new();
        {
            let mut writer = DcdWriter::new(&mut buffer, header.clone()).expect("Init writer");
            for f in 0..n_frames {
                let cell = UnitCell {
                    a: 50.0 + f as f64,
                    gamma: 90.0,
                    b: 50.0,
                    beta: 90.0,
                    alpha: 90.0,
                    c: 50.0,
                };
                let x: Vec<f32> = (0..n_atoms).map(|i| (i as f32) + f as f32).collect();
                let y: Vec<f32> = (0..n_atoms).map(|i| (i as f32) * 2.0).collect();
                let z: Vec<f32> = (0..n_atoms).map(|i| (i as f32) * 3.0).collect();

                let frame = DcdFrame {
                    frame_idx: f,
                    unit_cell: Some(cell),
                    x,
                    y,
                    z,
                };
                writer.write_frame(&frame).expect("Write frame");
            }
            assert_eq!(writer.frames_written(), n_frames);
        }

        // Now read back using DcdReader
        let cursor = Cursor::new(buffer);
        let mut reader = DcdReader::new(cursor).expect("Init reader");
        assert_eq!(reader.header().n_atoms, n_atoms);
        assert_eq!(reader.header().n_frames, n_frames);
        assert!(reader.header().has_unit_cell);

        for f in 0..n_frames {
            let frame = reader.next_frame().expect("Read frame").expect("Frame exists");
            assert_eq!(frame.frame_idx, f);
            assert_eq!(frame.x.len(), n_atoms);
            assert_eq!(frame.x[0], f as f32);
            assert_eq!(frame.x[1], 1.0 + f as f32);

            let cell = frame.unit_cell.expect("Unit cell");
            assert!((cell.a - (50.0 + f as f64)).abs() < 1e-6);

            // Test C-alpha backbone extraction
            let ca_indices = vec![0, 2, 4];
            let trace = frame.to_backbone_trace(&ca_indices).expect("Extract trace");
            assert_eq!(trace.len(), 3);
            assert_eq!(trace.coordinates()[0], Point3::new(f as f64, 0.0, 0.0));
            assert_eq!(trace.coordinates()[1], Point3::new(2.0 + f as f64, 4.0, 6.0));
        }

        // EOF check
        let eof_frame = reader.next_frame().expect("Read at EOF");
        assert!(eof_frame.is_none());
    }
}

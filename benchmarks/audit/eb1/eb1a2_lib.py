"""EB-1 Amendment 2 helpers (docs/EB1_AMENDMENT_2.md). Harness code, not production.

- `read_rows`: read contiguous row ranges of a .npy file straight into one preallocated owned array (no memory map,
  no intermediate full array); raises on any impossible range instead of truncating.
- `read_d2_stored`: first n stored D2 frames across the four trajectory files, read directly; raises if n exceeds the
  available frames.
- `measure_reps`: repeat a measured call, releasing each output before the next repetition.
"""
import gc
import time

import numpy as np


class SizeError(ValueError):
    pass


def _header(f):
    v = np.lib.format.read_magic(f)
    rd = np.lib.format.read_array_header_1_0 if v == (1, 0) else np.lib.format.read_array_header_2_0
    shape, fortran, dtype = rd(f)
    if fortran:
        raise ValueError('Fortran-ordered arrays are not supported')
    return shape, dtype, f.tell()


def npy_shape(path):
    with open(path, 'rb') as f:
        return _header(f)[:2]


def read_rows(path, ranges, out=None):
    """ranges: list of (start, n) contiguous row blocks. Returns an owned C-contiguous array of sum(n) rows, or fills
    `out` (a C-contiguous array of exactly that shape and dtype) in place without any intermediate array."""
    with open(path, 'rb', buffering=0) as f:                    # unbuffered: readinto writes straight into `out`
        shape, dtype, off = _header(f)
        row_bytes = int(np.prod(shape[1:], dtype=np.int64)) * dtype.itemsize
        total = sum(n for _, n in ranges)
        for s, n in ranges:
            if s < 0 or n < 0 or s + n > shape[0]:
                raise SizeError(f'requested rows {s}..{s + n} exceed the {shape[0]} rows of {path.name}')
        if out is None:
            out = np.empty((total,) + tuple(shape[1:]), dtype)
        elif out.shape != (total,) + tuple(shape[1:]) or out.dtype != dtype or not out.flags['C_CONTIGUOUS']:
            raise SizeError('output buffer does not match the requested rows')
        buf = memoryview(out).cast('B')
        pos = 0
        for s, n in ranges:
            f.seek(off + s * row_bytes)
            nb = n * row_bytes
            if f.readinto(buf[pos:pos + nb]) != nb:
                raise SizeError(f'short read from {path.name}')
            pos += nb
    return out


def read_d2_stored(cache, n):
    """First n stored (pre-burn-in) D2 frames across trajectories 0–3, with per-frame boxes."""
    lens = [npy_shape(cache / f'D2_stored_{i}.npy')[0][0] for i in range(4)]
    if n > sum(lens):
        raise SizeError(f'requested {n} D2 stored frames; {sum(lens)} available')
    L = npy_shape(cache / 'D2_stored_0.npy')[0][1]
    st = np.empty((n, L, 3), np.float32)
    bx = np.empty((n, 3))
    pos = 0
    for i, ln in enumerate(lens):
        take = min(ln, n - pos)
        if take <= 0:
            break
        read_rows(cache / f'D2_stored_{i}.npy', [(0, take)], out=st[pos:pos + take])
        bx[pos:pos + take] = np.load(cache / f'D2_box_{i}.npy')
        pos += take
    return st, bx


def validate_size(requested, available, what):
    if requested > available:
        raise SizeError(f'{what}: requested {requested} > available {available}')
    return requested


def measure_reps(fn, max_reps=3, long_s=5.0):
    """Run fn() up to max_reps times (once if the first run takes ≥ long_s); each output is released (and garbage
    collected) before the next repetition, so at most one output is alive at a time. Returns (last output, timings)."""
    walls, cpus, last = [], [], None
    for r in range(max_reps):
        if last is not None:
            del last
            last = None
            gc.collect()
        t0, c0 = time.perf_counter(), time.process_time()
        last = fn()
        walls.append(time.perf_counter() - t0)
        cpus.append(time.process_time() - c0)
        if r == 0 and walls[0] >= long_s:
            break
    return last, walls, cpus

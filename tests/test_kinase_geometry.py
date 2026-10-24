"""Pure geometry tests for topofold_inputs.kinase_geometry (no MDAnalysis needed): IUPAC dihedral reference cases, an
independent formula on SYNTHETIC points, and NaN (never zero) for degenerate or missing input."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.kinase_geometry import dihedral, distance  # noqa: E402

def test_dihedral_reference_cases_and_independent_formula():
    P = lambda *x: np.array(x, float)  # noqa: E731
    assert dihedral(P(1, 0, 0), P(0, 0, 0), P(0, 1, 0), P(-1, 1, 0))[0] == pytest.approx(math.pi)   # trans
    assert dihedral(P(1, 0, 0), P(0, 0, 0), P(0, 1, 0), P(1, 1, 0))[0] == pytest.approx(0.0)        # cis
    # IUPAC: front bond at 3 o'clock, back bond at 12 o'clock looking from p1 to p2 -> -90 degrees
    assert dihedral(P(1, 0, 0), P(0, 0, 0), P(0, 1, 0), P(0, 1, 1))[0] == pytest.approx(-math.pi / 2)
    rng = np.random.default_rng(7)
    pts = rng.normal(size=(2000, 4, 3))

    def ref(p0, p1, p2, p3):          # independent textbook form, opposite sign convention
        b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        return math.atan2(np.dot(m1, n2), np.dot(n1, n2))

    ours = dihedral(pts[:, 0], pts[:, 1], pts[:, 2], pts[:, 3])
    theirs = -np.array([ref(*q) for q in pts])
    d = np.abs(np.angle(np.exp(1j * (ours - theirs))))
    assert d.max() < 1e-12
    assert np.all((ours > -math.pi) & (ours <= math.pi))


def test_dihedral_degenerate_and_missing_give_nan_not_zero():
    a = np.array([[0.0, 0, 0]])
    assert np.isnan(dihedral(a, np.array([[1.0, 0, 0]]), np.array([[2.0, 0, 0]]), np.array([[3.0, 0, 0]]))[0])
    nan = np.full((1, 3), np.nan)
    assert np.isnan(dihedral(nan, a, a + 1, a + 2)[0])
    assert np.isnan(distance(nan, a)[0])
    assert distance(np.array([[0.0, 0, 0]]), np.array([[3.0, 4, 0]]))[0] == pytest.approx(5.0)



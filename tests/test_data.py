"""Load the real Pantheon+ files when they are present."""

from pathlib import Path

import numpy as np
import pytest

from likelihood import load_pantheon_plus

ROOT = Path(__file__).resolve().parents[1]
DAT = ROOT / "data" / "Pantheon+SH0ES.dat"
COV = ROOT / "data" / "Pantheon+SH0ES_STAT+SYS.cov"

pytestmark = pytest.mark.skipif(
    not (DAT.exists() and COV.exists()),
    reason="Pantheon+ data not downloaded",
)


def test_pantheon_plus_loads_and_applies_z_cut():
    data = load_pantheon_plus(DAT, COV, zmin=0.01)
    assert data.n_raw == 1701
    assert data.n < data.n_raw
    assert data.n > 1000
    assert data.z.min() > 0.01
    assert data.cov.shape == (data.n, data.n)
    assert np.all(np.isfinite(data.mu))
    assert np.all(data.mu > 0.0)
    # Covariance must be usable as a χ² metric.
    eig = np.linalg.eigvalsh(data.cov[:32, :32])
    assert np.all(eig > 0.0)

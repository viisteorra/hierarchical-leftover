"""DES-SN5YR loader and DESI DR2 sample, when the files are present."""

from pathlib import Path

import numpy as np
import pytest

from bao import DESI_DR2, n_bao
from likelihood import load_des_dovekie

ROOT = Path(__file__).resolve().parents[1]
HD = ROOT / "data" / "DES-Dovekie_HD.csv"
NPZ = ROOT / "data" / "DES-Dovekie_STAT+SYS.npz"


def test_desi_dr2_counts_nonoverlapping_ratios():
    # BGS DV + 6 anisotropic tracers.
    assert n_bao(DESI_DR2) == 13
    names = [p.name for p in DESI_DR2]
    assert names.count("LRG3+ELG1") == 1
    assert "LRG3" not in names


@pytest.mark.skipif(not (HD.exists() and NPZ.exists()), reason="DES-SN5YR not downloaded")
def test_des_dovekie_loads():
    data = load_des_dovekie(HD, NPZ)
    assert data.n == 1829 or data.n == 1820 or data.n > 1500
    assert data.cov.shape == (data.n, data.n)
    assert np.all(np.isfinite(data.mu))
    eig = np.linalg.eigvalsh(data.cov[:24, :24])
    assert np.all(eig > 0.0)

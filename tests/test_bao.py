"""DESI BAO distances recover a known rd; locked Ω_DE is a shape, not a scale."""

import numpy as np

from bao import BAOPoint, bao_chi2, fit_rd, n_bao
from cosmology import D_H, D_M, D_V
from geometry import OMEGA_DE_TODAY


def test_dv_identity():
    z, H0 = 0.7, 70.0
    dv = D_V(z, H0)
    expect = (z * D_M(z, H0) ** 2 * D_H(z, H0)) ** (1.0 / 3.0)
    assert abs(dv - expect) < 1e-9


def test_n_bao_counts_ratios_not_tracers():
    # 2 DV + 5 DM_DH pairs = 12 numbers (BGS, QSO + 5 anisotropic).
    assert n_bao() == 12


def test_perfect_theory_has_zero_chi2():
    H0, rd = 70.0, 147.0
    sample = (
        BAOPoint("DV", 0.30, "DV", (D_V(0.30, H0) / rd,), (0.15,)),
        BAOPoint(
            "aniso",
            0.71,
            "DM_DH",
            (D_M(0.71, H0) / rd, D_H(0.71, H0) / rd),
            (0.3, 0.6),
            -0.4,
        ),
    )
    assert bao_chi2(H0, rd, sample=sample) < 1e-10


def test_fit_rd_recovers_truth():
    H0, true_rd = 70.0, 147.09
    sample = (
        BAOPoint("BGS", 0.30, "DV", (D_V(0.30, H0, omega_de=OMEGA_DE_TODAY) / true_rd,), (0.05,)),
        BAOPoint(
            "LRG",
            0.71,
            "DM_DH",
            (
                D_M(0.71, H0, omega_de=OMEGA_DE_TODAY) / true_rd,
                D_H(0.71, H0, omega_de=OMEGA_DE_TODAY) / true_rd,
            ),
            (0.08, 0.10),
            -0.4,
        ),
        BAOPoint("QSO", 1.49, "DV", (D_V(1.49, H0, omega_de=OMEGA_DE_TODAY) / true_rd,), (0.10,)),
    )
    fit = fit_rd(H0, omega_de=OMEGA_DE_TODAY, sample=sample)
    assert fit["success"]
    assert abs(fit["rd"] - true_rd) < 0.05
    assert fit["chi2"] < 1e-6


def test_h0_rd_degeneracy_on_shape():
    """Rescaling H0 and rd together leaves BAO ratios unchanged at fixed Ω_DE."""
    z = 0.93
    a = D_M(z, 70.0, omega_de=OMEGA_DE_TODAY) / 147.0
    b = D_M(z, 80.0, omega_de=OMEGA_DE_TODAY) / (147.0 * 70.0 / 80.0)
    assert abs(a - b) < 1e-10

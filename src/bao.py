"""DESI DR1 BAO distances. Shape test of the locked Ω_m = 1 − 0.69014.

Numbers are Table 18 of DESI 2024 III (Adame et al., arXiv:2404.03000),
the configuration-space fiducial results used in the cosmology paper.
Lyman-α is DESI 2024 IV (arXiv:2404.03001).

BAO measures D/rd. With Ω_DE locked, E(z) is fixed, so the only free
scale is rd (or H0·rd). We never mix a SH0ES H0 with a Planck rd in a
joint χ² — that is the Hubble tension, not a test of 0.69014.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from cosmology import D_H, D_M, D_V
from geometry import OMEGA_DE_TODAY

# Planck 2018 TT,TE,EE+lowE+lensing sound horizon (Mpc). Used only for
# the inverse-distance-ladder H0 readout, not for the shape χ².
PLANCK_RD_MPC = 147.09


@dataclass(frozen=True)
class BAOPoint:
    name: str
    z: float
    kind: str  # "DV" or "DM_DH"
    value: tuple[float, ...]
    sigma: tuple[float, ...]
    corr: float = 0.0


# Cosmology sample: non-overlapping tracers as in DESI 2024 VI.
DESI_DR1 = (
    BAOPoint("BGS", 0.30, "DV", (7.93,), (0.15,)),
    BAOPoint("LRG1", 0.51, "DM_DH", (13.62, 20.98), (0.25, 0.61), -0.445),
    BAOPoint("LRG2", 0.71, "DM_DH", (16.85, 20.08), (0.32, 0.60), -0.420),
    BAOPoint("LRG3+ELG1", 0.93, "DM_DH", (21.71, 17.88), (0.28, 0.35), -0.389),
    BAOPoint("ELG2", 1.32, "DM_DH", (27.79, 13.82), (0.69, 0.42), -0.444),
    BAOPoint("QSO", 1.49, "DV", (26.07,), (0.67,)),
    BAOPoint("Lyα", 2.33, "DM_DH", (39.71, 8.52), (0.94, 0.17), -0.477),
)

# DESI DR2 Table 4 of Abdul-Karim et al. 2025, arXiv:2503.14738.
# Non-overlapping cosmology set (LRG3+ELG1, not LRG3 and ELG1 separately).
DESI_DR2 = (
    BAOPoint("BGS", 0.295, "DV", (7.942,), (0.075,)),
    BAOPoint("LRG1", 0.510, "DM_DH", (13.588, 21.863), (0.167, 0.425), -0.459),
    BAOPoint("LRG2", 0.706, "DM_DH", (17.351, 19.455), (0.177, 0.330), -0.404),
    BAOPoint("LRG3+ELG1", 0.934, "DM_DH", (21.576, 17.641), (0.152, 0.193), -0.416),
    BAOPoint("ELG2", 1.321, "DM_DH", (27.601, 14.176), (0.318, 0.221), -0.434),
    BAOPoint("QSO", 1.484, "DM_DH", (30.512, 12.817), (0.760, 0.516), -0.500),
    BAOPoint("Lyα", 2.330, "DM_DH", (38.988, 8.632), (0.531, 0.101), -0.431),
)


def _theory(point: BAOPoint, H0: float, rd: float, omega_de: float, w: float) -> np.ndarray:
    if point.kind == "DV":
        return np.array([D_V(point.z, H0, w=w, omega_de=omega_de) / rd], dtype=float)
    if point.kind == "DM_DH":
        return np.array(
            [
                D_M(point.z, H0, w=w, omega_de=omega_de) / rd,
                D_H(point.z, H0, w=w, omega_de=omega_de) / rd,
            ],
            dtype=float,
        )
    raise ValueError(point.kind)


def _cov(point: BAOPoint) -> np.ndarray:
    sig = np.array(point.sigma, dtype=float)
    if sig.size == 1:
        return np.array([[sig[0] ** 2]])
    cov = np.diag(sig**2)
    cov[0, 1] = cov[1, 0] = point.corr * sig[0] * sig[1]
    return cov


def bao_chi2(
    H0: float,
    rd: float,
    omega_de: float = OMEGA_DE_TODAY,
    w: float = -1.0,
    sample=DESI_DR1,
) -> float:
    total = 0.0
    for point in sample:
        delta = np.array(point.value, dtype=float) - _theory(point, H0, rd, omega_de, w)
        cov = _cov(point)
        total += float(delta @ np.linalg.solve(cov, delta))
    return total


def n_bao(sample=DESI_DR1) -> int:
    return int(sum(len(p.value) for p in sample))


def fit_rd(
    H0: float,
    omega_de: float = OMEGA_DE_TODAY,
    w: float = -1.0,
    sample=DESI_DR1,
    rd_bounds: tuple[float, float] = (100.0, 200.0),
) -> dict:
    """Shape-only: Ω_DE and H0 fixed, fit rd. H0 here is a dummy scale."""

    def objective(rd: float) -> float:
        return bao_chi2(H0, float(rd), omega_de=omega_de, w=w, sample=sample)

    result = minimize_scalar(objective, bounds=rd_bounds, method="bounded")
    rd = float(result.x)
    chi2 = float(result.fun)
    n = n_bao(sample)
    return {
        "H0": H0,
        "rd": rd,
        "rd_h": rd * H0 / 100.0,
        "omega_de": omega_de,
        "w": w,
        "chi2": chi2,
        "n": n,
        "npar": 1,
        "dof": n - 1,
        "chi2_dof": chi2 / (n - 1),
        "success": bool(result.success),
    }


def fit_h0_given_rd(
    rd: float = PLANCK_RD_MPC,
    omega_de: float = OMEGA_DE_TODAY,
    w: float = -1.0,
    sample=DESI_DR1,
    H0_bounds: tuple[float, float] = (50.0, 90.0),
) -> dict:
    """Inverse distance ladder: Planck rd, locked Ω_DE, fit H0 from BAO."""

    def objective(H0: float) -> float:
        return bao_chi2(float(H0), rd, omega_de=omega_de, w=w, sample=sample)

    result = minimize_scalar(objective, bounds=H0_bounds, method="bounded")
    H0 = float(result.x)
    chi2 = float(result.fun)
    n = n_bao(sample)
    return {
        "H0": H0,
        "rd": rd,
        "omega_de": omega_de,
        "w": w,
        "chi2": chi2,
        "n": n,
        "npar": 1,
        "dof": n - 1,
        "chi2_dof": chi2 / (n - 1),
        "success": bool(result.success),
        "note": "Planck rd + BAO; this is not the SH0ES H0",
    }


def fit_rd_omega_de(
    H0: float = 70.0,
    w: float = -1.0,
    sample=DESI_DR1,
    omega_bounds: tuple[float, float] = (0.20, 0.85),
) -> dict:
    """COMPARISON ONLY. Float Ω_DE against BAO shape; do not replace the lock."""

    def objective(omega_de: float) -> float:
        return fit_rd(H0, omega_de=float(omega_de), w=w, sample=sample)["chi2"]

    result = minimize_scalar(objective, bounds=omega_bounds, method="bounded")
    omega_de = float(result.x)
    inner = fit_rd(H0, omega_de=omega_de, w=w, sample=sample)
    inner["npar"] = 2
    inner["dof"] = n_bao(sample) - 2
    inner["chi2_dof"] = inner["chi2"] / inner["dof"]
    inner["comparison_only"] = True
    return inner

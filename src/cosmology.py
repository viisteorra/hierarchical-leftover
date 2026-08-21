"""Flat wCDM distances with Ω_DE locked by geometry unless explicitly overridden.

Radiation is omitted for SN-scale distances (z ≲ 2). It is included only in
the CMB shift-parameter helper, where it actually matters. Curvature is zero.
Present-day Ω_DE is read from geometry.py and is not a fit parameter.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

from geometry import OMEGA_DE_TODAY

# km/s. d_L then comes out in Mpc when H0 is in km/s/Mpc.
C_KMS = 299792.458

# μ = 5 log10(d_L / 10 pc) with d_L in Mpc is 5 log10(d_L) + 25.
MU_OFFSET = 25.0

# 1/H0 in Gyr when H0 is in km s^-1 Mpc^-1.
GYR_PER_HO_UNIT = 977.7922216804897

# Radiation density for CMB integrals only (T_CMB = 2.7255 K, N_eff = 3.046).
T_CMB = 2.7255
N_EFF = 3.046
_OMEGA_GAMMA_H2 = 2.469e-5 * (T_CMB / 2.725) ** 4
OMEGA_R_H2 = _OMEGA_GAMMA_H2 * (1.0 + N_EFF * (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0))

_GRID_POINTS_LOWZ = 10001
_GRID_POINTS_HIGHZ = 40001


def hubble_E(z, omega_de: float = OMEGA_DE_TODAY, w: float = -1.0, omega_r: float = 0.0):
    """E(z) = H(z)/H0 for a flat mix of dust, dark energy, and optional radiation."""
    z = np.asarray(z, dtype=float)
    omega_m = 1.0 - omega_de - omega_r
    de_exp = 3.0 * (1.0 + w)
    return np.sqrt(
        omega_m * (1.0 + z) ** 3
        + omega_r * (1.0 + z) ** 4
        + omega_de * (1.0 + z) ** de_exp
    )


def comoving_chi(
    z,
    omega_de: float = OMEGA_DE_TODAY,
    w: float = -1.0,
    omega_r: float = 0.0,
):
    """Dimensionless comoving distance χ(z) = ∫_0^z dz'/E(z')."""
    z = np.asarray(z, dtype=float)
    z_flat = np.atleast_1d(z)
    z_max = float(np.max(z_flat))
    if z_max <= 0.0:
        chi = np.zeros_like(z_flat, dtype=float)
        return chi.reshape(z.shape) if z.shape else chi[0]

    ngrid = _GRID_POINTS_HIGHZ if z_max > 10.0 else _GRID_POINTS_LOWZ
    zg = np.linspace(0.0, z_max, ngrid)
    inv_E = 1.0 / hubble_E(zg, omega_de=omega_de, w=w, omega_r=omega_r)
    chi_grid = cumulative_trapezoid(inv_E, zg, initial=0.0)
    chi = np.interp(z_flat, zg, chi_grid)
    return chi.reshape(z.shape) if z.shape else float(chi[0])


def luminosity_distance(
    z,
    H0: float,
    w: float = -1.0,
    omega_de: float = OMEGA_DE_TODAY,
    omega_r: float = 0.0,
):
    """Luminosity distance in Mpc for a flat universe.

    d_L(z) = (1+z) (c/H0) ∫_0^z dz'/E(z')
    """
    z = np.asarray(z, dtype=float)
    chi = comoving_chi(z, omega_de=omega_de, w=w, omega_r=omega_r)
    return (1.0 + z) * (C_KMS / H0) * chi


def distance_modulus(
    z,
    H0: float,
    w: float = -1.0,
    omega_de: float = OMEGA_DE_TODAY,
    omega_r: float = 0.0,
):
    """Distance modulus μ(z) = 5 log10(d_L/Mpc) + 25."""
    d_l = luminosity_distance(z, H0, w=w, omega_de=omega_de, omega_r=omega_r)
    d_l = np.asarray(d_l, dtype=float)
    if np.any(d_l <= 0.0):
        raise ValueError("luminosity distance must be positive (z > 0)")
    return 5.0 * np.log10(d_l) + MU_OFFSET


def D_M(z, H0: float, w: float = -1.0, omega_de: float = OMEGA_DE_TODAY):
    """Comoving transverse distance D_M in Mpc (flat: D_M = χ)."""
    return (C_KMS / H0) * comoving_chi(z, omega_de=omega_de, w=w)


def D_H(z, H0: float, w: float = -1.0, omega_de: float = OMEGA_DE_TODAY):
    """Hubble distance c/H(z) in Mpc."""
    return (C_KMS / H0) / hubble_E(z, omega_de=omega_de, w=w)


def D_V(z, H0: float, w: float = -1.0, omega_de: float = OMEGA_DE_TODAY):
    """Spherically averaged BAO distance [z D_M² D_H]^(1/3) in Mpc."""
    z = np.asarray(z, dtype=float)
    return (z * D_M(z, H0, w=w, omega_de=omega_de) ** 2 * D_H(z, H0, w=w, omega_de=omega_de)) ** (
        1.0 / 3.0
    )


def age_gyr(
    H0: float,
    w: float = -1.0,
    omega_de: float = OMEGA_DE_TODAY,
    omega_r: float = 0.0,
) -> float:
    """Age of the universe today in Gyr."""
    a = np.geomspace(1e-8, 1.0, 20001)
    z = 1.0 / a - 1.0
    ez = hubble_E(z, omega_de=omega_de, w=w, omega_r=omega_r)
    integral = trapezoid(1.0 / (a * ez), a)
    return float(GYR_PER_HO_UNIT / H0 * integral)


def z_star_hu_sugiyama(omega_b_h2: float, omega_m_h2: float) -> float:
    """Photon-decoupling redshift, Hu & Sugiyama (1996)."""
    g1 = 0.0783 * omega_b_h2 ** (-0.238) / (1.0 + 39.5 * omega_b_h2**0.763)
    g2 = 0.560 / (1.0 + 21.1 * omega_b_h2**1.81)
    return 1048.0 * (1.0 + 0.00124 * omega_b_h2 ** (-0.738)) * (1.0 + g1 * omega_m_h2**g2)


def shift_parameter_R(
    H0: float,
    omega_de: float = OMEGA_DE_TODAY,
    w: float = -1.0,
    omega_b_h2: float = 0.02236,
) -> float:
    """CMB shift parameter R = sqrt(Ω_m) ∫_0^{z*} dz/E(z), radiation included.

    R is nearly independent of H0 except through Ω_r and z*(ω_m). Use a
    Planck-like H0 when comparing to the Planck R prior.
    """
    h = H0 / 100.0
    omega_r = OMEGA_R_H2 / h**2
    omega_m = 1.0 - omega_de - omega_r
    omega_m_h2 = omega_m * h**2
    z_star = z_star_hu_sugiyama(omega_b_h2, omega_m_h2)
    chi = comoving_chi(z_star, omega_de=omega_de, w=w, omega_r=omega_r)
    return float(np.sqrt(omega_m) * chi)

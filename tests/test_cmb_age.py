"""CMB shift parameter and ages: cross-checks against astropy / Planck-ish R."""

import numpy as np
from astropy.cosmology import FlatwCDM

from cosmology import age_gyr, shift_parameter_R, z_star_hu_sugiyama
from geometry import OMEGA_DE_TODAY, OMEGA_M
from likelihood import PLANCK2018_H0, PLANCK2018_OMEGA_DE, PLANCK_R


def test_age_matches_astropy_no_radiation():
    H0 = 67.36
    got = age_gyr(H0, w=-1.0, omega_de=OMEGA_DE_TODAY, omega_r=0.0)
    ref = FlatwCDM(H0=H0, Om0=OMEGA_M, w0=-1.0, Tcmb0=0).age(0).value
    assert abs(got - ref) < 0.01  # 10 Myr


def test_higher_H0_makes_a_younger_universe():
    t_slow = age_gyr(67.0, omega_de=OMEGA_DE_TODAY)
    t_fast = age_gyr(73.5, omega_de=OMEGA_DE_TODAY)
    assert t_fast < t_slow


def test_z_star_is_around_1090():
    z = z_star_hu_sugiyama(0.02236, 0.143)
    assert 1080 < z < 1100


def test_shift_parameter_near_planck_for_planck_omega():
    """R at Planck Ω_DE and Planck H0 must sit on the Chen+2019 prior.

    This is a pipeline sanity check, not a test of the geometric lock.
    """
    R = shift_parameter_R(PLANCK2018_H0, omega_de=PLANCK2018_OMEGA_DE, w=-1.0)
    # Approximate R (no neutrinos in the late-time split, fitting formula z*):
    # should land within ~1% of the published 1.7502.
    assert abs(R - PLANCK_R) / PLANCK_R < 0.015


def test_geometric_R_is_finite_and_nearby():
    R = shift_parameter_R(PLANCK2018_H0, omega_de=OMEGA_DE_TODAY, w=-1.0)
    assert np.isfinite(R)
    assert abs(R - PLANCK_R) < 0.05

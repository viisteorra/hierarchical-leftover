"""Distance-module smoke tests against the Hubble law and astropy."""

import numpy as np
from astropy.cosmology import FlatwCDM

from cosmology import C_KMS, distance_modulus, hubble_E, luminosity_distance
from geometry import OMEGA_DE_TODAY, OMEGA_M


def test_friedmann_equation_today():
    assert abs(hubble_E(0.0) - 1.0) < 1e-12
    assert abs(hubble_E(0.0, omega_de=OMEGA_DE_TODAY, w=-1.0) - 1.0) < 1e-12


def test_low_z_hubble_law():
    H0 = 70.0
    z = np.array([1e-4, 3e-4, 1e-3])
    d_l = luminosity_distance(z, H0, w=-1.0)
    d_hubble = C_KMS * z / H0
    # First order: d_L = cz/H0. Relative error should be O(z).
    rel = np.abs(d_l / d_hubble - 1.0)
    assert np.all(rel < 5.0 * z)


def test_luminosity_distance_zero_at_z_zero():
    assert luminosity_distance(0.0, 70.0) == 0.0


def test_matches_astropy_flat_wcdm():
    H0 = 73.0
    w = -1.0
    z = np.array([0.01, 0.05, 0.1, 0.3, 0.7, 1.0, 1.5, 2.0])
    # Tcmb0=0 drops radiation so the comparison is the same model.
    ref = FlatwCDM(H0=H0, Om0=OMEGA_M, w0=w, Tcmb0=0)
    got = luminosity_distance(z, H0, w=w, omega_de=OMEGA_DE_TODAY)
    want = ref.luminosity_distance(z).value
    np.testing.assert_allclose(got, want, rtol=1e-5, atol=1e-4)

    mu_got = distance_modulus(z, H0, w=w, omega_de=OMEGA_DE_TODAY)
    mu_want = ref.distmod(z).value
    np.testing.assert_allclose(mu_got, mu_want, rtol=1e-5, atol=1e-4)


def test_matches_astropy_when_w_is_not_minus_one():
    H0 = 70.0
    w = -0.9
    z = np.array([0.1, 0.5, 1.0])
    ref = FlatwCDM(H0=H0, Om0=OMEGA_M, w0=w, Tcmb0=0)
    got = luminosity_distance(z, H0, w=w, omega_de=OMEGA_DE_TODAY)
    want = ref.luminosity_distance(z).value
    np.testing.assert_allclose(got, want, rtol=1e-5, atol=1e-4)


def test_scalar_and_array_agree():
    H0 = 70.0
    z = 0.42
    scalar = luminosity_distance(z, H0)
    vector = luminosity_distance(np.array([z]), H0)
    assert abs(float(scalar) - float(vector[0])) < 1e-12


def test_w_minus_one_is_a_cosmological_constant():
    z = np.array([0.2, 0.8, 1.5])
    H0 = 70.0
    a = luminosity_distance(z, H0, w=-1.0)
    b = luminosity_distance(z, H0, w=-1.0 + 1e-12)
    np.testing.assert_allclose(a, b, rtol=1e-9)


def test_de_domination_makes_distances_longer():
    z, H0 = 1.0, 70.0
    high_de = luminosity_distance(z, H0, omega_de=0.8)
    low_de = luminosity_distance(z, H0, omega_de=0.5)
    assert high_de > low_de


def test_hubble_E_positive_on_a_wide_grid():
    z = np.linspace(0.0, 3.0, 50)
    e = hubble_E(z, omega_de=OMEGA_DE_TODAY, w=-1.0)
    assert np.all(e > 0.0)
    assert np.all(np.isfinite(e))

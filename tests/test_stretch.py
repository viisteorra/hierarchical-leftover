"""Stop the cohesion campaign if stretch disagrees with the frozen table."""

import numpy as np

from geometry import OMEGA_DE_TODAY
from hde_stretch_reference import B, OMEGA_DE, OMEGA_DE_FLOAT, REFERENCE_STRETCH
from stretch import stretch_S, stretch_dmu


def test_lock_is_forty_nine_over_seventy_one():
    assert OMEGA_DE.numerator == 49
    assert OMEGA_DE.denominator == 71
    assert abs(OMEGA_DE_FLOAT - OMEGA_DE_TODAY) < 1e-15
    assert B == 12


def test_stretch_matches_reference_table():
    for z, delta0, s_ref, dmu_ref in REFERENCE_STRETCH:
        s = float(np.asarray(stretch_S(z, delta0)))
        dmu = float(np.asarray(stretch_dmu(z, delta0)))
        assert abs(s - s_ref) < 1e-12
        assert abs(dmu - dmu_ref) < 1e-12


def test_zero_stretch_is_identity():
    z = np.array([0.01, 0.5, 1.2])
    np.testing.assert_allclose(stretch_S(z, 0.0), 1.0)
    np.testing.assert_allclose(stretch_dmu(z, 0.0), 0.0)


def test_stretch_fades_toward_one_at_high_z():
    assert stretch_S(0.0, 0.05) > stretch_S(2.0, 0.05)
    assert abs(stretch_S(100.0, 0.05) - 1.0) < 1e-3

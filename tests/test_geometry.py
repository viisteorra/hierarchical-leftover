"""The 1:11 rule and T = 4.9/7.1 are frozen. Do not 'improve' them."""

from geometry import OMEGA_DE_TODAY, OMEGA_M, Q4_WEIGHT, Q5_WEIGHT, TOTAL, W, r


def test_weights_are_the_1_to_11_rule():
    assert Q4_WEIGHT == 1
    assert Q5_WEIGHT == 11
    assert TOTAL == 12


def test_generation_multiplier_matches_weighted_integer_rules():
    expected_r = (1 * (2 / 4) + 11 * (2 / 5)) / 12
    assert r == expected_r
    assert abs(r - 4.9 / 12) < 1e-15


def test_omega_de_is_the_derived_tail():
    assert OMEGA_DE_TODAY == r / (1 - r)
    assert abs(OMEGA_DE_TODAY - 4.9 / 7.1) < 1e-12
    # Five-decimal lock from the design document.
    assert abs(OMEGA_DE_TODAY - 0.69014) < 5e-6


def test_flatness_closes():
    assert OMEGA_M == 1.0 - OMEGA_DE_TODAY
    assert abs(OMEGA_M + OMEGA_DE_TODAY - 1.0) < 1e-15
    assert W == -1.0


def test_not_the_planck_comparison_value():
    # Nearby, but not the same number, and not a fit parameter.
    assert abs(OMEGA_DE_TODAY - 0.689) > 5e-4

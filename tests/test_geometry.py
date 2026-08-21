"""Continuous lock Ω_DE = ln(2). Discrete 12-fold is the rational approximation."""

import math

from geometry import (
    K_CURVATURE,
    LN2,
    OMEGA_DE_TODAY,
    OMEGA_M,
    Q4_WEIGHT,
    Q5_WEIGHT,
    R_DISCRETE,
    T_DISCRETE,
    TOTAL,
    W,
    r,
)


def test_weights_are_the_1_to_11_rule():
    assert Q4_WEIGHT == 1
    assert Q5_WEIGHT == 11
    assert TOTAL == 12


def test_continuous_residual_is_ln2():
    assert abs(LN2 - math.log(2.0)) < 1e-15
    assert abs(OMEGA_DE_TODAY - math.log(2.0)) < 1e-15
    assert abs(r - LN2 / (1.0 + LN2)) < 1e-15
    assert abs(OMEGA_DE_TODAY - r / (1.0 - r)) < 1e-15
    assert W == -1.0
    assert K_CURVATURE == 0


def test_discrete_12_fold_approximates_ln2():
    assert abs(R_DISCRETE - 49 / 120) < 1e-15
    assert abs(T_DISCRETE - 49 / 71) < 1e-15
    assert abs(T_DISCRETE - LN2) < 0.004
    assert abs(R_DISCRETE - r) < 0.002


def test_flatness_closes():
    assert abs(OMEGA_M - (1.0 - OMEGA_DE_TODAY)) < 1e-15
    assert abs(OMEGA_M + OMEGA_DE_TODAY - 1.0) < 1e-15


def test_not_the_planck_comparison_value():
    assert abs(OMEGA_DE_TODAY - 0.689) > 5e-4

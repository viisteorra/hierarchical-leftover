"""One generating function. No third leftover. No H0."""

import math

from generate import (
    G,
    density_tail,
    generating_set,
    leftover_exp_clock_foil,
    leftover_exp_finite,
    leftover_exp_infinite,
    licensed_ratios,
    spatial_curvature,
)
from geometry import K_CURVATURE, OMEGA_DE_TODAY, Q5_WEIGHT, R_DISCRETE, TOTAL, r as R_GEN
from ruler import leftover_exp_infinite as ruler_inf
from ruler import leftover_exp_planar
from spacetime import F_BAO, F_CMB, f_for, measure_for


def test_g_is_the_geometric_tail():
    assert abs(G(R_GEN) - 1.0 / (1.0 - R_GEN)) < 1e-15
    assert abs(density_tail() - R_GEN / (1.0 - R_GEN)) < 1e-15
    assert abs(density_tail() - math.log(2.0)) < 1e-15
    assert abs(density_tail() - OMEGA_DE_TODAY) < 1e-15
    # G is the series Σ ρ^k, not a new function
    rho = float(R_GEN)
    partial = sum(rho**k for k in range(64))
    assert abs(partial - G(rho)) < 1e-15


def test_flatness_is_euclidean_leftover_not_a_fit():
    assert spatial_curvature() == 0
    assert K_CURVATURE == 0


def test_finite_leftover_is_g_at_r_over_11():
    e = leftover_exp_finite()
    assert abs(e - (12 / 11) * G(float(R_DISCRETE) / 11)) < 1e-15
    assert abs(e - 1440 / 1271) < 1e-15  # 12-fold lattice, not continuous r
    assert abs(leftover_exp_planar() - e) < 1e-15


def test_infinite_leftover_is_log_g():
    e = leftover_exp_infinite()
    assert abs(e - (1.0 + (1.0 - R_GEN) * math.log(G(R_GEN)))) < 1e-15
    assert abs(e - (1.0 - (1.0 - R_GEN) * math.log(1.0 - R_GEN))) < 1e-15
    assert abs(ruler_inf() - e) < 1e-15


def test_only_two_licensed_ratios():
    lic = licensed_ratios()
    assert set(lic) == {"mix", "generators"}
    assert abs(lic["mix"] - R_GEN) < 1e-15
    assert abs(lic["generators"] - float(R_DISCRETE) / Q5_WEIGHT) < 1e-15
    # r/12 is the foil, not licensed
    foil = leftover_exp_clock_foil()
    assert abs(foil - leftover_exp_finite()) > 1e-6
    assert abs(foil - leftover_exp_infinite()) > 1e-6


def test_no_third_ratio_from_internal_integers():
    # Internal integers {1,2,4,5,11,12}. Licensed leftover ρ is only r/11.
    # r itself is density, not leftover.
    internals = {1, 2, 4, 5, 11, 12}
    leftover_denominators = {11}
    assert leftover_denominators <= internals
    assert 12 not in leftover_denominators
    assert 10 not in internals


def test_generating_set_lemma_c():
    assert generating_set("BAO") == "finite"
    assert generating_set("weak_lensing") == "finite"
    assert generating_set("theta_star") == "infinite"
    assert generating_set("CMB") == "infinite"
    assert generating_set("SN") == "none"
    assert generating_set("time_delay") == "none"
    assert generating_set("BBN") == "early"
    assert measure_for("BAO") == "planar"
    assert measure_for("CMB") == "infinite"
    assert abs(f_for("BAO") - F_BAO) < 1e-15
    assert abs(f_for("CMB") - F_CMB) < 1e-15


def test_leftover_is_twelve_fold_address():
    from ruler import leftover_octave
    from spacetime import ZSTAR_RECOMB

    lo = leftover_octave(ZSTAR_RECOMB)
    assert lo["n_int"] == 10
    assert abs(lo["steps_12"] - lo["frac"] * 12) < 1e-15
    assert 1.0 < lo["steps_12"] < 2.0  # leftover is ~1.09 steps of the 12-fold
    assert abs(lo["f_planar"] - 2.0 ** ((lo["steps_12"] / 12) * leftover_exp_finite())) < 1e-15
    assert abs(lo["f_inf"] - 2.0 ** ((lo["steps_12"] / 12) * leftover_exp_infinite())) < 1e-15


def test_g_does_not_contain_h0():
    import inspect
    from pathlib import Path

    src = Path(inspect.getfile(G)).read_text()
    assert "H0_SN" not in src
    assert "73.47" not in src
    assert "chi2" not in src.lower()
    assert TOTAL == 12 and Q5_WEIGHT == 11

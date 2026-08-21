"""No free cosmological knobs. Every symbol traces to P, elsewhere, or a measurement."""

from pathlib import Path

from geometry import OMEGA_DE_TODAY, Q4_WEIGHT, Q5_WEIGHT, TOTAL, W, r
from ruler import leftover_exp_infinite, leftover_exp_planar, leftover_f, leftover_octave
from spacetime import F_BAO, F_CMB, ZSTAR_RECOMB, f_for, measure_for
from uniqueness import (
    PERIOD_N,
    SEED,
    HYPER_STEPS,
    T_FRAC,
    derivation_chain,
    period_from_fifth_convergent,
    square_from_hypercubes,
)

ROOT = Path(__file__).resolve().parents[1]


def test_measure_c_is_observable_type_not_chi2():
    assert measure_for("BAO") == "planar"
    assert measure_for("DV") == "planar"
    assert measure_for("theta_star") == "infinite"
    assert measure_for("CMB") == "infinite"
    assert abs(f_for("BAO") - F_BAO) < 1e-15
    assert abs(f_for("theta_star") - F_CMB) < 1e-15
    assert abs(f_for("BAO") - f_for("CMB")) > 0.005


def test_swapped_measures_are_a_different_assignment():
    # Using the null-metric leftover on BAO is not Lemma C.
    assert abs(f_for("CMB") - F_BAO) > 1e-4
    assert abs(f_for("BAO") - F_CMB) > 1e-4


def test_geometry_and_ruler_are_not_fits():
    for rel in ("src/geometry.py", "src/ruler.py", "src/uniqueness.py", "src/generate.py"):
        text = (ROOT / rel).read_text().lower()
        assert "chi2" not in text
        assert "minimize" not in text
        assert "73.47" not in text
        assert "h0_sn" not in text


def test_no_h0_in_leftover_formula():
    a = leftover_f(ZSTAR_RECOMB)
    b = leftover_octave(ZSTAR_RECOMB)["f_planar"]
    # Changing a dummy H0 cannot exist; leftover args are z* only.
    assert leftover_f(ZSTAR_RECOMB) == a
    assert leftover_octave(ZSTAR_RECOMB)["f_planar"] == b
    # Completed octaves drop out.
    phi = leftover_octave(ZSTAR_RECOMB)["frac"]
    for n in range(8, 14):
        z = 2.0 ** (n + phi) - 1.0
        assert abs(leftover_f(z) - a) < 1e-12
        assert abs(leftover_octave(z)["f_planar"] - b) < 1e-12


def test_exponents_are_the_lemmas():
    from generate import G
    from geometry import R_DISCRETE as _RD

    assert abs(leftover_exp_planar() - (12 / 11) * G(float(_RD) / 11)) < 1e-15
    assert abs(leftover_exp_planar() - 1440 / 1271) < 1e-15
    e = leftover_exp_infinite()
    assert abs(e - (1.0 - (1.0 - r) * __import__("math").log(1.0 - r))) < 1e-15


def test_lock_is_the_u3_algebra():
    assert (Q4_WEIGHT, Q5_WEIGHT, TOTAL) == (1, 11, 12)
    assert (PERIOD_N, SEED, HYPER_STEPS) == (12, 1, 11)
    import math

    assert abs(float(T_FRAC) - 49 / 71) < 1e-15
    assert abs(OMEGA_DE_TODAY - math.log(2.0)) < 1e-15
    assert W == -1.0
    assert square_from_hypercubes() == (4, 5)
    assert period_from_fifth_convergent()[0] == 12


def test_first_principles_document():
    text = (ROOT / "docs" / "FIRST_PRINCIPLES.md").read_text()
    assert "resonat" in text.lower()
    assert "7/12" in text
    assert "1:11" in text
    assert "not arbitrary" in text.lower() or "not a fit" in text.lower()
    assert "How to take it further" in text
    assert "log_2(3/2)" in text or r"\log_2(3/2)" in text


def test_chain_has_one_primitive():
    chain = derivation_chain()
    assert sum(1 for row in chain if row["free"]) == 1
    assert chain[0]["symbol"] == "P"
    origins = " ".join(row["origin"] for row in chain)
    assert "χ²" not in origins and "chi2" not in origins.lower()
    assert "H0" not in "".join(str(row["value"]) for row in chain if row["symbol"] != "H0 map")

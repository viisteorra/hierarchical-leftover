"""The geometric core cannot silently become a fit parameter."""

from pathlib import Path

from geometry import OMEGA_DE_TODAY, Q4_WEIGHT, Q5_WEIGHT, r
from ruler import leftover_exp_infinite, leftover_f, leftover_octave
from spacetime import F_ANSATZ_T10, F_AXIOM, F_CLOCK, F_MIXED_2D, ZSTAR_RECOMB

SRC = Path(__file__).resolve().parents[1] / "src" / "geometry.py"
ST = Path(__file__).resolve().parents[1] / "src" / "spacetime.py"


def test_source_still_hardcodes_the_1_11_rule():
    text = SRC.read_text()
    assert "Q4_WEIGHT = 1" in text
    assert "Q5_WEIGHT = 11" in text
    assert "OMEGA_DE_TODAY = r / (1 - r)" in text
    assert "minimize" not in text
    assert "chi2" not in text.lower()


def test_exact_rational_tail():
    assert abs(OMEGA_DE_TODAY - 49 / 71) < 1e-12
    assert Q4_WEIGHT / Q5_WEIGHT == 1 / 11
    assert abs(r - 49 / 120) < 1e-12


def test_leftover_theorem_is_infinite_euclidean_tail():
    """A2: 2^{frac · [1-(1-r)ln(1-r)]}, not 4D cutoff, not T/10."""
    lo = leftover_octave(ZSTAR_RECOMB)
    assert abs(F_CLOCK - 2.0 ** lo["frac"]) < 1e-12
    assert abs(F_AXIOM - 2.0 ** (lo["frac"] * leftover_exp_infinite())) < 1e-15
    assert abs(F_AXIOM - leftover_f(ZSTAR_RECOMB)) < 1e-15
    assert abs(F_AXIOM - F_CLOCK) > 1e-4
    assert abs(F_AXIOM - F_MIXED_2D) > 1e-4
    st = ST.read_text()
    assert "12/11" in st
    assert "ln(1-r)" in st or "ln(1−r)" in st or "(1-r)" in st
    assert abs(F_ANSATZ_T10 - 710 / 661) < 1e-12
    assert abs(F_AXIOM - F_ANSATZ_T10) > 1e-4

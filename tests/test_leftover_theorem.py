"""Leftover theorem: infinite-D tail-weighted Euclidean mean. No H0, no T/N, no χ²."""

import math
from pathlib import Path

from geometry import OMEGA_DE_TODAY, r as R_GEN
from ruler import (
    DIM_EXP_1,
    DIM_EXP_2,
    DIM_EXP_3,
    DIM_EXP_4,
    leftover_exp_infinite,
    leftover_f,
    leftover_octave,
    mix_compress_f,
    mix_stretch_f,
    tail_per_octave_f,
)
from spacetime import (
    F_ANSATZ_T10,
    F_AXIOM,
    F_BAO,
    F_CLOCK,
    F_CMB,
    F_FINITE4,
    F_MIXED_2D,
    F_PLANCK_Z,
    PLANCK_ZSTAR,
    ZSTAR_RECOMB,
)

RULER = Path(__file__).resolve().parents[1] / "src" / "ruler.py"
ST = Path(__file__).resolve().parents[1] / "src" / "spacetime.py"
TH = Path(__file__).resolve().parents[1] / "docs" / "THEOREM.md"


def test_closed_form_matches_truncated_series():
    assert abs(DIM_EXP_1 - 1) < 1e-15
    assert abs(DIM_EXP_2 - 12 / 11) < 1e-15
    assert abs(DIM_EXP_3 - 3 / 2) < 1e-15
    assert abs(DIM_EXP_4 - 4 / 3) < 1e-15
    e = leftover_exp_infinite()
    assert abs(e - (1.0 - (1.0 - R_GEN) * math.log(1.0 - R_GEN))) < 1e-15
    num = den = 0.0
    for d in range(1, 64):
        w = R_GEN ** (d - 1)
        ed = 1.0 if d == 1 else d / (d - 1)
        num += w * ed
        den += w
    assert abs(e - num / den) < 1e-12


def test_two_measure_bao_is_tailed_planar_cmb_is_infinite():
    lo = leftover_octave(ZSTAR_RECOMB)
    from ruler import leftover_exp_planar
    assert abs(F_BAO - 2.0 ** (lo["frac"] * leftover_exp_planar())) < 1e-15
    from generate import leftover_exp_finite

    assert abs(leftover_exp_planar() - leftover_exp_finite()) < 1e-15
    assert abs(F_BAO - mix_stretch_f(lo["frac"])) > 1e-4  # tail is not bare 12/11
    assert abs(F_CMB - leftover_f(ZSTAR_RECOMB)) < 1e-15
    assert abs(F_BAO - F_CMB) > 1e-4
    assert abs(F_AXIOM - F_CMB) < 1e-15


def test_clock_2d_finite4_and_infinite():
    lo = leftover_octave(ZSTAR_RECOMB)
    assert lo["n_int"] == 10
    assert abs(F_CLOCK - 1090.84 / 1024) < 1e-15
    assert abs(F_MIXED_2D - mix_stretch_f(lo["frac"])) < 1e-15
    assert abs(F_FINITE4 - 2.0 ** (lo["frac"] * 325 / 264)) < 1e-15
    assert abs(F_AXIOM - 2.0 ** (lo["frac"] * leftover_exp_infinite())) < 1e-15
    assert abs(F_AXIOM - leftover_f(ZSTAR_RECOMB)) < 1e-15
    assert abs(F_AXIOM - lo["f_inf"]) < 1e-15
    assert abs(F_AXIOM - F_CLOCK) > 1e-4
    assert abs(F_AXIOM - F_MIXED_2D) > 1e-4
    assert abs(F_AXIOM - F_FINITE4) > 1e-4


def test_direction_inverses_are_not_the_theorem():
    lo = leftover_octave(ZSTAR_RECOMB)
    assert abs(leftover_f(ZSTAR_RECOMB) - mix_compress_f(lo["frac"])) > 1e-4
    assert mix_compress_f(lo["frac"]) < F_CLOCK < F_MIXED_2D


def test_independent_of_completed_octaves():
    frac = leftover_octave(ZSTAR_RECOMB)["frac"]
    fs = [leftover_f(2.0 ** (n + frac) - 1.0) for n in range(8, 14)]
    assert max(fs) - min(fs) < 1e-12
    assert abs(fs[0] - F_AXIOM) < 1e-12


def test_leftover_f_source_has_no_h0_or_omega():
    src = RULER.read_text().split("def leftover_f")[1].split("def remainder_k")[0]
    assert "OMEGA" not in src
    assert "H0" not in src
    assert leftover_f.__code__.co_varnames[: leftover_f.__code__.co_argcount] == ("z_star",)
    assert leftover_f(ZSTAR_RECOMB) != OMEGA_DE_TODAY


def test_ansatz_is_not_the_theorem():
    assert abs(F_ANSATZ_T10 - 1.0 / (1.0 - OMEGA_DE_TODAY / 10.0)) < 1e-12
    assert abs(F_ANSATZ_T10 - tail_per_octave_f()) < 1e-15
    assert abs(F_AXIOM - F_ANSATZ_T10) > 1e-4


def test_recombination_not_bao_hubble_ratio():
    empirical = 73.47106273986928 / 68.38338180546
    assert abs(F_AXIOM - empirical) > 1e-3
    assert abs(F_PLANCK_Z - leftover_f(PLANCK_ZSTAR)) < 1e-15


def test_theorem_document_and_axiom_lock():
    th = TH.read_text()
    st = ST.read_text()
    assert r"1-(1-r)\ln(1-r)" in th or "1 − (1−r)" in th or "1-(1-r)" in th
    assert "710/661" in th
    assert "in the proof" in th and "H_0" in th
    assert "12/11" in st
    assert "ln(1-r)" in st or "ln(1−r)" in st or "(1-r)" in st

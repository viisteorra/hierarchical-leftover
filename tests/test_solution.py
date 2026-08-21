"""Frozen solution numbers. Given P; no χ² in the formulas."""

from pathlib import Path

from geometry import OMEGA_DE_TODAY, W
from hubble import solution
from spacetime import F_BAO, F_CMB, ZSTAR_RECOMB
from uniqueness import T_FRAC


def test_density_and_w():
    import math

    assert abs(OMEGA_DE_TODAY - math.log(2.0)) < 1e-15
    assert abs(float(T_FRAC) - 49 / 71) < 1e-15
    assert W == -1.0


def test_two_measure_maps_land_on_sn():
    s = solution()
    assert abs(F_BAO * s["h0_bao"] - s["h0_sn"]) < 0.25
    assert abs(F_CMB * s["h0_th"] - s["h0_sn"]) < 0.25
    assert abs(F_BAO - F_CMB) > 0.005
    assert s["map_bao_sigma"] < 0.3
    assert s["map_th_sigma"] < 0.3
    assert s["early_early_sigma"] < 1.0


def test_solved_document():
    text = (Path(__file__).resolve().parents[1] / "docs" / "SOLVED.md").read_text()
    assert "49" in text and "71" in text
    assert "1440" in text and "1271" in text
    assert "one leftover" in text.lower() or "One leftover" in text
    assert str(ZSTAR_RECOMB) in text or "1089.84" in text

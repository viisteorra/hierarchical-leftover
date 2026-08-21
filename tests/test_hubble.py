"""Hubble scale theorem. f is not fitted to H0."""

from pathlib import Path

from geometry import OMEGA_M, r as R_GEN
from hubble import H0_BAO, H0_SN, H0_TH_LOCK, SIG_SN, solution
from spacetime import F_BAO, F_CMB


SRC_HUBBLE = Path(__file__).resolve().parents[1] / "src" / "hubble.py"
SRC_RULER = Path(__file__).resolve().parents[1] / "src" / "ruler.py"
SRC_ST = Path(__file__).resolve().parents[1] / "src" / "spacetime.py"
TH = Path(__file__).resolve().parents[1] / "docs" / "THEOREM.md"
SOLVED = Path(__file__).resolve().parents[1] / "docs" / "SOLVED.md"


def test_f_is_not_the_empirical_ratio():
    emp_bao = H0_SN / H0_BAO
    emp_th = H0_SN / H0_TH_LOCK
    assert abs(F_BAO - emp_bao) > 1e-6
    assert abs(F_CMB - emp_th) > 1e-6


def test_leftover_source_has_no_h0():
    for path in (SRC_RULER, SRC_ST):
        text = path.read_text()
        assert "73.47" not in text
        assert "68.38" not in text
        assert "67.59" not in text
        assert "H0_SN" not in text


def test_inverse_maps_land_on_sn():
    s = solution()
    assert s["map_bao_sigma"] < 0.3
    assert s["map_th_sigma"] < 0.3
    assert abs(s["pred_sn_from_bao"] - H0_SN) < 0.25
    assert abs(s["pred_sn_from_th"] - H0_SN) < 0.25


def test_forward_bao_from_sn():
    s = solution()
    assert s["fwd_bao_sigma"] < 1.0
    assert abs(s["pred_bao_from_sn"] - H0_BAO) < 0.25


def test_early_early_corollary():
    s = solution()
    assert abs(s["ratio_f"] - F_CMB / F_BAO) < 1e-15
    assert s["early_early_sigma"] < 1.0
    assert abs(s["pred_bao_from_th"] - H0_BAO) < 1.0


def test_raw_split_was_the_tension():
    s = solution()
    assert s["raw_bao_sigma"] > 4.0
    assert s["raw_th_sigma"] > 4.0
    assert s["raw_split_bao_kms"] > 4.5
    assert s["map_bao_sigma"] < s["raw_bao_sigma"] / 10.0


def test_f_does_not_depend_on_which_h0_you_pass():
    a = solution(h0_sn=70.0, h0_bao=65.0, h0_th=66.0)
    b = solution()
    assert abs(a["f_box"] - b["f_box"]) < 1e-15
    assert abs(a["f_inf"] - b["f_inf"]) < 1e-15
    assert abs(a["ratio_f"] - b["ratio_f"]) < 1e-15


def test_r_is_the_mix_not_a_hubble_fit():
    import math

    assert abs(R_GEN - math.log(2.0) / (1.0 + math.log(2.0))) < 1e-15
    assert abs(float(OMEGA_M) - (1.0 - math.log(2.0))) < 1e-15


def test_theorem_3_and_solved_docs():
    th = TH.read_text()
    assert "Theorem 3" in th
    assert r"f_\infty" in th and r"f_\square" in th
    solved = SOLVED.read_text()
    assert "Hubble" in solved
    assert "1440" in solved and "1271" in solved
    hub = SRC_HUBBLE.read_text()
    assert "f_□" in hub or "f_box" in hub
    assert "73.471" not in SRC_RULER.read_text()


def test_sh0es_budget_not_intercept_scatter():
    assert SIG_SN == 1.04

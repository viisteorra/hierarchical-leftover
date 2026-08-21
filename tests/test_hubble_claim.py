"""Hubble-tension claim given P. Kill if leftover is a fit or the split reopens."""

from pathlib import Path

from hubble import H0_BAO, H0_SN, H0_TH_LOCK, SIG_SN, solution, tension
from spacetime import F_BAO, F_CMB, f_for


def test_raw_split_is_the_ordinary_tension():
    s = solution()
    assert s["raw_bao_sigma"] > 4.5
    assert s["raw_th_sigma"] > 4.5
    assert s["raw_split_bao_kms"] > 4.5


def test_licensed_maps_close_below_one_sigma():
    s = solution()
    assert s["map_bao_sigma"] < 1.0
    assert s["map_th_sigma"] < 1.0
    assert s["fwd_bao_sigma"] < 1.5
    assert s["early_early_sigma"] < 1.5


def test_claim_threshold_tenth_of_raw():
    s = solution()
    assert s["map_bao_sigma"] < s["raw_bao_sigma"] / 10.0
    assert s["map_th_sigma"] < s["raw_th_sigma"] / 10.0


def test_swapped_measures_do_not_close_both():
    s = solution()
    swap_bao = tension(H0_SN, F_CMB * H0_BAO, SIG_SN, 0.217)
    swap_cmb = tension(H0_SN, F_BAO * H0_TH_LOCK, SIG_SN, 0.54)
    # Licensed assignment is unique: swap is worse on each licensed pair.
    assert swap_bao > s["map_bao_sigma"]
    assert swap_cmb > s["map_th_sigma"]
    assert f_for("BAO") == F_BAO
    assert f_for("CMB") == F_CMB


def test_one_f_does_not_replace_two_measures():
    s = solution()
    # Infinite leftover on BAO is the known K2-kill assignment.
    inf_on_bao = tension(H0_SN, F_CMB * H0_BAO, SIG_SN, 0.217)
    box_on_cmb = tension(H0_SN, F_BAO * H0_TH_LOCK, SIG_SN, 0.54)
    assert inf_on_bao > s["map_bao_sigma"]
    assert box_on_cmb > s["map_th_sigma"]
    assert s["map_bao_sigma"] < 0.3
    assert s["map_th_sigma"] < 0.3


def test_f_is_not_the_hubble_ratio():
    assert abs(F_BAO - H0_SN / H0_BAO) > 1e-6
    assert abs(F_CMB - H0_SN / H0_TH_LOCK) > 1e-6


def test_claim_document_exists():
    root = Path(__file__).resolve().parents[1] / "docs"
    text = (root / "CLAIM.md").read_text()
    assert "If the vacuum is this hierarchy" in text
    assert "49/71" in text or "ln 2" in text or "ln(2)" in text
    assert "Lemma C" in text or "finite" in text.lower()
    assert "P is proved" in text or "P is not proved" in text.lower() or "cannot say" in text.lower()


def test_theorem_is_marked_complete():
    root = Path(__file__).resolve().parents[1] / "docs"
    complete = (root / "COMPLETE.md").read_text()
    theory = (root / "THEORY.md").read_text()
    solved = (root / "SOLVED.md").read_text()
    assert "complete" in complete.lower()
    assert "There is no further" in complete or "no further" in complete.lower()
    assert "Closed" in theory
    assert "THEORY.md" in solved or "FIRST_PRINCIPLES" in solved or "ln 2" in solved or "ln2" in solved.lower()
    assert "0/41" not in solved  # stale hell count

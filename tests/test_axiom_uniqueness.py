"""U1–U3: 12, {4,5}, one seed → 1:11. No cosmology."""

import math

from geometry import OMEGA_DE_TODAY, Q4_WEIGHT, Q5_WEIGHT, TOTAL, r
from uniqueness import (
    FIFTH_CENTS,
    FIFTH_STEPS,
    HYPER_STEPS,
    PERIOD_N,
    Q_EUCLID,
    Q_HYPER_FIRST,
    R_FRAC,
    SEED,
    T_FRAC,
    continued_fraction_fifth,
    derivation_chain,
    fifth_cents,
    fifth_cut_basin,
    mixed_fills,
    mix_r,
    one_seed_fill,
    period_from_fifth_convergent,
    smallest_nfold_fifth,
    square_from_hypercubes,
    tail,
)


def test_u1_twelve_is_smallest_nfold_fifth():
    n, k, c = smallest_nfold_fifth()
    assert (n, k) == (12, 7)
    assert abs(c - FIFTH_CENTS) < 1e-9
    assert abs(c) < 2.0
    for m in range(2, 12):
        assert abs(fifth_cents(m)) >= 2.0


def test_u1_is_convergent_not_a_cents_cut():
    n, k, c = period_from_fifth_convergent()
    assert (n, k) == (12, 7)
    assert (n, k) == (PERIOD_N, FIFTH_STEPS)
    conv = continued_fraction_fifth()
    assert (3, 5) in conv and (7, 12) in conv
    # 4/7 is a semi-convergent, not a best approximation — not the period
    assert (4, 7) not in conv
    lo, hi = fifth_cut_basin()
    assert lo < 2.0 < hi
    assert 15.0 < hi < 17.0  # 4/7 ≈ 16.2 cents; 3/5 is 18 and is not the next n-fold
    n_lo, _, _ = smallest_nfold_fifth(cut=lo + 1e-6)
    n_mid, _, _ = smallest_nfold_fifth(cut=0.5 * (lo + hi))
    n_hi, _, _ = smallest_nfold_fifth(cut=hi - 1e-6)
    assert n_lo == n_mid == n_hi == 12


def test_u1_seven_twelve_is_a_convergent():
    conv = continued_fraction_fifth()
    assert (7, 12) in conv
    # previous convergent is worse than 2 cents
    assert abs(fifth_cents(5, 3)) > 2.0


def test_u2_square_boundary():
    assert Q_EUCLID == 4
    assert Q_HYPER_FIRST == 5
    # {4,q} Euclidean iff 2(q-2)==4
    assert 2 * (Q_EUCLID - 2) == 4
    assert 2 * (Q_HYPER_FIRST - 2) > 4
    assert 2 * (Q_HYPER_FIRST - 1 - 2) == 4  # q=4 is the boundary, q=5 first over
    assert square_from_hypercubes() == (4, 5)


def test_u3_one_seed_is_1_11():
    assert one_seed_fill(12) == (1, 11)
    assert (SEED, HYPER_STEPS) == (1, 11)
    assert PERIOD_N == 12 == TOTAL
    assert Q4_WEIGHT == 1 and Q5_WEIGHT == 11


def test_u3_tail_is_49_over_71():
    assert T_FRAC == 49 / 71 or T_FRAC == __import__("fractions").Fraction(49, 71)
    assert abs(float(T_FRAC) - 49 / 71) < 1e-15
    assert abs(float(R_FRAC) - 49 / 120) < 1e-15  # 1:11 fill algebra, not Ω_DE
    import math

    assert abs(OMEGA_DE_TODAY - math.log(2.0)) < 1e-15
    assert abs(r - math.log(2.0) / (1.0 + math.log(2.0))) < 1e-15
    fills = mixed_fills(12)
    assert fills[0] == (1, 11)
    assert (2, 10) in fills and (6, 6) in fills
    assert len(fills) == 11
    # other mixed fills exist and are different tails — not selected by χ²
    assert tail(2, 10) != tail(1, 11)
    assert tail(6, 6) != tail(1, 11)


def test_axioms_document_exists():
    from pathlib import Path

    ax = Path(__file__).resolve().parents[1] / "docs" / "AXIOMS.md"
    text = ax.read_text()
    assert "exactly one" in text.lower() or "Exactly one" in text
    assert "7/12" in text
    assert "{4,5}" in text or "q=5" in text
    assert "one primitive" in text.lower() or "The one primitive" in text
    assert "sound horizon" in text or "sound-horizon" in text
    assert "resonat" in text.lower()


def test_only_primitive_is_p():
    chain = derivation_chain()
    free = [row for row in chain if row["free"]]
    assert len(free) == 1
    assert free[0]["symbol"] == "P"
    assert any(row["symbol"] == "period N" and row["value"] == 12 for row in chain)
    assert any(row["symbol"] == "T=Ω_DE" and "ln(2)" in str(row["value"]) for row in chain)
    assert any(row["symbol"] == "G" for row in chain)
    assert any(row["symbol"] == "k" and row["value"] == 0 for row in chain)

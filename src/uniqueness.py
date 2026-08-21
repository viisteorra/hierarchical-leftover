"""Uniqueness lemmas for 12, {4,5}, and one-seed → 1:11.

No H0, no χ², no leftover. Combinatorics + tilings + continued fractions.
The 2-cent figure is a *check* inside a basin, not an input. U1 is the
continued-fraction convergent of log2(3/2) after 3/5, which is 7/12.
"""

from __future__ import annotations

import math
from fractions import Fraction

LOG2_FIFTH = math.log2(3 / 2)  # log2(3/2)
# Check only: any naive cut in (err(7/12), next n-fold) ≈ (1.96, 16.2) yields n=12.
# Not an input to PERIOD_N.
FIFTH_CENTS_CUT = 2.0

# Square tilings {4,q}: q squares meet at a vertex.
# Why squares: leftover uses Euclidean hypercubes q=2d; 2D ⇒ q=4.
Q_EUCLID = 4  # {4,4} Euclidean
Q_HYPER_FIRST = 5  # {4,5} first hyperbolic


def fifth_cents(n: int, k: int | None = None) -> float:
    """Cents of the interval 2^{k/n} relative to 3:2."""
    if k is None:
        k = int(round(n * LOG2_FIFTH))
    return 1200.0 * (k / n - LOG2_FIFTH)


def smallest_nfold_fifth(cut: float = FIFTH_CENTS_CUT) -> tuple[int, int, float]:
    """Smallest n>1 with some k/n within `cut` cents of the fifth.

    Returns (n, k, cents). Unique at cut=2: n=12, k=7, −1.955 cents.
    """
    for n in range(2, 200):
        k = int(round(n * LOG2_FIFTH))
        c = fifth_cents(n, k)
        if abs(c) < cut:
            return n, k, c
    raise RuntimeError("no n-fold fifth under cut")


def period_from_fifth_convergent(terms: int = 12) -> tuple[int, int, float]:
    """U1. Period is the denominator of the first CF convergent of log2(3/2) after 3/5.

    Convergents of log2(3/2): 1/1, 1/2, 3/5, 7/12, 24/41, 31/53, …
    Semi-convergents (4/7, …) are not best approximations and are not used.
    No cents threshold is an input. Returns (n, k, cents) = (12, 7, −1.955…).
    """
    conv = continued_fraction_fifth(terms=terms)
    seen_three_five = False
    for k, n in conv:
        if (k, n) == (3, 5):
            seen_three_five = True
            continue
        if seen_three_five and n > 5:
            return n, k, fifth_cents(n, k)
    raise RuntimeError("7/12 convergent missing")


def fifth_cut_basin() -> tuple[float, float]:
    """Open interval of naive cents-cuts that all select n=12.

    Upper edge is the next-best n-fold among n=2..11 (4/7 at ~16.2 cents),
    not 3/5, because 4/7 is a semi-convergent the CF proof already discards.
    2.0 sits inside; it is not fitted. U1 does not use this interval as input.
    """
    lo = abs(fifth_cents(12, 7))
    hi = min(abs(fifth_cents(n)) for n in range(2, 12))
    return lo, hi


def square_from_hypercubes() -> tuple[int, int]:
    """U2. 2D Euclidean hypercube has q=2d=4. First hyperbolic square tiling is q=5."""
    d = 2
    q4 = 2 * d
    q5 = q4 + 1
    return q4, q5


def continued_fraction_fifth(terms: int = 8) -> list[tuple[int, int]]:
    """Convergents of log2(3/2). 7/12 is the first with |cents|<2."""
    x = LOG2_FIFTH
    a: list[int] = []
    for _ in range(terms):
        ai = math.floor(x)
        a.append(ai)
        frac = x - ai
        if frac < 1e-15:
            break
        x = 1.0 / frac
    nums = [1, a[0]]
    dens = [0, 1]
    out = [(a[0], 1)]
    for i in range(1, len(a)):
        nums.append(a[i] * nums[-1] + nums[-2])
        dens.append(a[i] * dens[-1] + dens[-2])
        out.append((nums[-1], dens[-1]))
    return out


def one_seed_fill(period: int = 12) -> tuple[int, int]:
    """Exactly one Euclidean seed in a period of `period` parts: (1, period-1)."""
    if period < 2:
        raise ValueError("period")
    return 1, period - 1


def mix_r(n4: int, n5: int, q4: int = Q_EUCLID, q5: int = Q_HYPER_FIRST) -> Fraction:
    tot = n4 + n5
    return (n4 * Fraction(2, q4) + n5 * Fraction(2, q5)) / tot


def tail(n4: int, n5: int) -> Fraction:
    rr = mix_r(n4, n5)
    return rr / (1 - rr)


def mixed_fills(period: int = 12) -> list[tuple[int, int]]:
    """All mixed fills n4:n5 of a period. One-seed is (1, period-1) only."""
    return [(n4, period - n4) for n4 in range(1, period)]


PERIOD_N, FIFTH_STEPS, FIFTH_CENTS = period_from_fifth_convergent()
SEED, HYPER_STEPS = one_seed_fill(PERIOD_N)
R_FRAC = mix_r(SEED, HYPER_STEPS)
T_FRAC = tail(SEED, HYPER_STEPS)


def derivation_chain() -> list[dict]:
    """Every symbol in the theory. free=False except the named primitive P.

    z* is a measurement (recombination), not a hierarchy integer and not H0.
    """
    q4, q5 = square_from_hypercubes()
    return [
        {"symbol": "P", "value": "vacuum = this hierarchy; tail = Ω_DE", "origin": "primitive", "free": True},
        {"symbol": "octave", "value": "2:1", "origin": "elsewhere: linear waves, given P", "free": False},
        {"symbol": "fifth", "value": "3:2", "origin": "elsewhere: next integer ratio", "free": False},
        {"symbol": "period N", "value": PERIOD_N, "origin": "U1: CF convergent of log2(3/2) after 3/5", "free": False},
        {"symbol": "q_euclid", "value": q4, "origin": "U2: 2D hypercube q=2d", "free": False},
        {"symbol": "q_hyper", "value": q5, "origin": "U2: first hyperbolic square {4,5}", "free": False},
        {"symbol": "fill", "value": f"{SEED}:{HYPER_STEPS}", "origin": "U3: one Euclidean seed", "free": False},
        {"symbol": "r", "value": "ln2/(1+ln2)", "origin": "continuous: T=ln2 ⇒ r=T/(1+T)", "free": False},
        {"symbol": "G", "value": "1/(1−ρ)", "origin": "geometric tail of the hierarchy", "free": False},
        {"symbol": "T=Ω_DE", "value": "ln(2)", "origin": "continuous binary-scale residual", "free": False},
        {"symbol": "12-fold", "value": "49/71", "origin": "U1–U3 rational approximation to ln(2)", "free": False},
        {"symbol": "k", "value": 0, "origin": "U2: metric leftover is Euclidean hypercubes", "free": False},
        {"symbol": "Ω_m", "value": "1-ln(2)", "origin": "flatness k=0 ⇒ 1−T", "free": False},
        {"symbol": "w", "value": -1, "origin": "lemma: period does not drift", "free": False},
        {"symbol": "z*", "value": "photon decoupling", "origin": "B: measured second scale, not H0", "free": False},
        {"symbol": "φ", "value": "{log2(1+z*)}", "origin": "clock leftover identity", "free": False},
        {"symbol": "e_□", "value": "(12/11)G(r/11)", "origin": "finite generating set", "free": False},
        {"symbol": "e_∞", "value": "1+(1-r)ln G(r)", "origin": "infinite generating set", "free": False},
        {"symbol": "f_□", "value": "2^{φ e_□}", "origin": "C: finite set (spatial slices / BAO)", "free": False},
        {"symbol": "f_∞", "value": "2^{φ e_∞}", "origin": "C: infinite set (null metric / θ*)", "free": False},
        {"symbol": "H0 map", "value": "H0_L = f H0_E", "origin": "Theorem 3: two readings, one metric", "free": False},
    ]

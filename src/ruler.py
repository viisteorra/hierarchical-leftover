"""Octave-remainder ruler. Does not retune ln(2) or put H0 into f.

Density and leftover are evaluations of G(ρ)=1/(1−ρ) in generate.py.
This module turns φ={log2(1+z*)} into f=2^{φ e}.
"""

from __future__ import annotations

import math

from generate import leftover_exp_finite, leftover_exp_infinite as _e_inf
from geometry import OMEGA_DE_TODAY, Q5_WEIGHT, TOTAL, r as R_GEN

PERIOD = int(TOTAL)          # 12
Q_EUCLID = 4                 # square coordination, not the 1:11 *weight*
Q_HYPER = 5
UNFOLD_B = PERIOD * Q_EUCLID  # 48

# Euclidean hypercube coordination q=2d. Consecutive ratio q_d/q_{d-1}=d/(d-1).
# Leftover exponents: 1D clock, 2D mix 12/11, 3D space 3/2, 4D spacetime 4/3.
# Mean is 325/264. G1 is logarithmic ⇒ f = 2^{φ · ⟨e⟩} = geom mean of the four f_d.
DIM_EXP_1 = 1                      # clock
DIM_EXP_2 = PERIOD / Q5_WEIGHT     # 12/11
DIM_EXP_3 = 3 / 2                  # q_3/q_2 = 6/4
DIM_EXP_4 = 4 / 3                  # q_4/q_3 = 8/6
DIM_EXP_MEAN_4 = (DIM_EXP_1 + DIM_EXP_2 + DIM_EXP_3 + DIM_EXP_4) / 4  # 325/264, 4D cutoff, not A2


def leftover_exp_infinite() -> float:
    """Infinite generating set: e_∞ = 1 + (1−r) ln G(r). Same number as before."""
    return _e_inf()


def leftover_exp_planar() -> float:
    """Finite generating set: e_□ = (12/11) G(r/11) with continuous r.

    Clock-unit rival G(r/12) is the inverse-direction tail. Do not use T/N.
    """
    return leftover_exp_finite()


def leftover_octave(z_star: float) -> dict:
    """Integer octaves live in a(t). Leftover is the 12-fold remainder of that address.

    Binary address: log2(1+z) = n + φ, φ ∈ [0,1).
    12-fold reading: s = 12 φ steps of the period (not a χ² rounding).
    Then f = 2^{(s/12) e} = 2^{φ e}. Continuous density; discrete leftover lattice.
    """
    zs = float(z_star)
    octaves = math.log2(1.0 + zs)
    n_int = int(math.floor(octaves))
    frac = octaves - n_int
    steps_12 = frac * PERIOD
    f_clock = (1.0 + zs) / (2.0 ** n_int)
    f_mixed = 2.0 ** ((steps_12 / PERIOD) * DIM_EXP_2)
    f_planar = 2.0 ** ((steps_12 / PERIOD) * leftover_exp_planar())
    f_space = 2.0 ** ((steps_12 / PERIOD) * DIM_EXP_3)
    f_st = 2.0 ** ((steps_12 / PERIOD) * DIM_EXP_4)
    f_finite4 = 2.0 ** ((steps_12 / PERIOD) * DIM_EXP_MEAN_4)
    f_inf = 2.0 ** ((steps_12 / PERIOD) * leftover_exp_infinite())
    return {
        "z_star": zs,
        "octaves": octaves,
        "n_int": n_int,
        "frac": frac,
        "steps_12": steps_12,
        "f_exact": f_clock,
        "f_mixed": f_mixed,
        "f_planar": f_planar,
        "f_space": f_space,
        "f_spacetime": f_st,
        "f_finite4": f_finite4,
        "f_inf": f_inf,
    }


def leftover_f(z_star: float) -> float:
    """Theorem: f = 2^{ {log2(1+z*)} · [1 − (1−r) ln(1−r)] }.

    Infinite Euclidean dimensions, tail-weighted. 4D cutoff is not this theorem.
    z* is photon decoupling. Does not partition T.
    """
    return leftover_octave(z_star)["f_inf"]


def remainder_k(B: int, frac: float) -> int:
    return int(round(B * frac))


def remainder_f(B: int, frac: float) -> float:
    k = remainder_k(B, frac)
    return 2.0 ** (k / B) if B else 1.0


def mix_stretch_f(frac: float) -> float:
    """Clock leftover referred to the 11 hyperbolic generators: 2^{frac · 12/11}."""
    return 2.0 ** (frac * PERIOD / Q5_WEIGHT)


def mix_compress_f(frac: float) -> float:
    """Inverse measure 11/12. Not the theorem: φ is not in generator units."""
    return 2.0 ** (frac * Q5_WEIGHT / PERIOD)


def continuum_equal_f(frac: float) -> float:
    """B→∞ equal-tempered leftover. The remainder family's limit."""
    return 2.0 ** frac


def continuum_mixed_f(frac: float) -> float:
    """B→∞ leftover with the 1:11 ratio kept on the leftover itself."""
    return mix_stretch_f(frac)


def tail_per_octave_f(z_star: float | None = None) -> float:
    """Ansatz leftover f = 1/(1 − T/N_int). N is not a continuum theorem."""
    n = leftover_octave(z_star)["n_int"] if z_star is not None else 10
    return 1.0 / (1.0 - float(OMEGA_DE_TODAY) / n)


def cliff_frac(B: int, k: int) -> tuple[float, float]:
    """frac in [k - 1/2, k + 1/2) / B rounds to k."""
    lo = (k - 0.5) / B
    hi = (k + 0.5) / B
    return lo, hi


def z_from_frac_above_n(n_int: int, frac: float) -> float:
    return 2.0 ** (n_int + frac) - 1.0

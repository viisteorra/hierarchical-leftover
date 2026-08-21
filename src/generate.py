"""One generating function. Density and leftover are evaluations of G.

G(ρ) = 1/(1−ρ)

    T      = G(r) − 1                         Theorem 1
    e_□    = (12/11) G(r/11)                  finite generating set (11 steps)
    e_∞    = 1 + (1−r) ln G(r)                infinite generating set (dimensions)
    e_foil = (12/11) G(r/12)                  clock-unit; wrong direction (Lemma 3)

No third leftover. U1–U3 license two ratios only: r (the mix) and r/11
(leftover on the generators). r/12 is the inverse-direction foil.
The Hubble constant does not enter. Likelihood does not enter.
"""

from __future__ import annotations

import math

from geometry import Q5_WEIGHT, TOTAL, r as R_GEN

PERIOD = int(TOTAL)  # 12
GENERATORS = int(Q5_WEIGHT)  # 11


def G(rho: float) -> float:
    """Geometric tail 1/(1−ρ). This is Theorem 1's series, not a new knob.

    Σ_{k=0}^∞ ρ^k = 1/(1−ρ). Density is that series at the mix: T = G(r)−1.
    """
    rho = float(rho)
    if rho >= 1.0:
        raise ValueError("G(ρ) needs ρ < 1")
    return 1.0 / (1.0 - rho)


def spatial_curvature() -> int:
    """FRW k = 0.

    Metric leftover is Euclidean hypercubes q=2d (Theorem 2). Hyperbolic q=5
    lives in the 2D vacuum mix (density), not in the 3-geometry. Flat is a
    lemma from U2 + leftover, not a fit of Ω_k.
    """
    return 0


def density_tail() -> float:
    """T = G(r) − 1 = r/(1−r)."""
    return G(R_GEN) - 1.0


def leftover_exp_finite() -> float:
    """Planar leftover exponent: clock referred to 11 generators, then G(r/11)."""
    return (PERIOD / GENERATORS) * G(R_GEN / GENERATORS)


def leftover_exp_infinite() -> float:
    """Infinite-D leftover exponent: 1 + (1−r) ln G(r) = 1 − (1−r) ln(1−r)."""
    return 1.0 + (1.0 - R_GEN) * math.log(G(R_GEN))


def leftover_exp_clock_foil() -> float:
    """r/12 tail. Lemma 3 forbids it: leftover lives on the 11, not the clock."""
    return (PERIOD / GENERATORS) * G(R_GEN / PERIOD)


def licensed_ratios() -> dict[str, float]:
    """The only ρ U1–U3 put in G for physics. Not a scan."""
    return {
        "mix": float(R_GEN),  # density
        "generators": float(R_GEN) / GENERATORS,  # planar leftover
    }


def generating_set(observable: str) -> str:
    """Lemma C, simplified: finite set → e_□, infinite set → e_∞, none → f=1.

    Spatial slice rulers couple to the finite 11 generators.
    Null angles in the metric couple to the unbounded dimension tail.
    This octave (SN, time delays) has no leftover.
    Early microphysics (BBN, recombination internals) is Einstein-frame.
    """
    key = (
        observable.strip()
        .lower()
        .replace("*", "_star")
        .replace("θ", "theta")
        .replace(" ", "_")
    )
    finite = {
        "bao",
        "dm",
        "dh",
        "dv",
        "spatial",
        "slice",
        "d_m",
        "d_h",
        "d_v",
        "weak_lensing",
        "galaxy_shape",
        "clustering",
    }
    infinite = {
        "theta_star",
        "thetastar",
        "cmb",
        "null",
        "metric",
        "theta",
    }
    none = {
        "sn",
        "sn_local",
        "time_delay",
        "td",
        "h0_td",
        "local",
    }
    early = {
        "bbn",
        "cmb_r",
        "shift_r",
        "recombination",
        "omega_b",
    }
    if key in finite:
        return "finite"
    if key in infinite:
        return "infinite"
    if key in none:
        return "none"
    if key in early:
        return "early"
    raise ValueError(f"unknown probe: {observable!r}")

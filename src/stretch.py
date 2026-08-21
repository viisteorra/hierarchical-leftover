"""Local SN-only distance stretch. Does not touch 49/71, BAO, or CMB.

S(z) = 1 + δ0 / (1+z)

This is one redshift octave of fade: 2^{-log2(1+z)} = 1/(1+z).
At z=0, S=1+δ0. At high z, S→1. BAO and CMB never call this function.

δ0 is phenomenological. Stage A of the cohesion test sets it to zero.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np

from geometry import OMEGA_DE_TODAY, TOTAL

# Locked geometric density. Do not change.
OMEGA_DE = Fraction(49, 71)

# Octave period — same 12 as the 1:11 mix.
B = int(TOTAL)

assert abs(float(OMEGA_DE) - OMEGA_DE_TODAY) < 1e-15
assert B == 12


def stretch_S(z, delta0: float):
    """Multiplicative stretch on SN luminosity distance."""
    z = np.asarray(z, dtype=float)
    if delta0 == 0.0:
        return np.ones_like(z, dtype=float) if z.shape else 1.0
    return 1.0 + float(delta0) / (1.0 + z)


def stretch_dmu(z, delta0: float):
    """Additive stretch on the distance modulus: 5 log10(S)."""
    s = np.asarray(stretch_S(z, delta0), dtype=float)
    if np.any(s <= 0.0):
        raise ValueError("stretch_S must stay positive")
    return 5.0 * np.log10(s)

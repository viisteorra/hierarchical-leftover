"""Local SN-only distance stretch. Does not touch ln(2), BAO, or CMB.

S(z) = 1 + δ0 / (1+z)

This is one redshift octave of fade: 2^{-log2(1+z)} = 1/(1+z).
At z=0, S=1+δ0. At high z, S→1. BAO and CMB never call this function.

δ0 is phenomenological. Stage A of the cohesion test sets it to zero.
"""

from __future__ import annotations

import numpy as np

from geometry import OMEGA_DE_TODAY, TOTAL

B = int(TOTAL)  # leftover 12-fold addressing
assert B == 12
assert abs(OMEGA_DE_TODAY - __import__("math").log(2.0)) < 1e-15


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

"""Continuous binary-scale residual. Not a fit parameter.

The vacuum resonates on the octave 2:1. In continuous log-space the
Haar measure of one doubling is ln(2). That residual is Ω_DE today:

    Ω_DE = ln(2),   r = ln(2)/(1+ln(2)),   T = r/(1-r) = ln(2).

Density is ln(2). The old rationals are the same quantities:

    49/71  → ln(2)           (T, Ω_DE)
    49/120 → ln(2)/(1+ln(2)) (r)
    71/120 → 1/(1+ln(2))     (1−r)

12 and 12/11 remain leftover *addressing* (period / generators), not
those rationals. Do not put H0 into leftover f.
"""

from __future__ import annotations

import math

LN2 = math.log(2.0)
OMEGA_DE_TODAY = LN2
r = LN2 / (1.0 + LN2)
OMEGA_M = 1.0 - OMEGA_DE_TODAY
W = -1.0  # period does not drift
K_CURVATURE = 0  # Euclidean metric leftover

# Leftover addressing lattice only (period 12, 11 generators). Not Ω_DE.
Q4_WEIGHT = 1
Q5_WEIGHT = 11
TOTAL = Q4_WEIGHT + Q5_WEIGHT  # 12

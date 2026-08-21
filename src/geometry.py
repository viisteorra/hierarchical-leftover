"""Continuous binary-scale residual. Not a fit parameter.

The vacuum resonates on the octave 2:1. In continuous log-space the
Haar measure of one doubling is ln(2). That residual is Ω_DE today:

    Ω_DE = ln(2),   r = ln(2)/(1+ln(2)),   T = r/(1-r) = ln(2).

Density is continuous. Leftover is still 12-fold: the resonator is
*addressed* in equal 12-fold steps of the binary octave (U1–U3). The
1:11 mix (r=49/120, T=49/71) is that lattice, not a χ² retune of ln(2).
Planar leftover lives on the 11 generators of that lattice. Infinite
leftover uses the continuous mix r. Do not put H0 into leftover f.
"""

from __future__ import annotations

import math

LN2 = math.log(2.0)
OMEGA_DE_TODAY = LN2
r = LN2 / (1.0 + LN2)
OMEGA_M = 1.0 - OMEGA_DE_TODAY
W = -1.0  # period does not drift
K_CURVATURE = 0  # Euclidean metric leftover

# 12-fold leftover lattice (U1–U3). Density lock is ln(2), not this T.
Q4_WEIGHT = 1
Q5_WEIGHT = 11
TOTAL = Q4_WEIGHT + Q5_WEIGHT  # 12
R_DISCRETE = (Q4_WEIGHT * (2 / 4) + Q5_WEIGHT * (2 / 5)) / TOTAL  # 49/120
T_DISCRETE = R_DISCRETE / (1.0 - R_DISCRETE)  # 49/71

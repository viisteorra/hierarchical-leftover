"""Locked 1:11 mix on a 12-fold period. Not a fit parameter.

Uniqueness (src/uniqueness.py, docs/AXIOMS.md):
  U1  12 = denominator of the first CF convergent of log2(3/2) after 3/5.
  U2  q=4 Euclidean hypercube in 2D; q=5 first hyperbolic square tiling.
  U3  exactly one Euclidean seed per period ⇒ 1:11.
Then r=49/120, T=49/71. Do not retune from χ².
"""

# 1 part q=4  (r = 2/4 = 0.5)  +  11 parts q=5  (r = 2/5 = 0.4)
Q4_WEIGHT = 1
Q5_WEIGHT = 11
TOTAL = Q4_WEIGHT + Q5_WEIGHT          # 12

r = (Q4_WEIGHT * (2 / 4) + Q5_WEIGHT * (2 / 5)) / TOTAL   # 4.9/12
OMEGA_DE_TODAY = r / (1 - r)           # 4.9/7.1 ≈ 0.69014
OMEGA_M = 1.0 - OMEGA_DE_TODAY  # k=0: Euclidean metric leftover, not a fit of Ω_k
W = -1.0  # period does not drift (lemma from P + U1)
K_CURVATURE = 0  # lemma: q=2d Euclidean hypercubes, not hyperbolic 3-geometry

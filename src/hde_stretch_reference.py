"""Frozen identities for stretch unit tests. Do not edit to chase χ²."""

from fractions import Fraction

OMEGA_DE = Fraction(49, 71)
OMEGA_DE_FLOAT = 49.0 / 71.0
B = 12

# z, delta0, S, dmu = 5 log10(S)
# S = 1 + delta0/(1+z)
REFERENCE_STRETCH = (
    (0.0, 0.0, 1.0, 0.0),
    (1.0, 0.0, 1.0, 0.0),
    (3.0, 0.0, 1.0, 0.0),
    (0.0, 0.02, 1.02, 0.043000858809588),
    (1.0, 0.02, 1.01, 0.021606868913213),
    (3.0, 0.02, 1.005, 0.010830308782538),
)

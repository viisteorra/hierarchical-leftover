"""Single continuum. One metric. Same job Einstein’s metric does.

A0  Analogous axioms A–D and uniqueness U1–U3: docs/AXIOMS.md.
A1  Density: Ω_DE = 49/71, w = −1, flat. geometry.py. 2D vacuum mix.

A2  Leftover — one generating function G(ρ)=1/(1−ρ), two evaluations:
        T   = G(r)−1
        e_□ = (12/11) G(r/11)     finite set (BAO / spatial slices)
        e_∞ = 1+(1-r) ln G(r) = 1-(1-r)ln(1-r)   infinite set (θ* / metric)
        f   = 2^{φ e}
    No third leftover (r/12 is the Lemma-3 foil). Do not mix on one pair.
    Proof: docs/THEOREM.md, src/generate.py.  z* is recombination, not H0.

A3  One FRW metric: ds² = −c² dt² + a(t)² dχ². E(z) from A1.

A4  Two readings, measure matches the observable (Theorem 3):
        BAO:   H0_L = f_□ × H0_E(BAO+rd) ,  rd_L = rd_E / f_□
        CMB θ*: H0_L = f_∞ × H0_E(θ*)
    SN tests both. Early–early: H0_E(BAO)/H0_E(θ*) = f_∞/f_□.
    Hybrid (f_□ H0 with rd/f_∞) is mixed frames.

A5  Hierarchical time: t_U = t_FRW(H0_L) × f (f of that reading).

A6  Time delays in this octave: H0_TD = H0_L (SN).
"""

from __future__ import annotations

from generate import generating_set
from ruler import leftover_f, leftover_octave, tail_per_octave_f

# Photon-decoupling redshift: recombination measurement / CAMB at lock.
# Not SH0ES. Not chosen to make f = 1.074.
ZSTAR_RECOMB = 1089.84  # CAMB at lock Ω; Planck catalog 1089.80
ZSTAR_HS_LOCK = 1092.0806837280463
PLANCK_ZSTAR = 1089.80
HS_BIAS = 2.12
ZSTAR_CAMB = 1089.8389489807696

F_BAO = leftover_octave(ZSTAR_RECOMB)["f_planar"]  # tailed planar, BAO
F_CMB = leftover_f(ZSTAR_RECOMB)  # infinite Euclidean spacetime
F_AXIOM = F_CMB  # spacetime default; BAO must use F_BAO
F_CLOCK = leftover_octave(ZSTAR_RECOMB)["f_exact"]
F_MIXED_2D = leftover_octave(ZSTAR_RECOMB)["f_mixed"]  # bare 12/11 lemma
F_FINITE4 = leftover_octave(ZSTAR_RECOMB)["f_finite4"]
F_HS = leftover_f(ZSTAR_HS_LOCK)
F_CORR = leftover_f(ZSTAR_HS_LOCK - HS_BIAS)
F_PLANCK_Z = leftover_f(PLANCK_ZSTAR)
F_CAMB = leftover_f(ZSTAR_CAMB)
F_ANSATZ_T10 = tail_per_octave_f()


def measure_for(observable: str) -> str:
    """Lemma C. Finite generating set → planar; infinite set → f_∞.

    Same G. Not chosen from χ². SN / time delays have no leftover (this octave).
    """
    kind = generating_set(observable)
    if kind == "finite":
        return "planar"
    if kind == "infinite":
        return "infinite"
    raise ValueError(
        f"{observable!r} generating set is {kind!r}, not a leftover measure"
    )


def f_for(observable: str, z_star: float | None = None) -> float:
    """Leftover f for an observable. Measure is Lemma C, not a fit."""
    zs = ZSTAR_RECOMB if z_star is None else float(z_star)
    kind = measure_for(observable)
    if kind == "planar":
        return leftover_octave(zs)["f_planar"]
    return leftover_f(zs)


def f_leftover(z_star: float) -> float:
    """Spacetime leftover f_∞ at a given recombination z*."""
    return leftover_f(z_star)


def f_bao(z_star: float) -> float:
    """Tailed planar leftover f_□ at a given recombination z*."""
    return leftover_octave(z_star)["f_planar"]

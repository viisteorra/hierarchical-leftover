#!/usr/bin/env python3
"""Double-check after every enhancement. Exit 0 only if the geometric core holds."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def fail(msg: str) -> None:
    print(f"VERIFY FAIL: {msg}")
    sys.exit(1)


def main() -> int:
    from geometry import K_CURVATURE, OMEGA_DE_TODAY, Q4_WEIGHT, Q5_WEIGHT, TOTAL, W, r
    from uniqueness import PERIOD_N, SEED, HYPER_STEPS, T_FRAC, period_from_fifth_convergent, smallest_nfold_fifth
    from generate import G, density_tail, leftover_exp_finite, leftover_exp_infinite as g_inf
    from ruler import leftover_exp_infinite, leftover_f, leftover_octave, mix_compress_f
    from spacetime import F_ANSATZ_T10, F_AXIOM, F_BAO, F_CLOCK, F_CMB, F_FINITE4, F_MIXED_2D, ZSTAR_RECOMB

    geo = (ROOT / "src" / "geometry.py").read_text()
    if "Q4_WEIGHT = 1" not in geo or "Q5_WEIGHT = 11" not in geo:
        fail("1:11 weights changed in geometry.py")
    if "chi2" in geo.lower() or "minimize" in geo:
        fail("geometry.py must not contain a fit")
    import math as _math
    if abs(OMEGA_DE_TODAY - _math.log(2.0)) >= 1e-15:
        fail(f"Ω_DE drifted from ln(2): {OMEGA_DE_TODAY}")
    if W != -1.0:
        fail("w drifted from -1")
    if K_CURVATURE != 0:
        fail("spatial curvature drifted from k=0")
    if not (ROOT / "docs" / "SOLVED.md").exists():
        fail("SOLVED.md missing")
    if abs(r - _math.log(2.0) / (1.0 + _math.log(2.0))) >= 1e-15:
        fail(f"r drifted: {r}")
    if TOTAL != 12 or Q4_WEIGHT / Q5_WEIGHT != 1 / 11:
        fail("period/mix drifted")
    n, k, _c = smallest_nfold_fifth()
    if (n, k) != (12, 7) or PERIOD_N != 12:
        fail("U1: 12-fold fifth uniqueness drifted")
    n_cf, k_cf, _ = period_from_fifth_convergent()
    if (n_cf, k_cf) != (12, 7):
        fail("U1: CF period drifted from 7/12")
    if (SEED, HYPER_STEPS) != (1, 11):
        fail("U3: one-seed fill drifted")
    if abs(float(T_FRAC) - 49 / 71) >= 1e-15:
        fail("U3 tail drifted")
    if abs(density_tail() - _math.log(2.0)) >= 1e-15:
        fail("G(r)-1 drifted from ln(2)")
    if abs(leftover_exp_finite() - (12 / 11) * G(r / 11)) >= 1e-15:
        fail("G planar leftover drifted from (12/11)G(r/11)")
    if abs(g_inf() - leftover_exp_infinite()) >= 1e-15:
        fail("G infinite leftover drifted from ruler")
    if abs(g_inf() - (1.0 + (1.0 - r) * __import__("math").log(G(r)))) >= 1e-15:
        fail("e_∞ is not 1+(1-r)ln G(r)")
    ax = (ROOT / "docs" / "AXIOMS.md").read_text()
    if "7/12" not in ax or "exactly one" not in ax.lower():
        fail("AXIOMS.md missing uniqueness lemmas")
    if "one primitive" not in ax.lower() and "The one primitive" not in ax:
        fail("AXIOMS.md missing the single primitive P")

    lo = leftover_octave(ZSTAR_RECOMB)
    e_inf = leftover_exp_infinite()
    if abs(e_inf - (1.0 - (1.0 - r) * __import__("math").log(1.0 - r))) >= 1e-15:
        fail("infinite leftover exponent drifted from 1-(1-r)ln(1-r)")
    if abs(F_CLOCK - 2.0 ** lo["frac"]) >= 1e-12:
        fail("F_CLOCK is not 2^{frac}")
    if abs(F_CLOCK - 1090.84 / 1024) >= 1e-15:
        fail("clock leftover drifted from 1090.84/1024")
    if abs(F_AXIOM - 2.0 ** (lo["frac"] * e_inf)) >= 1e-15:
        fail("F_AXIOM is not 2^{frac · ⟨e⟩_∞}")
    if abs(F_AXIOM - leftover_f(ZSTAR_RECOMB)) >= 1e-15:
        fail("F_AXIOM != leftover_f(z*)")
    if abs(F_CMB - leftover_f(ZSTAR_RECOMB)) >= 1e-15:
        fail("F_CMB is not infinite leftover")
    from ruler import leftover_exp_planar
    if abs(leftover_exp_planar() - leftover_exp_finite()) >= 1e-15:
        fail("planar leftover exponent drifted from generate")
    if abs(F_BAO - 2.0 ** (lo["frac"] * leftover_exp_planar())) >= 1e-15:
        fail("F_BAO is not tailed planar leftover")
    if abs(F_BAO - F_MIXED_2D) <= 1e-4:
        fail("F_BAO collapsed onto bare 12/11")
    if abs(F_BAO - F_CMB) <= 1e-4:
        fail("BAO and CMB leftover collapsed")
    if abs(F_AXIOM - mix_compress_f(lo["frac"])) <= 1e-4:
        fail("theorem collapsed onto inverse 11/12")
    if abs(F_AXIOM - F_CLOCK) <= 1e-4:
        fail("infinite leftover collapsed onto clock")
    if abs(F_AXIOM - F_MIXED_2D) <= 1e-4:
        fail("infinite leftover collapsed onto 2D mix")
    if abs(F_AXIOM - F_FINITE4) <= 1e-4:
        fail("infinite leftover collapsed onto 4D cutoff")
    if abs(F_ANSATZ_T10 - 1.0 / (1.0 - OMEGA_DE_TODAY / 10.0)) >= 1e-12:
        fail("ansatz T/10 drifted")
    if abs(F_AXIOM - F_ANSATZ_T10) <= 1e-4:
        fail("theorem leftover collapsed onto T/10 ansatz")
    th = (ROOT / "docs" / "THEOREM.md").read_text()
    if "ln(1-r)" not in th and r"\ln(1-r)" not in th:
        fail("THEOREM.md missing infinite-D leftover statement")
    if "Theorem 3" not in th:
        fail("THEOREM.md missing Theorem 3 (Hubble scale)")
    if "One generating function" not in th:
        fail("THEOREM.md missing generating-function simplification")
    if not (ROOT / "docs" / "THEORY.md").exists():
        fail("THEORY.md missing")
    if "Closed" not in (ROOT / "docs" / "THEORY.md").read_text():
        fail("THEORY.md missing Closed section")
    if not (ROOT / "docs" / "COMPLETE.md").exists():
        fail("COMPLETE.md missing")
    if not (ROOT / "src" / "generate.py").exists():
        fail("src/generate.py missing")
    if not (ROOT / "src" / "hubble.py").exists():
        fail("src/hubble.py missing")
    ruler = (ROOT / "src" / "ruler.py").read_text()
    if "73.47" in ruler or "H0_SN" in ruler:
        fail("ruler.py must not contain H0")
    from hubble import solution as hubble_solution
    hs = hubble_solution()
    if hs["map_bao_sigma"] >= 3.0 or hs["map_th_sigma"] >= 3.0:
        fail("Hubble maps reopened")
    if abs(hs["f_box"] - F_BAO) >= 1e-15 or abs(hs["f_inf"] - F_CMB) >= 1e-15:
        fail("hubble.solution leftover drifted from spacetime")

    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", str(ROOT / "tests")],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr)
        fail("pytest failed")
    print("VERIFY OK", proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else "")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

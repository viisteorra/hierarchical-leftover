#!/usr/bin/env python3
"""Empirical score of the continuous ln(2) lock. No free cosmological parameters."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bao import DESI_DR2, PLANCK_RD_MPC, bao_chi2, fit_h0_given_rd  # noqa: E402
from cosmology import shift_parameter_R  # noqa: E402
from geometry import LN2, OMEGA_DE_TODAY, OMEGA_M, R_DISCRETE, T_DISCRETE, W, r as R_GEN  # noqa: E402
from hubble import SIG_BAO, SIG_PLANCK_H0, SIG_SN, solution, tension  # noqa: E402
from likelihood import (  # noqa: E402
    PLANCK_R,
    PLANCK_R_SIGMA,
    SNLikelihood,
    load_pantheon_plus,
)
from spacetime import F_BAO, F_CMB, ZSTAR_RECOMB  # noqa: E402

RESULTS = ROOT / "results"
DATA = ROOT / "data"


def main() -> int:
    like = SNLikelihood(
        load_pantheon_plus(DATA / "Pantheon+SH0ES.dat", DATA / "Pantheon+SH0ES_STAT+SYS.cov")
    )
    sn = like.fit_h0(omega_de=OMEGA_DE_TODAY)
    sn_free = like.fit_h0_omega_de(H0_guess=sn["H0"], omega_guess=OMEGA_DE_TODAY)
    bao = fit_h0_given_rd(PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    s = solution(h0_sn=float(sn["H0"]), h0_bao=float(bao["H0"]))
    dchi_sn = sn["chi2"] - sn_free["chi2"]
    dchi_bao = bao_chi2(sn["H0"], PLANCK_RD_MPC / F_BAO, sample=DESI_DR2) - bao["chi2"]
    R = shift_parameter_R(sn["H0"] / F_CMB, omega_de=OMEGA_DE_TODAY)
    R_pull = (R - PLANCK_R) / PLANCK_R_SIGMA

    print("=" * 72)
    print("CONTINUOUS LOCK  Ω_DE = ln(2)  w=-1  k=0  official DESI DR2 cov")
    print("=" * 72)
    print(f"  Ω_DE={OMEGA_DE_TODAY:.12f}  Ω_m={OMEGA_M:.12f}  r={R_GEN:.12f}")
    print(f"  discrete 12-fold T={T_DISCRETE:.12f}=49/71  Δ={T_DISCRETE-LN2:+.6f}")
    print(f"  discrete r={R_DISCRETE:.12f}=49/120")
    print(f"  z*={ZSTAR_RECOMB}  f_□={F_BAO:.6f}  f_∞={F_CMB:.6f}")
    print(f"  w={W}")
    print()
    print("--- SN Pantheon+ SH0ES ---")
    print(f"  lock H0={sn['H0']:.4f}  χ²={sn['chi2']:.2f}  n={sn['n']}")
    print(f"  free Ω={sn_free['omega_de']:.4f}  χ²={sn_free['chi2']:.2f}  Δχ²={dchi_sn:.2f}")
    print()
    print("--- DESI DR2 official 13×13 cov + Planck rd ---")
    print(f"  H0={bao['H0']:.4f}  χ²={bao['chi2']:.3f}  n={bao['n']}  dof={bao['dof']}")
    print(f"  leftover pair (H0_SN, rd/f_□) Δχ²={dchi_bao:+.3f}")
    print()
    print("--- leftover maps ---")
    print(
        f"  BAO  {s['h0_bao']:.3f} × f_□ = {s['pred_sn_from_bao']:.3f}  vs SN {s['h0_sn']:.3f}  "
        f"{s['map_bao_sigma']:.3f}σ  raw {s['raw_bao_sigma']:.2f}σ"
    )
    print(
        f"  θ*   {s['h0_th']:.3f} × f_∞ = {s['pred_sn_from_th']:.3f}  vs SN {s['h0_sn']:.3f}  "
        f"{s['map_th_sigma']:.3f}σ  raw {s['raw_th_sigma']:.2f}σ"
    )
    print(f"  CMB R={R:.4f}  pull {R_pull:+.2f}σ")
    print(f"  SH0ES cal ±{SIG_SN}  BAO ±{SIG_BAO}  Planck H0 ±{SIG_PLANCK_H0}")
    print("=" * 72)

    out = {
        "omega_de": OMEGA_DE_TODAY,
        "omega_m": OMEGA_M,
        "r": R_GEN,
        "T_discrete": T_DISCRETE,
        "f_box": F_BAO,
        "f_inf": F_CMB,
        "sn": sn,
        "sn_free": sn_free,
        "bao": bao,
        "dchi_sn": dchi_sn,
        "dchi_bao_k2": dchi_bao,
        "R": R,
        "R_pull": R_pull,
        "solution": s,
    }
    path = RESULTS / "continuous.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

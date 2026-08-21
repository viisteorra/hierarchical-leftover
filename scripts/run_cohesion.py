#!/usr/bin/env python3
"""Cohesion campaign. Do not change 49/71. Do not stretch BAO or CMB.

Stage A — lock Ω_DE, no stretch (the original claim).
Stage B — one SN-only parameter δ0 (local octave layer, phenomenological).
Stage C — free Ω_DE comparison (Δχ² of the lock). Does not free w(z).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from bao import DESI_DR2, fit_rd  # noqa: E402
from cosmology import shift_parameter_R  # noqa: E402
from geometry import OMEGA_DE_TODAY, Q4_WEIGHT, Q5_WEIGHT  # noqa: E402
from likelihood import (  # noqa: E402
    PLANCK2018_H0,
    PLANCK2018_OMEGA_DE,
    PLANCK_R,
    PLANCK_R_SIGMA,
    SNLikelihood,
    load_des_dovekie,
    load_pantheon_plus,
)
from stretch import B, OMEGA_DE  # noqa: E402

DATA = ROOT / "data"
RESULTS = ROOT / "results"
PANTHEON_DAT = DATA / "Pantheon+SH0ES.dat"
PANTHEON_COV = DATA / "Pantheon+SH0ES_STAT+SYS.cov"
DES_HD = DATA / "DES-Dovekie_HD.csv"
DES_NPZ = DATA / "DES-Dovekie_STAT+SYS.npz"


def _R_chi2(omega_de: float) -> float:
    R = shift_parameter_R(PLANCK2018_H0, omega_de=omega_de, w=-1.0)
    return float(((R - PLANCK_R) / PLANCK_R_SIGMA) ** 2)


def main() -> int:
    if not (PANTHEON_DAT.exists() and DES_HD.exists()):
        raise SystemExit("missing data; run python scripts/download_data.py")
    RESULTS.mkdir(parents=True, exist_ok=True)

    print("=" * 64)
    print("Cohesion — Stages A, B, C")
    print("=" * 64)
    print(f"  lock Ω_DE = {float(OMEGA_DE):.12f}  = 49/71")
    print(f"  B (octave period) = {B}")
    print(f"  mix = {Q4_WEIGHT}:{Q5_WEIGHT}")
    print("  SN stretch: S=1+δ0/(1+z)   BAO/CMB: never stretched")
    print("  w is locked at -1")
    print()

    like_p = SNLikelihood(load_pantheon_plus(PANTHEON_DAT, PANTHEON_COV))
    like_d = SNLikelihood(load_des_dovekie(DES_HD, DES_NPZ))

    print("Stage A: lock, no stretch …")
    a_p = like_p.fit_h0(w=-1.0, omega_de=OMEGA_DE_TODAY)
    a_d = like_d.fit_h0(w=-1.0, omega_de=OMEGA_DE_TODAY)
    a_bao = fit_rd(70.0, omega_de=OMEGA_DE_TODAY, w=-1.0, sample=DESI_DR2)
    a_R = _R_chi2(OMEGA_DE_TODAY)
    print(f"  Pantheon+  H0={a_p['H0']:.4f}  χ²={a_p['chi2']:.3f}  χ²/dof={a_p['chi2_dof']:.4f}")
    print(f"  DES-SN5YR  H0={a_d['H0']:.4f}  χ²={a_d['chi2']:.3f}  χ²/dof={a_d['chi2_dof']:.4f}")
    print(f"  DESI DR2   rd(h=0.7)={a_bao['rd']:.2f}  χ²={a_bao['chi2']:.3f} / {a_bao['dof']}")
    print(f"  CMB R      χ²={a_R:.3f}")
    print()

    print("Stage B: lock, SN-only δ0 …")
    b_p = like_p.fit_h0_delta0(omega_de=OMEGA_DE_TODAY, H0_guess=a_p["H0"])
    b_d = like_d.fit_h0_delta0(omega_de=OMEGA_DE_TODAY, H0_guess=a_d["H0"])
    print(
        f"  Pantheon+  H0={b_p['H0']:.4f}  δ0={b_p['delta0']:+.5f}  "
        f"χ²={b_p['chi2']:.3f}  Δχ² vs A={b_p['chi2']-a_p['chi2']:+.3f}"
    )
    print(
        f"  DES-SN5YR  H0={b_d['H0']:.4f}  δ0={b_d['delta0']:+.5f}  "
        f"χ²={b_d['chi2']:.3f}  Δχ² vs A={b_d['chi2']-a_d['chi2']:+.3f}"
    )
    print("  BAO and CMB unchanged from Stage A (not stretched)")
    print()

    print("Stage C: free Ω_DE, no stretch, w=-1 …")
    c_p = like_p.fit_h0_omega_de(w=-1.0, H0_guess=a_p["H0"], omega_guess=OMEGA_DE_TODAY)
    c_d = like_d.fit_h0_omega_de(w=-1.0, H0_guess=a_d["H0"], omega_guess=OMEGA_DE_TODAY)
    # BAO: scan Ω_DE, profile rd
    omega_grid = np.linspace(0.62, 0.76, 29)
    chi2_bao = np.array(
        [fit_rd(70.0, omega_de=float(o), w=-1.0, sample=DESI_DR2)["chi2"] for o in omega_grid]
    )
    i_bao = int(np.argmin(chi2_bao))
    c_bao_omega = float(omega_grid[i_bao])
    c_bao_chi2 = float(chi2_bao[i_bao])
    print(
        f"  Pantheon+  Ω_DE={c_p['omega_de']:.4f}  H0={c_p['H0']:.4f}  "
        f"χ²={c_p['chi2']:.3f}  lock Δχ²={a_p['chi2']-c_p['chi2']:+.3f}"
    )
    print(
        f"  DES-SN5YR  Ω_DE={c_d['omega_de']:.4f}  H0={c_d['H0']:.4f}  "
        f"χ²={c_d['chi2']:.3f}  lock Δχ²={a_d['chi2']-c_d['chi2']:+.3f}"
    )
    print(
        f"  DESI DR2   Ω_DE={c_bao_omega:.4f}  χ²={c_bao_chi2:.3f}  "
        f"lock Δχ²={a_bao['chi2']-c_bao_chi2:+.3f}"
    )
    print()
    print("Not claimed: Hubble tension is not solved. w(z) was not freed.")
    print("Stage B stretch is SN-only and phenomenological.")

    payload = {
        "lock": float(OMEGA_DE),
        "B": B,
        "stage_A": {
            "pantheon": a_p,
            "des": a_d,
            "bao": a_bao,
            "cmb_R_chi2": a_R,
        },
        "stage_B": {"pantheon": b_p, "des": b_d},
        "stage_C": {
            "pantheon": c_p,
            "des": c_d,
            "bao_omega_de": c_bao_omega,
            "bao_chi2": c_bao_chi2,
            "bao_delta_chi2_lock": a_bao["chi2"] - c_bao_chi2,
        },
    }
    out = RESULTS / "cohesion.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

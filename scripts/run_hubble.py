#!/usr/bin/env python3
"""Finish the Hubble scale: leftover predicts early H0, CAMB scores θ*.

f is not fitted to H0. Inverse maps and the forward CAMB prediction
are both printed. Early–early H0 ratio is the two-measure corollary.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bao import DESI_DR2, PLANCK_RD_MPC, bao_chi2, fit_h0_given_rd  # noqa: E402
from cosmology import age_gyr  # noqa: E402
from geometry import OMEGA_DE_TODAY, OMEGA_M  # noqa: E402
from hubble import (  # noqa: E402
    PLANCK_100THETA,
    PLANCK_100THETA_SIG,
    PLANCK_OMMH2,
    PLANCK_OMMH2_SIG,
    PLANCK_RD,
    PLANCK_RD_SIG,
    SIG_BAO,
    SIG_SN,
    omega_m_h2,
    solution,
    tension,
)
from likelihood import SNLikelihood, load_pantheon_plus  # noqa: E402
from spacetime import F_BAO, F_CMB  # noqa: E402

RESULTS = ROOT / "results"
DATA = ROOT / "data"
OMBH2 = 0.02237


def camb_at(H0: float) -> dict:
    import camb

    h = H0 / 100.0
    omch2 = OMEGA_M * h**2 - OMBH2
    p = camb.CAMBparams()
    p.set_cosmology(H0=H0, ombh2=OMBH2, omch2=omch2, omk=0.0, tau=0.0544)
    p.InitPower.set_params(As=1e-10 * np.exp(3.044), ns=0.9649)
    p.set_for_lmax(2500, lens_potential_accuracy=1)
    r = camb.get_results(p)
    d = r.get_derived_params()
    cls = r.get_lensed_scalar_cls(lmax=2500)
    tt = cls[:2501, 0]
    ell = np.arange(tt.size)
    peak = int(ell[80:400][np.argmax(tt[80:400])])
    ts = float(d["thetastar"])
    th = 100.0 * ts if ts < 1 else ts
    return {
        "H0": H0,
        "100theta": th,
        "rd": float(d["rdrag"]),
        "zstar": float(d["zstar"]),
        "age": float(d["age"]),
        "peak_l": peak,
        "theta_pull": (th - PLANCK_100THETA) / PLANCK_100THETA_SIG,
        "rd_pull": (float(d["rdrag"]) - PLANCK_RD) / PLANCK_RD_SIG,
        "zstar_pull": (float(d["zstar"]) - 1089.80) / 0.21,
        "wm": float(OMEGA_M) * h**2,
        "wm_pull": (float(OMEGA_M) * h**2 - PLANCK_OMMH2) / PLANCK_OMMH2_SIG,
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    like = SNLikelihood(
        load_pantheon_plus(DATA / "Pantheon+SH0ES.dat", DATA / "Pantheon+SH0ES_STAT+SYS.cov")
    )
    sn = like.fit_h0(omega_de=OMEGA_DE_TODAY)
    bao = fit_h0_given_rd(PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    h0_sn = float(sn["H0"])
    h0_bao = float(bao["H0"])
    sol = solution(h0_sn=h0_sn, h0_bao=h0_bao)
    chi0 = float(bao["chi2"])
    chi_map = bao_chi2(h0_sn, PLANCK_RD_MPC / F_BAO, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    chi_hyb = bao_chi2(h0_sn, PLANCK_RD_MPC / F_CMB, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)

    print("=" * 72)
    print("HUBBLE SCALE  given P  (f not fitted to H0)")
    print("=" * 72)
    print(f"  Ω_DE={OMEGA_DE_TODAY:.6f}=49/71  w=-1")
    print(f"  f_□={sol['f_box']:.6f}  f_∞={sol['f_inf']:.6f}  f_∞/f_□={sol['ratio_f']:.6f}")
    print(f"  SN H0={h0_sn:.3f}  DESI+rd={h0_bao:.3f}  θ*-lock={sol['h0_th']:.4f}")

    print("\n--- raw split (no leftover) ---")
    print(
        f"  SN vs BAO   {sol['raw_split_bao_kms']:+.2f} km/s  {sol['raw_bao_sigma']:.2f}σ"
    )
    print(
        f"  SN vs θ*-H0 {sol['raw_split_th_kms']:+.2f} km/s  {sol['raw_th_sigma']:.2f}σ"
    )

    print("\n--- inverse maps (early → local) ---")
    print(
        f"  BAO  {h0_bao:.3f} × f_□ = {sol['pred_sn_from_bao']:.3f}  vs SN {h0_sn:.3f}  "
        f"{sol['map_bao_sigma']:.3f}σ  Δ={sol['split_bao_kms']:+.3f}"
    )
    print(
        f"  θ*   {sol['h0_th']:.3f} × f_∞ = {sol['pred_sn_from_th']:.3f}  vs SN {h0_sn:.3f}  "
        f"{sol['map_th_sigma']:.3f}σ  Δ={sol['split_th_kms']:+.3f}"
    )
    print(f"  BAO χ² early={chi0:.2f}  (H0_SN, rd/f_□)={chi_map:.2f}  Δ={chi_map-chi0:+.3f}")
    print(f"  hybrid (H0_SN, rd/f_∞) χ²={chi_hyb:.2f}  unused")

    print("\n--- forward (local SN → early) ---")
    print(
        f"  SN/f_□ = {sol['pred_bao_from_sn']:.3f}  vs DESI+rd {h0_bao:.3f}  "
        f"{sol['fwd_bao_sigma']:.3f}σ"
    )
    print(f"  SN/f_∞ = {sol['pred_th_from_sn']:.4f}  (CAMB input, not the θ* bisection)")

    print("\n--- corollary: early–early H0 ratio = f_∞/f_□ ---")
    print(
        f"  predicted DESI/θ* = {sol['ratio_f']:.6f}  measured {sol['ratio_h']:.6f}  "
        f"H0_BAO vs θ*×ratio {sol['early_early_sigma']:.3f}σ"
    )

    print("\n--- CAMB at leftover-predicted H0_E = SN/f_∞ ---")
    try:
        pred = camb_at(sol["pred_th_from_sn"])
        lock = camb_at(sol["h0_th"])
        mixed = camb_at(h0_bao)
        camb_ok = True
    except Exception as exc:
        print(f"  CAMB failed: {exc}")
        pred = lock = mixed = {}
        camb_ok = False

    if camb_ok:
        print(
            f"  predicted  H0={pred['H0']:.4f}  100θ*={pred['100theta']:.5f}  "
            f"{pred['theta_pull']:+.2f}σ  rd={pred['rd']:.2f} {pred['rd_pull']:+.2f}σ  "
            f"ℓ={pred['peak_l']}"
        )
        print(
            f"  θ*-lock    H0={lock['H0']:.4f}  100θ*={lock['100theta']:.5f}  "
            f"{lock['theta_pull']:+.2f}σ  (control: bisection, not leftover)"
        )
        print(
            f"  mixed BAO  H0={mixed['H0']:.4f}  100θ*={mixed['100theta']:.5f}  "
            f"{mixed['theta_pull']:+.1f}σ  (wrong frame; leftover does not eat this)"
        )
        wm_cmb = omega_m_h2(sol["pred_th_from_sn"], float(OMEGA_M))
        wm_pull = (wm_cmb - PLANCK_OMMH2) / PLANCK_OMMH2_SIG
        print(
            f"  ω_m CMB-frame={wm_cmb:.5f}  {wm_pull:+.2f}σ  "
            f"CAMB ω_m={pred['wm']:.5f} {pred['wm_pull']:+.2f}σ"
        )
        t_u = age_gyr(h0_sn) * F_CMB
        print(f"  t_U = t_FRW(H0_SN)×f_∞ = {t_u:.3f} Gyr  CAMB age={pred['age']:.3f}")

    print("\n--- Hubble tension residual ---")
    print(f"  BAO scale     {sol['raw_bao_sigma']:.2f}σ → {sol['map_bao_sigma']:.3f}σ")
    if camb_ok:
        print(
            f"  Planck 100θ*  {sol['raw_th_sigma']:.2f}σ (as H0) → "
            f"{pred['theta_pull']:+.2f}σ (forward CAMB)"
        )
        print(f"  inverse θ* map (H0 units, SH0ES ±{SIG_SN}) {sol['map_th_sigma']:.3f}σ")
    print(f"  SH0ES cal ±{SIG_SN}  BAO+rd ±{SIG_BAO}")
    print("=" * 72)

    out = {
        "solution": sol,
        "chi2": {"early": chi0, "mapped": chi_map, "hybrid": chi_hyb, "dchi": chi_map - chi0},
        "camb_predicted": pred,
        "camb_theta_lock": lock,
        "camb_mixed_bao": mixed,
        "camb_ok": camb_ok,
    }
    path = RESULTS / "hubble.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"Wrote {path}")

    # Fail the script only if the scale split reopened.
    bad = sol["map_bao_sigma"] >= 3.0 or sol["map_th_sigma"] >= 3.0
    if camb_ok:
        bad = bad or abs(pred["theta_pull"]) >= 3.0
        bad = bad or pred["peak_l"] not in range(215, 226)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())

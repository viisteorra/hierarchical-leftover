#!/usr/bin/env python3
"""Joint SN + BAO + CMB-R test of the locked Ω_DE = 49/71."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from bao import DESI_DR1, DESI_DR2, fit_rd  # noqa: E402
from cosmology import shift_parameter_R  # noqa: E402
from geometry import OMEGA_DE_TODAY, OMEGA_M, Q4_WEIGHT, Q5_WEIGHT, r  # noqa: E402
from history import HISTORY  # noqa: E402
from likelihood import (  # noqa: E402
    PLANCK2018_H0,
    PLANCK2018_OMEGA_DE,
    PLANCK2018_OMEGA_DE_SIGMA,
    PLANCK_R,
    PLANCK_R_SIGMA,
    SNLikelihood,
    dchi2_interval,
    load_des_dovekie,
    load_pantheon_plus,
)

DATA = ROOT / "data"
RESULTS = ROOT / "results"
PANTHEON_DAT = DATA / "Pantheon+SH0ES.dat"
PANTHEON_COV = DATA / "Pantheon+SH0ES_STAT+SYS.cov"
DES_HD = DATA / "DES-Dovekie_HD.csv"
DES_NPZ = DATA / "DES-Dovekie_STAT+SYS.npz"


def mix_tail(n4: int, n5: int) -> float:
    tot = n4 + n5
    rr = (n4 * 0.5 + n5 * 0.4) / tot
    return rr / (1.0 - rr)


def _sn_profile(like: SNLikelihood, omega: np.ndarray) -> np.ndarray:
    return like.profile(omega, "omega_de", w=-1.0)


def _bao_profile(sample, omega: np.ndarray) -> np.ndarray:
    return np.array([fit_rd(70.0, omega_de=float(o), sample=sample)["chi2"] for o in omega])


def _R_profile(omega: np.ndarray) -> np.ndarray:
    R = np.array([shift_parameter_R(PLANCK2018_H0, omega_de=float(o), w=-1.0) for o in omega])
    return ((R - PLANCK_R) / PLANCK_R_SIGMA) ** 2


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not (PANTHEON_DAT.exists() and DES_HD.exists()):
        raise SystemExit("missing data; run scripts/download_data.py and the DES download")

    print("=" * 72)
    print("Closeout — finish the experiment with data that already exists")
    print("=" * 72)
    print(f"  lock  Ω_DE = {OMEGA_DE_TODAY:.12f}   1:{Q5_WEIGHT}   r={r:.10f}")
    print()

    print("Loading Pantheon+ …")
    like_p = SNLikelihood(load_pantheon_plus(PANTHEON_DAT, PANTHEON_COV))
    print("Loading DES-SN5YR / Dovekie …")
    like_d = SNLikelihood(load_des_dovekie(DES_HD, DES_NPZ))
    print(f"  Pantheon+ N={like_p.data.n}   DES-SN5YR N={like_d.data.n}")

    omega = np.linspace(0.62, 0.76, 29)

    print("Profiling Ω_DE (each probe, H0 or rd maximized) …")
    chi2_p = _sn_profile(like_p, omega)
    chi2_d = _sn_profile(like_d, omega)
    chi2_b1 = _bao_profile(DESI_DR1, omega)
    chi2_b2 = _bao_profile(DESI_DR2, omega)
    chi2_R = _R_profile(omega)

    # Independent-probe combination used for the closeout (not a full MCMC).
    joint = (chi2_p - chi2_p.min()) + (chi2_d - chi2_d.min()) + (chi2_b2 - chi2_b2.min()) + chi2_R
    # chi2_R is already a Δχ² from the Planck prior mean.

    iv_p = dchi2_interval(omega, chi2_p)
    iv_d = dchi2_interval(omega, chi2_d)
    iv_b1 = dchi2_interval(omega, chi2_b1)
    iv_b2 = dchi2_interval(omega, chi2_b2)
    iv_j = dchi2_interval(omega, joint)

    def at_lock(chi2):
        return float(np.interp(OMEGA_DE_TODAY, omega, chi2 - chi2.min()))

    print()
    print("Per-probe Ω_DE (w=-1, scale profiled)")
    print(f"  Pantheon+     {iv_p['best']:.4f} +{iv_p['plus']:.4f}/-{iv_p['minus']:.4f}   lock Δχ²={at_lock(chi2_p):.2f}")
    print(f"  DES-SN5YR     {iv_d['best']:.4f} +{iv_d['plus']:.4f}/-{iv_d['minus']:.4f}   lock Δχ²={at_lock(chi2_d):.2f}")
    print(f"  DESI DR1 BAO  {iv_b1['best']:.4f} +{iv_b1['plus']:.4f}/-{iv_b1['minus']:.4f}   lock Δχ²={at_lock(chi2_b1):.2f}")
    print(f"  DESI DR2 BAO  {iv_b2['best']:.4f} +{iv_b2['plus']:.4f}/-{iv_b2['minus']:.4f}   lock Δχ²={at_lock(chi2_b2):.2f}")
    print(f"  CMB R (Planck H0)  lock Δχ²={at_lock(chi2_R):.2f}")
    print(f"  JOINT (P+ + DES + DR2 + R)  {iv_j['best']:.4f} +{iv_j['plus']:.4f}/-{iv_j['minus']:.4f}")
    print(f"      lock Δχ²={at_lock(joint):.2f}   Planck18 Δχ²={float(np.interp(PLANCK2018_OMEGA_DE, omega, joint-joint.min())):.2f}")

    print()
    print("Historical published Ω_DE vs 4.9/7.1")
    hist_rows = []
    for h in HISTORY:
        pull = (OMEGA_DE_TODAY - h.omega_de) / h.sigma
        hist_rows.append({**h.__dict__, "pull": pull})
        print(f"  {h.year:6.1f}  {h.name:24s}  {h.omega_de:.4f}±{h.sigma:.4f}   pull {pull:+.2f}σ")

    print()
    print("Small integer q=4 / q=5 mixes (comparison only — 1:11 stays locked)")
    mixes = []
    for n4 in range(0, 4):
        for n5 in range(1, 13):
            if n4 == 0 and n5 != 1:
                continue
            T = mix_tail(n4, n5)
            mixes.append((n4, n5, n4 + n5, T, abs(T - OMEGA_DE_TODAY), abs(T - iv_j["best"])))
    mixes.sort(key=lambda t: (t[2], t[4]))
    for n4, n5, tot, T, dlock, djoint in mixes[:12]:
        mark = "  <-- LOCK" if (n4, n5) == (1, 11) else ""
        print(f"  {n4}:{n5:2d}  parts={tot:2d}  T={T:.5f}  |T-lock|={dlock:.5f}  |T-joint|={djoint:.5f}{mark}")

    print()
    print("Writing figures …")
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.plot(omega, chi2_p - chi2_p.min(), label="Pantheon+", lw=2)
    ax.plot(omega, chi2_d - chi2_d.min(), label="DES-SN5YR", lw=2)
    ax.plot(omega, chi2_b2 - chi2_b2.min(), label="DESI DR2 BAO", lw=2)
    ax.plot(omega, chi2_R - chi2_R.min(), label="CMB R", lw=2, ls="--")
    ax.plot(omega, joint - joint.min(), label="joint", color="k", lw=2.4)
    ax.axhline(1.0, color="0.6", ls=":", lw=1)
    ax.axvline(OMEGA_DE_TODAY, color="C3", ls="--", lw=1.6, label="1:11 lock")
    ax.axvline(PLANCK2018_OMEGA_DE, color="C2", ls=":", lw=1.4, label="Planck 2018")
    ax.set_xlabel(r"$\Omega_{DE}$ today")
    ax.set_ylabel(r"$\Delta\chi^2$ (scale profiled)")
    ax.set_ylim(-0.3, 10)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(RESULTS / "closeout_joint_omega.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    years = [h.year for h in HISTORY]
    vals = [h.omega_de for h in HISTORY]
    sigs = [h.sigma for h in HISTORY]
    ax.errorbar(years, vals, yerr=sigs, fmt="o", color="0.2", capsize=3)
    ax.axhline(OMEGA_DE_TODAY, color="C3", ls="--", lw=1.6, label=r"$4.9/7.1=0.69014$")
    ax.axhspan(
        PLANCK2018_OMEGA_DE - PLANCK2018_OMEGA_DE_SIGMA,
        PLANCK2018_OMEGA_DE + PLANCK2018_OMEGA_DE_SIGMA,
        color="C2",
        alpha=0.15,
        label="Planck 2018 ±1σ",
    )
    ax.set_xlabel("year")
    ax.set_ylabel(r"published $\Omega_{DE}$ (flat $\Lambda$CDM)")
    ax.set_title("previous days — no waiting for the next survey")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(RESULTS / "closeout_history.png", dpi=140)
    plt.close(fig)

    payload = {
        "lock": OMEGA_DE_TODAY,
        "omega_m": OMEGA_M,
        "probes": {
            "pantheon_plus": iv_p,
            "des_sn5yr": iv_d,
            "desi_dr1": iv_b1,
            "desi_dr2": iv_b2,
            "joint": iv_j,
            "dchi2_lock": {
                "pantheon_plus": at_lock(chi2_p),
                "des_sn5yr": at_lock(chi2_d),
                "desi_dr1": at_lock(chi2_b1),
                "desi_dr2": at_lock(chi2_b2),
                "cmb_R": at_lock(chi2_R),
                "joint": at_lock(joint),
            },
        },
        "history": hist_rows,
        "note_r_of_a": (
            "The given axioms fix today's partition T=r/(1-r). They do not "
            "map generation index to scale factor, so r(a) is not derived. "
            "The parameter-free completion is constant r, i.e. w=-1."
        ),
    }
    (RESULTS / "closeout.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {RESULTS / 'closeout.json'}")
    print(f"wrote {RESULTS / 'closeout_joint_omega.png'}")
    print(f"wrote {RESULTS / 'closeout_history.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Full viability campaign for a public write-up.

Locked: Ω_DE = 4.9/7.1 ≈ 0.69014 from the 1:11 rule.
Floated only as a *comparison* so we can quote the χ² cost of the lock.

Tests
  1. Pantheon+ SH0ES Hubble diagram, w = -1, fit H0 + 1σ profile
  2. Same, constant w free
  3. Float Ω_DE (comparison only) and Δχ² of the lock
  4. Planck 2018 Ω_Λ Gaussian pull
  5. DESI DR1 BAO shape (fit rd); inverse-ladder H0 with Planck rd
  6. CMB shift parameter R at a Planck-like H0
  7. Age of the universe at SH0ES and Planck H0
  8. Redshift-binned SN residuals
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from bao import (  # noqa: E402
    DESI_DR1,
    PLANCK_RD_MPC,
    fit_h0_given_rd,
    fit_rd,
    fit_rd_omega_de,
    n_bao,
)
from cosmology import (  # noqa: E402
    D_H,
    D_M,
    D_V,
    age_gyr,
    distance_modulus,
    shift_parameter_R,
)
from geometry import OMEGA_DE_TODAY, OMEGA_M, Q4_WEIGHT, Q5_WEIGHT, TOTAL, r  # noqa: E402
from likelihood import (  # noqa: E402
    PLANCK2018_H0,
    PLANCK2018_OMEGA_DE,
    PLANCK2018_OMEGA_DE_SIGMA,
    PLANCK2018_OMEGA_M,
    PLANCK2018_OMEGA_M_SIGMA,
    PLANCK_OMEGA_DE,
    PLANCK_R,
    PLANCK_R_SIGMA,
    SNLikelihood,
    dchi2_interval,
    load_pantheon_plus,
)

DATA_DIR = ROOT / "data"
RESULTS = ROOT / "results"
DAT_PATH = DATA_DIR / "Pantheon+SH0ES.dat"
COV_PATH = DATA_DIR / "Pantheon+SH0ES_STAT+SYS.cov"

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.bbox": "tight",
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
    }
)


def _ensure_data() -> None:
    if DAT_PATH.exists() and COV_PATH.exists():
        return
    from download_data import main as download_main

    if download_main() != 0:
        raise SystemExit("failed to download Pantheon+ data")


def _fmt_pm(interval: dict, digits: int = 3) -> str:
    return (
        f"{interval['best']:.{digits}f} "
        f"+{interval['plus']:.{digits}f}/-{interval['minus']:.{digits}f}"
    )


def _plot_hubble(data, fit_geom, fit_w, path: Path) -> None:
    z, mu = data.z, data.mu
    order = np.argsort(z)
    z_grid = np.linspace(z.min(), z.max(), 400)
    mu_geom = distance_modulus(z_grid, fit_geom["H0"], w=-1.0, omega_de=OMEGA_DE_TODAY)
    mu_w = distance_modulus(z_grid, fit_w["H0"], w=fit_w["w"], omega_de=OMEGA_DE_TODAY)
    mu_at = distance_modulus(z, fit_geom["H0"], w=-1.0, omega_de=OMEGA_DE_TODAY)

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.2), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    axes[0].scatter(z[order], mu[order], s=7, alpha=0.32, c="0.35", label="Pantheon+ $MU_{SH0ES}$")
    axes[0].plot(z_grid, mu_geom, c="C0", lw=2.0, label=rf"lock $w=-1$, $H_0={fit_geom['H0']:.2f}$")
    axes[0].plot(
        z_grid,
        mu_w,
        c="C1",
        lw=2.0,
        ls="--",
        label=rf"lock, $w={fit_w['w']:.3f}$, $H_0={fit_w['H0']:.2f}$",
    )
    axes[0].set_ylabel(r"distance modulus $\mu$")
    axes[0].legend(loc="lower right", fontsize=9)
    axes[0].set_title(r"geometric $\Omega_{DE}=4.9/7.1\approx0.69014$ on Pantheon+")
    axes[1].axhline(0.0, color="0.5", lw=0.8)
    axes[1].scatter(z, mu - mu_at, s=7, alpha=0.32, c="C0")
    axes[1].set_xlabel(r"$z_{\mathrm{HD}}$")
    axes[1].set_ylabel(r"$\mu-\mu_{w=-1}$")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_profile(x, chi2, interval, xlabel, lock_x, lock_label, path, lock_chi2=None) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(x, chi2 - chi2.min(), color="C0", lw=2)
    ax.axhline(1.0, color="0.6", ls=":", lw=1, label=r"$\Delta\chi^2=1$")
    ax.axvline(interval["best"], color="C0", ls="--", lw=1)
    if lock_x is not None:
        ax.axvline(lock_x, color="C3", ls="--", lw=1.4, label=lock_label)
        if lock_chi2 is not None:
            ax.scatter([lock_x], [lock_chi2 - chi2.min()], color="C3", zorder=5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\Delta\chi^2$")
    ax.set_ylim(-0.3, 8.0)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_omega_profile(omega, chi2, lock, free, path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    y = chi2 - chi2.min()
    ax.plot(omega, y, color="C0", lw=2, label="Pantheon+ (H0 profiled)")
    ax.axhline(1.0, color="0.6", ls=":", lw=1)
    ax.axvline(lock, color="C3", ls="--", lw=1.6, label=rf"geometry $4.9/7.1={lock:.5f}$")
    ax.axvline(PLANCK2018_OMEGA_DE, color="C2", ls="--", lw=1.4, label=rf"Planck 2018 ${PLANCK2018_OMEGA_DE:.4f}$")
    ax.axvspan(
        PLANCK2018_OMEGA_DE - PLANCK2018_OMEGA_DE_SIGMA,
        PLANCK2018_OMEGA_DE + PLANCK2018_OMEGA_DE_SIGMA,
        color="C2",
        alpha=0.12,
        label=r"Planck $\pm1\sigma$",
    )
    ax.axvline(free["omega_de"], color="C1", ls=":", lw=1.4, label=rf"SN free ${free['omega_de']:.3f}$")
    ax.set_xlabel(r"$\Omega_{DE}$ today")
    ax.set_ylabel(r"$\Delta\chi^2$ (Pantheon+)")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(-0.3, 8.0)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_binned_residuals(data, fit, path: Path) -> None:
    edges = np.array([0.01, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.4])
    mu_m = distance_modulus(data.z, fit["H0"], w=-1.0, omega_de=OMEGA_DE_TODAY)
    res = data.mu - mu_m
    # Diagonal-only visual errors; the fit itself uses the full covariance.
    diag = np.sqrt(np.clip(np.diag(data.cov), 1e-12, None))
    centers, means, errs, counts = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (data.z >= lo) & (data.z < hi)
        if m.sum() < 8:
            continue
        w = 1.0 / diag[m] ** 2
        mean = float(np.sum(w * res[m]) / np.sum(w))
        err = float(1.0 / np.sqrt(np.sum(w)))
        centers.append(0.5 * (lo + hi))
        means.append(mean)
        errs.append(err)
        counts.append(int(m.sum()))
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.errorbar(centers, means, yerr=errs, fmt="o", color="C0", capsize=3)
    for x, n in zip(centers, counts):
        ax.text(x, 0.012, str(n), ha="center", fontsize=8, color="0.4")
    ax.set_xlabel(r"$z_{\mathrm{HD}}$ bin centre")
    ax.set_ylabel(r"weighted mean $\mu-\mu_{\rm model}$ (diag only)")
    ax.set_title("binned SN residuals, locked geometry (visual; not the χ² metric)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_bao(fit, path: Path) -> None:
    H0, rd, omega_de = fit["H0"], fit["rd"], fit["omega_de"]
    z_grid = np.linspace(0.05, 2.5, 200)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3), sharex=True)
    for point in DESI_DR1:
        if point.kind == "DV":
            axes[0].errorbar(
                point.z, point.value[0], yerr=point.sigma[0], fmt="o", color="C0", capsize=3
            )
        else:
            axes[1].errorbar(
                point.z,
                point.value[0],
                yerr=point.sigma[0],
                fmt="o",
                color="C3",
                capsize=3,
            )
            axes[1].errorbar(
                point.z,
                point.value[1],
                yerr=point.sigma[1],
                fmt="s",
                color="C1",
                capsize=3,
            )
    axes[0].plot(z_grid, D_V(z_grid, H0, omega_de=omega_de) / rd, color="C3", lw=2, label="geometry")
    axes[1].plot(z_grid, D_M(z_grid, H0, omega_de=omega_de) / rd, color="C3", lw=2, label=r"$D_M/r_d$")
    axes[1].plot(z_grid, D_H(z_grid, H0, omega_de=omega_de) / rd, color="C1", lw=2, ls="--", label=r"$D_H/r_d$")
    axes[0].legend(fontsize=8)
    axes[0].set_ylabel(r"$D_V/r_d$")
    axes[1].set_ylabel(r"$D_M/r_d$, $D_H/r_d$")
    axes[0].set_xlabel(r"$z$")
    axes[1].set_xlabel(r"$z$")
    axes[0].set_title("DESI DR1 BAO, locked Ω_DE")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_contour(like, H0_grid, w_grid, path: Path, best: dict) -> None:
    chi2 = np.empty((len(w_grid), len(H0_grid)))
    for i, w in enumerate(w_grid):
        for j, H0 in enumerate(H0_grid):
            chi2[i, j] = like.chi2(H0, w=w, omega_de=OMEGA_DE_TODAY)
    d = chi2 - chi2.min()
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    # 1σ, 2σ for 2 parameters: Δχ² = 2.30, 6.17
    cs = ax.contour(H0_grid, w_grid, d, levels=[2.30, 6.17], colors=["C0", "C0"], linewidths=[2, 1])
    ax.clabel(cs, fmt={2.30: r"$1\sigma$", 6.17: r"$2\sigma$"}, fontsize=9)
    ax.plot(best["H0"], best["w"], "o", color="C1", label="best fit")
    ax.axhline(-1.0, color="C3", ls="--", lw=1.2, label=r"$w=-1$")
    ax.set_xlabel(r"$H_0$ (km s$^{-1}$ Mpc$^{-1}$)")
    ax.set_ylabel(r"$w$")
    ax.set_title(r"Pantheon+, $\Omega_{DE}$ locked")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> int:
    _ensure_data()
    RESULTS.mkdir(parents=True, exist_ok=True)

    data = load_pantheon_plus(DAT_PATH, COV_PATH)
    like = SNLikelihood(data)

    print("=" * 72)
    print("Hierarchical dark energy — full viability campaign")
    print("=" * 72)
    print(f"  1:11 geometry     r = {r:.12f}   Ω_DE = {OMEGA_DE_TODAY:.12f}")
    print(f"  Ω_M               = {OMEGA_M:.12f}")
    print(f"  Pantheon+         N = {data.n} / {data.n_raw}   z = {data.z.min():.4f}–{data.z.max():.4f}")
    print(f"  DESI DR1 BAO      N = {n_bao()} distance ratios")
    print()

    print("SN 1: lock Ω_DE, w=-1, fit H0 …")
    fit1 = like.fit_h0(w=-1.0, omega_de=OMEGA_DE_TODAY)
    H0_grid = np.linspace(fit1["H0"] - 3.0, fit1["H0"] + 3.0, 41)
    chi2_H0 = like.profile(H0_grid, "H0", w=-1.0, omega_de=OMEGA_DE_TODAY)
    H0_int = dchi2_interval(H0_grid, chi2_H0)
    fit1["H0_err_plus"] = H0_int["plus"]
    fit1["H0_err_minus"] = H0_int["minus"]
    print(f"  H0  = {_fmt_pm(H0_int, 3)} km/s/Mpc")
    print(f"  χ²  = {fit1['chi2']:.3f}  / dof {fit1['dof']}  = {fit1['chi2_dof']:.4f}")

    print("SN 2: lock Ω_DE, fit H0 and w …")
    fit2 = like.fit_h0_w(omega_de=OMEGA_DE_TODAY, H0_guess=fit1["H0"])
    w_grid = np.linspace(-1.4, -0.6, 33)
    chi2_w = like.profile(w_grid, "w", omega_de=OMEGA_DE_TODAY)
    w_int = dchi2_interval(w_grid, chi2_w)
    fit2["w_err_plus"] = w_int["plus"]
    fit2["w_err_minus"] = w_int["minus"]
    print(f"  w   = {_fmt_pm(w_int, 3)}")
    print(f"  H0  = {fit2['H0']:.3f}")
    print(f"  χ²  = {fit2['chi2']:.3f}   Δχ² vs w=-1 = {fit2['chi2']-fit1['chi2']:+.3f}")

    print("SN comparison: float Ω_DE (NOT a retune) …")
    fit_free = like.fit_h0_omega_de(w=-1.0, H0_guess=fit1["H0"], omega_guess=OMEGA_DE_TODAY)
    omega_grid = np.linspace(0.55, 0.80, 41)
    chi2_om = like.profile(omega_grid, "omega_de", w=-1.0)
    om_int = dchi2_interval(omega_grid, chi2_om)
    dchi2_lock = fit1["chi2"] - fit_free["chi2"]
    print(f"  free Ω_DE = {_fmt_pm(om_int, 4)}")
    print(f"  lock cost Δχ² = {dchi2_lock:.3f}   ({np.sqrt(max(dchi2_lock, 0)):.2f}σ equivalent, 1 dof)")

    print("Planck 2018 Gaussian pull on the derived density …")
    pull_planck = (OMEGA_DE_TODAY - PLANCK2018_OMEGA_DE) / PLANCK2018_OMEGA_DE_SIGMA
    pull_m = (OMEGA_M - PLANCK2018_OMEGA_M) / PLANCK2018_OMEGA_M_SIGMA
    print(f"  Ω_DE geom {OMEGA_DE_TODAY:.5f} vs Planck {PLANCK2018_OMEGA_DE:.4f} ± {PLANCK2018_OMEGA_DE_SIGMA:.4f}")
    print(f"  pull = {pull_planck:+.2f}σ  (Ω_M pull {pull_m:+.2f}σ)")

    print("DESI DR1 BAO shape, dummy H0=70, fit rd …")
    bao_geom = fit_rd(70.0, omega_de=OMEGA_DE_TODAY, w=-1.0)
    bao_planck = fit_rd(70.0, omega_de=PLANCK2018_OMEGA_DE, w=-1.0)
    bao_free = fit_rd_omega_de(H0=70.0, w=-1.0)
    print(f"  geom  χ² = {bao_geom['chi2']:.2f} / {bao_geom['dof']}   rd(h=0.7) = {bao_geom['rd']:.2f} Mpc")
    print(f"  Planck χ² = {bao_planck['chi2']:.2f}")
    print(f"  free Ω_DE = {bao_free['omega_de']:.4f}  χ² = {bao_free['chi2']:.2f}  Δχ² lock = {bao_geom['chi2']-bao_free['chi2']:+.2f}")

    print("Inverse distance ladder: Planck rd, locked Ω_DE, fit H0 from BAO …")
    ladder = fit_h0_given_rd(rd=PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY)
    print(f"  H0(BAO+rd) = {ladder['H0']:.2f} km/s/Mpc   χ² = {ladder['chi2']:.2f}")
    print("  (this is the CMB/BBN ruler; it is not the SH0ES intercept)")

    print("CMB shift parameter R at Planck H0 …")
    R_geom = shift_parameter_R(PLANCK2018_H0, omega_de=OMEGA_DE_TODAY, w=-1.0)
    R_planck_om = shift_parameter_R(PLANCK2018_H0, omega_de=PLANCK2018_OMEGA_DE, w=-1.0)
    R_pull = (R_geom - PLANCK_R) / PLANCK_R_SIGMA
    print(f"  R_geom = {R_geom:.5f}   R_Planck prior = {PLANCK_R:.4f} ± {PLANCK_R_SIGMA:.4f}")
    print(f"  pull = {R_pull:+.2f}σ   (R at Planck Ω_DE = {R_planck_om:.5f})")

    print("Age of the universe …")
    age_shoes = age_gyr(fit1["H0"], w=-1.0, omega_de=OMEGA_DE_TODAY)
    age_planck_h0 = age_gyr(PLANCK2018_H0, w=-1.0, omega_de=OMEGA_DE_TODAY)
    print(f"  t0(H0=SN {fit1['H0']:.2f})   = {age_shoes:.3f} Gyr")
    print(f"  t0(H0=Planck {PLANCK2018_H0:.2f}) = {age_planck_h0:.3f} Gyr")

    print()
    print("Writing figures …")
    _plot_hubble(data, fit1, fit2, RESULTS / "hubble_diagram.png")
    _plot_profile(
        H0_grid,
        chi2_H0,
        H0_int,
        r"$H_0$ (km s$^{-1}$ Mpc$^{-1}$)",
        None,
        "",
        RESULTS / "profile_H0.png",
    )
    _plot_profile(
        w_grid,
        chi2_w,
        w_int,
        r"$w$",
        -1.0,
        r"$w=-1$",
        RESULTS / "profile_w.png",
        lock_chi2=fit1["chi2"],
    )
    _plot_omega_profile(omega_grid, chi2_om, OMEGA_DE_TODAY, fit_free, RESULTS / "profile_omega_de.png")
    _plot_binned_residuals(data, fit1, RESULTS / "residuals_binned.png")
    _plot_bao(bao_geom, RESULTS / "bao_desi.png")
    H0c = np.linspace(fit2["H0"] - 1.5, fit2["H0"] + 1.5, 25)
    wc = np.linspace(-1.25, -0.70, 25)
    _plot_contour(like, H0c, wc, RESULTS / "contour_H0_w.png", fit2)

    payload = {
        "geometry": {
            "Q4_WEIGHT": Q4_WEIGHT,
            "Q5_WEIGHT": Q5_WEIGHT,
            "TOTAL": TOTAL,
            "r": r,
            "OMEGA_DE_TODAY": OMEGA_DE_TODAY,
            "OMEGA_M": OMEGA_M,
            "locked": True,
        },
        "sn": {
            "n": data.n,
            "locked_w_minus_1": fit1,
            "locked_free_w": fit2,
            "free_omega_de_comparison_only": fit_free,
            "H0_1sigma": H0_int,
            "w_1sigma": w_int,
            "omega_de_1sigma_if_floated": om_int,
            "delta_chi2_lock_vs_free_omega": dchi2_lock,
            "delta_chi2_free_w": fit2["chi2"] - fit1["chi2"],
        },
        "planck2018": {
            "omega_de": PLANCK2018_OMEGA_DE,
            "omega_de_sigma": PLANCK2018_OMEGA_DE_SIGMA,
            "omega_m": PLANCK2018_OMEGA_M,
            "pull_sigma": pull_planck,
            "design_doc_comparison_value": PLANCK_OMEGA_DE,
        },
        "bao_desi_dr1": {
            "locked": bao_geom,
            "planck_omega": bao_planck,
            "free_omega_comparison_only": bao_free,
            "delta_chi2_lock_vs_free": bao_geom["chi2"] - bao_free["chi2"],
            "inverse_ladder": ladder,
        },
        "cmb_shift_R": {
            "R_geom": R_geom,
            "R_planck_prior": PLANCK_R,
            "R_sigma": PLANCK_R_SIGMA,
            "pull_sigma": R_pull,
            "H0_used": PLANCK2018_H0,
        },
        "age_gyr": {
            "at_sn_H0": age_shoes,
            "at_planck_H0": age_planck_h0,
        },
    }
    (RESULTS / "validation.json").write_text(json.dumps(payload, indent=2) + "\n")

    verdict_lines = [
        "Viability (present-day density only, not a dynamical DE theory)",
        f"  Planck 2018 Ω_Λ pull: {pull_planck:+.2f}σ",
        f"  SN lock cost vs free Ω_DE: Δχ² = {dchi2_lock:.2f}",
        f"  SN w vs -1: Δχ² = {fit2['chi2']-fit1['chi2']:.2f}",
        f"  DESI BAO lock cost vs free Ω_DE: Δχ² = {bao_geom['chi2']-bao_free['chi2']:.2f}",
        f"  CMB R pull at Planck H0: {R_pull:+.2f}σ",
        "  The 1:11 ratio was not varied.",
    ]
    report = "\n".join(
        [
            "Hierarchical DE — full validation",
            f"Ω_DE = {OMEGA_DE_TODAY:.12f}  (locked, 1 part q=4 + 11 parts q=5)",
            "",
            f"SN H0 (w=-1)     {_fmt_pm(H0_int, 3)} km/s/Mpc   χ²/dof = {fit1['chi2_dof']:.4f}",
            f"SN w             {_fmt_pm(w_int, 3)}   Δχ²(w=-1) = {fit2['chi2']-fit1['chi2']:+.3f}",
            f"SN free Ω_DE     {_fmt_pm(om_int, 4)}   lock Δχ² = {dchi2_lock:+.3f}",
            f"Planck 2018 pull {pull_planck:+.2f}σ",
            f"DESI BAO geom χ² {bao_geom['chi2']:.2f}/{bao_geom['dof']}  free Ω_DE = {bao_free['omega_de']:.4f}",
            f"BAO+Planck rd H0 {ladder['H0']:.2f} km/s/Mpc",
            f"CMB R            {R_geom:.5f}  ({R_pull:+.2f}σ vs {PLANCK_R:.4f}±{PLANCK_R_SIGMA:.4f})",
            f"t0               {age_shoes:.2f} Gyr (SN H0) / {age_planck_h0:.2f} Gyr (Planck H0)",
            "",
            *verdict_lines,
            "",
        ]
    )
    (RESULTS / "validation.txt").write_text(report)
    print()
    print(report)
    print(f"wrote {RESULTS / 'validation.json'}")
    print(f"wrote {RESULTS / 'validation.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

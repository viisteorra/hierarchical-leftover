#!/usr/bin/env python3
"""SN Hubble-diagram fits with Ω_DE locked to 49/71.

1. w = -1, fit H0.
2. constant w free, still lock Ω_DE.
3. Planck Ω_DE ≈ 0.689 reference, w = -1, fit H0.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cosmology import distance_modulus  # noqa: E402
from geometry import (  # noqa: E402
    OMEGA_DE_TODAY,
    OMEGA_M,
    Q4_WEIGHT,
    Q5_WEIGHT,
    TOTAL,
    r,
)
from likelihood import (  # noqa: E402
    PLANCK_OMEGA_DE,
    SNLikelihood,
    load_pantheon_plus,
)

DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
DAT_PATH = DATA_DIR / "Pantheon+SH0ES.dat"
COV_PATH = DATA_DIR / "Pantheon+SH0ES_STAT+SYS.cov"


def _ensure_data() -> None:
    if DAT_PATH.exists() and COV_PATH.exists():
        return
    sys.path.insert(0, str(ROOT / "scripts"))
    from download_data import main as download_main

    if download_main() != 0:
        raise SystemExit("failed to download Pantheon+ data")


def _fmt_fit(label: str, fit: dict, chi2_ref: float | None = None) -> str:
    lines = [
        f"{label}",
        f"  Ω_DE        = {fit['omega_de']:.8f}",
        f"  w           = {fit['w']:.6f}",
        f"  H0          = {fit['H0']:.4f} km/s/Mpc",
        f"  χ²          = {fit['chi2']:.3f}",
        f"  N, npar, dof = {fit['n']}, {fit['npar']}, {fit['dof']}",
        f"  χ²/dof      = {fit['chi2_dof']:.4f}",
        f"  success     = {fit['success']}",
    ]
    if chi2_ref is not None:
        lines.append(f"  Δχ²         = {fit['chi2'] - chi2_ref:+.3f}")
    return "\n".join(lines)


def _plot_hubble(data, fit_geom, fit_w, out_path: Path) -> None:
    z = data.z
    mu = data.mu
    order = np.argsort(z)
    z_s = z[order]
    mu_s = mu[order]
    z_grid = np.linspace(z.min(), z.max(), 400)

    mu_geom = distance_modulus(z_grid, fit_geom["H0"], w=-1.0, omega_de=OMEGA_DE_TODAY)
    mu_w = distance_modulus(z_grid, fit_w["H0"], w=fit_w["w"], omega_de=OMEGA_DE_TODAY)
    mu_geom_at_z = distance_modulus(z, fit_geom["H0"], w=-1.0, omega_de=OMEGA_DE_TODAY)

    fig, axes = plt.subplots(
        2, 1, figsize=(8.5, 7.0), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    axes[0].scatter(z_s, mu_s, s=8, alpha=0.35, color="0.3", label="Pantheon+ MU_SH0ES")
    axes[0].plot(
        z_grid,
        mu_geom,
        color="C0",
        lw=2,
        label=rf"geom $w=-1$, $H_0={fit_geom['H0']:.2f}$",
    )
    axes[0].plot(
        z_grid,
        mu_w,
        color="C1",
        lw=2,
        ls="--",
        label=rf"geom $w={fit_w['w']:.3f}$, $H_0={fit_w['H0']:.2f}$",
    )
    axes[0].set_ylabel(r"distance modulus $\mu$")
    axes[0].legend(loc="lower right", fontsize=9)
    axes[0].set_title(r"locked $\Omega_{DE}=0.69014$ on Pantheon+")

    axes[1].axhline(0.0, color="0.5", lw=0.8)
    axes[1].scatter(z, mu - mu_geom_at_z, s=8, alpha=0.35, color="C0")
    axes[1].set_xlabel(r"$z_{\mathrm{HD}}$")
    axes[1].set_ylabel(r"$\mu - \mu_{w=-1}$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def main() -> int:
    sys.path.insert(0, str(ROOT / "scripts"))
    _ensure_data()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    data = load_pantheon_plus(DAT_PATH, COV_PATH)
    like = SNLikelihood(data)

    print("=" * 64)
    print("Hierarchical dark energy — locked geometric core")
    print("=" * 64)
    print(f"  q4:q5 weights     = {Q4_WEIGHT}:{Q5_WEIGHT}  (total {TOTAL})")
    print(f"  r                 = {r:.12f}   (= 4.9/12)")
    print(f"  Ω_DE today        = {OMEGA_DE_TODAY:.12f}   (= 4.9/7.1 ≈ 0.69014)")
    print(f"  Ω_M               = {OMEGA_M:.12f}")
    print(f"  sample            = Pantheon+ SH0ES, zHD > {data.zmin}")
    print(f"  N (after z cut)   = {data.n} / {data.n_raw}")
    print(f"  z range           = {data.z.min():.4f} … {data.z.max():.4f}")
    print()

    print("Running experiment 1: w = -1, Ω_DE locked, fit H0 …")
    fit1 = like.fit_h0(w=-1.0, omega_de=OMEGA_DE_TODAY)
    print(_fmt_fit("Experiment 1  (geometric Λ, fit H0)", fit1))
    print()

    print("Running experiment 2: Ω_DE locked, fit H0 and constant w …")
    fit2 = like.fit_h0_w(
        omega_de=OMEGA_DE_TODAY,
        H0_guess=fit1["H0"],
        w_guess=-1.0,
    )
    print(_fmt_fit("Experiment 2  (geometric Ω_DE, free w)", fit2, chi2_ref=fit1["chi2"]))
    print()

    print("Running Planck reference: Ω_DE = 0.689, w = -1, fit H0 …")
    fit_planck = like.fit_h0(w=-1.0, omega_de=PLANCK_OMEGA_DE)
    print(
        _fmt_fit(
            "Reference     (Planck Ω_DE ≈ 0.689, w = -1)",
            fit_planck,
            chi2_ref=fit1["chi2"],
        )
    )
    print()

    payload = {
        "geometry": {
            "Q4_WEIGHT": Q4_WEIGHT,
            "Q5_WEIGHT": Q5_WEIGHT,
            "TOTAL": TOTAL,
            "r": r,
            "OMEGA_DE_TODAY": OMEGA_DE_TODAY,
            "OMEGA_M": OMEGA_M,
            "note": "1:11 ratio and Ω_DE are geometrically derived and were not fitted",
        },
        "data": {
            "catalogue": "Pantheon+SH0ES",
            "observable": "MU_SH0ES",
            "covariance": "STAT+SYS",
            "zmin": data.zmin,
            "n": data.n,
            "n_raw": data.n_raw,
            "z_min_used": float(data.z.min()),
            "z_max_used": float(data.z.max()),
        },
        "experiment_1_w_minus_1": fit1,
        "experiment_2_free_w": fit2,
        "reference_planck_omega_de": fit_planck,
        "delta_chi2_free_w_minus_locked_w": fit2["chi2"] - fit1["chi2"],
        "delta_chi2_planck_minus_geometric": fit_planck["chi2"] - fit1["chi2"],
    }

    json_path = RESULTS_DIR / "fit_results.json"
    txt_path = RESULTS_DIR / "fit_results.txt"
    json_path.write_text(json.dumps(payload, indent=2) + "\n")

    report = "\n".join(
        [
            "Hierarchical DE — first campaign",
            "",
            f"r = {r:.12f}",
            f"Ω_DE = {OMEGA_DE_TODAY:.12f}  (locked, 1:11 geometry)",
            f"Ω_M  = {OMEGA_M:.12f}",
            "",
            _fmt_fit("Experiment 1  (geometric Λ, fit H0)", fit1),
            "",
            _fmt_fit("Experiment 2  (geometric Ω_DE, free w)", fit2, chi2_ref=fit1["chi2"]),
            "",
            _fmt_fit(
                "Reference     (Planck Ω_DE ≈ 0.689, w = -1)",
                fit_planck,
                chi2_ref=fit1["chi2"],
            ),
            "",
            "The 1:11 ratio was not varied. Δχ² vs Planck tests whether the",
            "derived density is viable; it is not a license to retune T.",
            "",
        ]
    )
    txt_path.write_text(report)

    plot_path = RESULTS_DIR / "hubble_diagram.png"
    print("Writing Hubble diagram …")
    _plot_hubble(data, fit1, fit2, plot_path)

    print(f"wrote {json_path}")
    print(f"wrote {txt_path}")
    print(f"wrote {plot_path}")
    print()
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())

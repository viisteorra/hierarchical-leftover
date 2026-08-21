#!/usr/bin/env python3
"""Nails: first-principles audit + Hubble claim + swapped-measure kills.

Does not retune 1:11. Does not put H0 into f.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bao import DESI_DR2, PLANCK_RD_MPC, bao_chi2, fit_h0_given_rd  # noqa: E402
from geometry import OMEGA_DE_TODAY, W  # noqa: E402
from hubble import SIG_BAO, SIG_SN, solution, tension  # noqa: E402
from likelihood import SNLikelihood, load_pantheon_plus  # noqa: E402
from ruler import leftover_octave  # noqa: E402
from spacetime import F_BAO, F_CMB, ZSTAR_RECOMB, f_for, measure_for  # noqa: E402
from uniqueness import derivation_chain, period_from_fifth_convergent  # noqa: E402

RESULTS = ROOT / "results"
DATA = ROOT / "data"
nails: list[dict] = []


def nail(kid: str, fired: bool, detail: str) -> None:
    nails.append({"id": kid, "fired": bool(fired), "detail": detail})
    print(f"  [{'KILL' if fired else 'live'}] {kid}: {detail}")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    print("=" * 72)
    print("NAILS  first principles + Hubble claim")
    print("=" * 72)

    print("\n--- N0  derivation chain (only P is unproved) ---")
    chain = derivation_chain()
    free = [row for row in chain if row["free"]]
    for row in chain:
        tag = "PRIMITIVE" if row["free"] else "lemma"
        print(f"  [{tag}] {row['symbol']}: {row['value']}  ← {row['origin']}")
    nail("N0_only_P", len(free) != 1 or free[0]["symbol"] != "P", f"{len(free)} free")
    n, k, _ = period_from_fifth_convergent()
    nail("N0_u1_cf", (n, k) != (12, 7), f"period {n}/{k}")
    nail("N0_w", W != -1.0, f"w={W}")
    nail("N0_T", abs(OMEGA_DE_TODAY - 49 / 71) >= 1e-12, "T")

    print("\n--- N1  leftover source has no H0, no χ² ---")
    geo = (ROOT / "src" / "geometry.py").read_text().lower()
    rul = (ROOT / "src" / "ruler.py").read_text().lower()
    nail("N1_geo_fit", "chi2" in geo or "minimize" in geo or "73.47" in geo, "geometry fit")
    nail("N1_rul_h0", "h0_sn" in rul or "73.47" in rul, "ruler H0")

    print("\n--- N2  Lemma C: measure from observable type ---")
    nail("N2_bao", measure_for("BAO") != "planar" or abs(f_for("BAO") - F_BAO) > 1e-15, "BAO")
    nail("N2_cmb", measure_for("CMB") != "infinite" or abs(f_for("CMB") - F_CMB) > 1e-15, "CMB")

    like = SNLikelihood(
        load_pantheon_plus(DATA / "Pantheon+SH0ES.dat", DATA / "Pantheon+SH0ES_STAT+SYS.cov")
    )
    sn = like.fit_h0(omega_de=OMEGA_DE_TODAY)
    sn_free = like.fit_h0_omega_de(H0_guess=sn["H0"], omega_guess=OMEGA_DE_TODAY)
    bao = fit_h0_given_rd(PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    h0_sn, h0_bao, chi0 = float(sn["H0"]), float(bao["H0"]), float(bao["chi2"])
    s = solution(h0_sn=h0_sn, h0_bao=h0_bao)

    print("\n--- N3  49/71 is not the SN χ² minimum ---")
    print(f"  lock Ω={OMEGA_DE_TODAY:.5f}  SN-free Ω={sn_free['omega_de']:.5f}  "
          f"Δχ²={sn['chi2']-sn_free['chi2']:.2f}")
    nail(
        "N3_not_sn_fit",
        abs(sn_free["omega_de"] - 49 / 71) < 0.002,
        "SN χ² min landed on 49/71 — would look like a fit",
    )

    print("\n--- N4  Hubble scale ---")
    print(f"  raw BAO {s['raw_bao_sigma']:.2f}σ → map {s['map_bao_sigma']:.3f}σ")
    print(f"  raw θ*  {s['raw_th_sigma']:.2f}σ → map {s['map_th_sigma']:.3f}σ")
    print(f"  forward SN/f_□ {s['pred_bao_from_sn']:.3f} vs DESI {h0_bao:.3f}  {s['fwd_bao_sigma']:.3f}σ")
    print(f"  corollary f_∞/f_□ {s['ratio_f']:.6f} vs DESI/θ* {s['ratio_h']:.6f}  {s['early_early_sigma']:.3f}σ")
    nail("N4_bao", s["map_bao_sigma"] >= 1.0, f"{s['map_bao_sigma']:.3f}σ")
    nail("N4_cmb_map", s["map_th_sigma"] >= 1.0, f"{s['map_th_sigma']:.3f}σ")
    nail("N4_raw_was_tension", s["raw_bao_sigma"] < 4.0, f"raw {s['raw_bao_sigma']:.2f}σ")
    nail("N4_closed_vs_raw", s["map_bao_sigma"] > s["raw_bao_sigma"] / 10.0, "not 10× closed")

    print("\n--- N5  swapped measures must fail ---")
    swap_bao = tension(h0_sn, F_CMB * h0_bao, SIG_SN, SIG_BAO)
    swap_cmb = tension(h0_sn, F_BAO * s["h0_th"], SIG_SN, 0.54)
    dchi_lic = bao_chi2(h0_sn, PLANCK_RD_MPC / F_BAO, sample=DESI_DR2) - chi0
    dchi_swp = bao_chi2(h0_sn, PLANCK_RD_MPC / F_CMB, sample=DESI_DR2) - chi0
    print(f"  licensed BAO {s['map_bao_sigma']:.3f}σ  Δχ²={dchi_lic:+.2f}")
    print(f"  f_∞ on BAO  {swap_bao:.3f}σ  Δχ²={dchi_swp:+.2f}")
    print(f"  licensed θ* {s['map_th_sigma']:.3f}σ")
    print(f"  f_□ on θ*   {swap_cmb:.3f}σ")
    nail("N5_swap_bao_k2", dchi_swp < 10.0, f"inf-on-BAO Δχ²={dchi_swp:.2f} (must be ugly)")
    nail("N5_swap_bao_worse", swap_bao <= s["map_bao_sigma"], "inf-on-BAO not worse")
    nail("N5_swap_cmb_worse", swap_cmb <= s["map_th_sigma"], "planar-on-CMB not worse")
    nail("N5_k2_licensed", dchi_lic >= 10.0, f"licensed Δχ²={dchi_lic:.2f}")

    print("\n--- N6  z_d foil (photon z* is B; drag is not leftover's scale) ---")
    z_d = 1059.87
    f_d = leftover_octave(z_d)["f_planar"]
    sig_d = tension(h0_sn, f_d * h0_bao, SIG_SN, SIG_BAO)
    print(f"  z* licensed {ZSTAR_RECOMB}  {s['map_bao_sigma']:.3f}σ")
    print(f"  z_d foil    {z_d}  f_□={f_d:.6f}  {sig_d:.2f}σ")
    nail("N6_zd_beats", sig_d < s["map_bao_sigma"], "drag leftover beat photon leftover")

    print("\n--- N7  forward CAMB (from hubble.json if present) ---")
    hubp = RESULTS / "hubble.json"
    camb = {}
    if hubp.exists():
        hub = json.loads(hubp.read_text())
        camb = hub.get("camb_predicted") or {}
        if camb:
            tp = float(camb["theta_pull"])
            print(f"  SN/f_∞ H0={camb['H0']:.4f}  100θ*={camb['100theta']:.5f}  {tp:+.2f}σ  ℓ={camb['peak_l']}")
            nail("N7_theta", abs(tp) >= 3.0, f"forward θ* {tp:+.2f}σ")
            nail("N7_peak", camb.get("peak_l") not in range(215, 226), f"ℓ={camb.get('peak_l')}")
        else:
            nail("N7_theta", False, "no CAMB block (run scripts/run_hubble.py)")
    else:
        nail("N7_theta", False, "hubble.json missing (run scripts/run_hubble.py)")

    n_kill = sum(1 for k in nails if k["fired"])
    print("\n" + "=" * 72)
    print(f"NAILS  {n_kill}/{len(nails)} fired")
    if n_kill == 0:
        print("CLAIM (given P): Hubble scale split is leftover of one FRW.")
        print("  0 free cosmological parameters in f and Ω_DE.")
        print("  Two measures are Lemma C (observable type), not a fit.")
        print("  P itself is not proved. Survive ≠ P true.")
    else:
        print("Fired:")
        for k in nails:
            if k["fired"]:
                print(f"  - {k['id']}: {k['detail']}")
    print("=" * 72)
    out = {
        "n_fired": n_kill,
        "n_nails": len(nails),
        "nails": nails,
        "solution": s,
        "sn_free_omega": sn_free["omega_de"],
        "swap": {"bao_inf": swap_bao, "cmb_box": swap_cmb, "dchi_swap": dchi_swp, "dchi_lic": dchi_lic},
        "z_d": {"sigma": sig_d, "f": f_d},
        "camb": camb,
        "chain": chain,
    }
    path = RESULTS / "nails.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"Wrote {path}")
    return 0 if n_kill == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

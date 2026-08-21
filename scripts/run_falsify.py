#!/usr/bin/env python3
"""Exhaustive attempt to disprove the closed spacetime.

Kill criteria (any one is enough):
  K1  leftover f leaves SN vs BAO+rd split ≥ 3σ (SH0ES cal budget)
  K2  BAO χ² at the (H0_L, rd_L) pair worse than lock by Δχ² ≥ 10
  K3  CMB R in the Einstein frame ≥ 3σ from Planck
  K4  early ω_m ≥ 3σ from Planck 0.1430±0.0011
  K5  baryon t_U ≥ 3σ from Valcin 13.57±0.27
  K6  t_U < t_GC 13.39 (clusters older than the universe) at > 2σ
  K7  SN lock Δχ² vs free Ω_DE ≥ 9 (≈3σ)
  K8  time-delay H0_L vs delay-only past samples ≥ 3σ
      (H0_TD = H0_L in one continuum; 67 from non-delay lenses is not K8)
  K9  mixed frames (SN H0 + Planck rd, no f) looking GOOD would kill
      the claim that f is doing work; it must look BAD
  K10 f from HS z* vs bias-corr z* disagree enough to push K1 over 3σ
      at BOTH z*  (z* fragility)

Survive = none of K1–K10 fire.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bao import DESI_DR2, PLANCK_RD_MPC, bao_chi2, fit_h0_given_rd, fit_rd  # noqa: E402
from cosmology import (  # noqa: E402
    C_KMS,
    OMEGA_R_H2,
    age_gyr,
    comoving_chi,
    shift_parameter_R,
    z_star_hu_sugiyama,
)
from geometry import OMEGA_DE_TODAY, OMEGA_M  # noqa: E402
from likelihood import (  # noqa: E402
    PLANCK2018_H0,
    PLANCK2018_OMEGA_DE,
    PLANCK_R,
    PLANCK_R_SIGMA,
    SNLikelihood,
    load_des_dovekie,
    load_pantheon_plus,
)
from spacetime import (  # noqa: E402
    F_ANSATZ_T10,
    F_AXIOM,
    F_BAO,
    F_CAMB,
    F_CLOCK,
    F_CMB,
    F_FINITE4,
    F_HS,
    F_MIXED_2D,
    ZSTAR_HS_LOCK,
)
from ruler import leftover_f, mix_stretch_f, leftover_octave  # noqa: E402

RESULTS = ROOT / "results"
DATA = ROOT / "data"

H0_SN = 73.47106273986928
SIG_SN = 1.04
SIG_BAO = 0.217
PLANCK_OMMH2 = 0.1430
PLANCK_OMMH2_SIGMA = 0.0011
PLANCK_RD_SIGMA = 0.26
VALCIN_TU = 13.57
VALCIN_ERR = math.sqrt(0.15**2 + 0.23**2)
VALCIN_TGC = 13.39
VALCIN_TGC_ERR = math.sqrt(0.10**2 + 0.23**2)
TDCOSMO_H0 = 72.1
TDCOSMO_SIG = 3.85  # (4.0+3.7)/2
H0LICOW_H0 = 73.3
H0LICOW_SIG = 1.75
OMEGA_B_H2 = 0.02236


def tens(a, b, sa, sb):
    return abs(a - b) / math.sqrt(sa**2 + sb**2)


def kill(kills, kid, fired, detail):
    kills.append({"id": kid, "fired": bool(fired), "detail": detail})
    flag = "KILL" if fired else "live"
    print(f"  [{flag}] {kid}: {detail}")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    print("=" * 72)
    print("FALSIFY  closed spacetime  (try to kill it)")
    print("=" * 72)
    print(f"  Ω_DE={OMEGA_DE_TODAY:.6f}  H0_SN={H0_SN:.3f}")
    print(f"  A2 BAO  f_□=2^{{φ·(12/11)/(1-r/11)}}={F_BAO:.6f}")
    print(f"  A2 CMB  f_∞=2^{{φ·⟨e⟩}}={F_CMB:.6f}  (CAMB {F_CAMB:.6f})")
    print(f"  clock={F_CLOCK:.6f}  4D cutoff={F_FINITE4:.6f}  ansatz={F_ANSATZ_T10:.6f}")

    like = SNLikelihood(
        load_pantheon_plus(DATA / "Pantheon+SH0ES.dat", DATA / "Pantheon+SH0ES_STAT+SYS.cov")
    )
    sn_lock = like.fit_h0(omega_de=OMEGA_DE_TODAY)
    sn_free = like.fit_h0_omega_de(H0_guess=H0_SN, omega_guess=OMEGA_DE_TODAY)
    bao_early = fit_h0_given_rd(PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    H0_E_data = bao_early["H0"]  # inverse ladder, no f
    chi_bao_early = bao_early["chi2"]

    kills = []
    scoreboard = {}

    # --- K1, K10: BAO uses f_□, CMB θ* uses f_∞. Do not mix. ---
    print("\n--- K1 / K10  intercept split (measure matches observable) ---")
    H0_TH = 67.5991
    splits = {}
    F_BAO_HS = leftover_octave(ZSTAR_HS_LOCK)["f_planar"]
    F_CMB_HS = leftover_f(ZSTAR_HS_LOCK)
    for name, f, H0E, sE in (
        ("bao_planar", F_BAO, H0_E_data, SIG_BAO),
        ("cmb_inf", F_CMB, H0_TH, 0.54),
        ("bao_planar_HS", F_BAO_HS, H0_E_data, SIG_BAO),
        ("cmb_inf_HS", F_CMB_HS, H0_TH, 0.54),
        ("clock_on_bao", F_CLOCK, H0_E_data, SIG_BAO),
        ("ansatz_on_bao", F_ANSATZ_T10, H0_E_data, SIG_BAO),
        ("finite4_on_bao", F_FINITE4, H0_E_data, SIG_BAO),
    ):
        H0_L_pred = f * H0E
        sp = H0_SN - H0_L_pred
        sig = tens(H0_SN, H0_L_pred, SIG_SN, sE)
        chi = bao_chi2(H0_SN, PLANCK_RD_MPC / f, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
        splits[name] = {
            "f": f, "H0_E": H0E, "H0_L_pred": H0_L_pred, "split": sp, "sig": sig,
            "bao_chi2": chi, "rd_L": PLANCK_RD_MPC / f,
        }
        print(f"  {name:16s} f={f:.6f}  {H0E:.2f}→{H0_L_pred:.3f}  split={sp:+.3f}  {sig:.2f}σ  "
              f"BAO χ²(H0_SN,rd/f)={chi:.2f}")

    k1 = splits["bao_planar"]["sig"] >= 3.0 or splits["cmb_inf"]["sig"] >= 3.0
    kill(
        kills,
        "K1",
        k1,
        f"BAO f_□ {splits['bao_planar']['sig']:.2f}σ  CMB f_∞ {splits['cmb_inf']['sig']:.2f}σ (kill if either ≥3)",
    )
    k10 = (splits["bao_planar"]["sig"] >= 3.0 and splits["bao_planar_HS"]["sig"] >= 3.0) or (
        splits["cmb_inf"]["sig"] >= 3.0 and splits["cmb_inf_HS"]["sig"] >= 3.0
    )
    kill(
        kills,
        "K10",
        k10,
        f"BAO {splits['bao_planar']['sig']:.2f}/{splits['bao_planar_HS']['sig']:.2f}σ  "
        f"CMB {splits['cmb_inf']['sig']:.2f}/{splits['cmb_inf_HS']['sig']:.2f}σ both-z* ≥3",
    )

    # --- K2 BAO χ² at the BAO measure, not f_∞ ---
    print("\n--- K2  BAO χ² (planar f_□ only) ---")
    dchi = splits["bao_planar"]["bao_chi2"] - chi_bao_early
    chi_hyb = bao_chi2(H0_SN, PLANCK_RD_MPC / F_CMB, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    print(f"  early (H0_E, rd_P) χ²={chi_bao_early:.2f}")
    print(f"  local (H0_SN, rd/f_□) χ²={splits['bao_planar']['bao_chi2']:.2f}  Δ={dchi:+.2f}")
    print(f"  hybrid (H0_SN, rd/f_∞) χ²={chi_hyb:.2f}  (must stay unused)")
    kill(kills, "K2", dchi >= 10.0, f"Δχ²={dchi:.2f} (kill if ≥10)")

    # --- K9 mixed frames must be BAD ---
    print("\n--- K9  mixed frames (no f) must be bad ---")
    chi_mixed = bao_chi2(H0_SN, PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    print(f"  SN H0 + Planck rd χ²={chi_mixed:.1f}  (tension as χ²)")
    kill(kills, "K9", chi_mixed < 30.0,
         f"mixed χ²={chi_mixed:.1f} (kill if this looked fine, <30)")

    # --- K3 R ---
    print("\n--- K3  CMB R in Einstein/early frame ---")
    # early frame H0 for R: H0_E = H0_SN/f  (weak H0 dep)
    H0_E_hs = H0_SN / F_CMB
    R_e = shift_parameter_R(H0_E_hs, omega_de=OMEGA_DE_TODAY)
    R_p = shift_parameter_R(PLANCK2018_H0, omega_de=OMEGA_DE_TODAY)
    pull_e = (R_e - PLANCK_R) / PLANCK_R_SIGMA
    print(f"  R(H0_E={H0_E_hs:.2f})={R_e:.4f}  pull={pull_e:+.2f}σ")
    print(f"  R(Planck H0, lock Ω)={R_p:.4f}")
    kill(kills, "K3", abs(pull_e) >= 3.0, f"R pull {pull_e:+.2f}σ (kill if |σ|≥3)")

    # --- K4 ω_m ---
    print("\n--- K4  early ω_m ---")
    h_e = H0_E_hs / 100.0
    # radiation-corrected Ω_m
    omr = OMEGA_R_H2 / h_e**2
    omm = 1.0 - OMEGA_DE_TODAY - omr
    wm = omm * h_e**2
    pull_wm = (wm - PLANCK_OMMH2) / PLANCK_OMMH2_SIGMA
    wm_data = (1.0 - OMEGA_DE_TODAY) * (H0_E_data / 100.0) ** 2
    pull_wm_d = (wm_data - PLANCK_OMMH2) / PLANCK_OMMH2_SIGMA
    print(f"  ω_m (H0_SN/f_∞, CMB frame)={wm:.5f}  {pull_wm:+.2f}σ")
    print(f"  ω_m (BAO H0_E, mixed frame)={wm_data:.5f}  {pull_wm_d:+.2f}σ")
    print(f"  local ω_m (H0_SN)={(1-OMEGA_DE_TODAY)*(H0_SN/100)**2:.5f}  (not CAMB input)")
    kill(kills, "K4", abs(pull_wm) >= 3.0, f"CMB-frame ω_m {pull_wm:+.2f}σ (kill if |σ|≥3)")

    # --- K5 K6 ages ---
    print("\n--- K5 / K6  baryon ages ---")
    t_L = age_gyr(H0_SN)
    t_E = age_gyr(H0_E_data)
    t_b = t_L * F_BAO
    t_b_cmb = t_L * F_CMB
    t_b_corr = t_L * F_CAMB
    print(f"  t_γ(H0_SN)={t_L:.3f}")
    print(f"  t_□=f_□ t_γ={t_b:.3f}  Valcin {(t_b-VALCIN_TU)/VALCIN_ERR:+.2f}σ  "
          f"t_GC {(t_b-VALCIN_TGC)/VALCIN_TGC_ERR:+.2f}σ")
    print(f"  t_∞=f_∞ t_γ={t_b_cmb:.3f}  Valcin {(t_b_cmb-VALCIN_TU)/VALCIN_ERR:+.2f}σ")
    print(f"  t_E=t_FRW(H0_E)={t_E:.3f}  Valcin {(t_E-VALCIN_TU)/VALCIN_ERR:+.2f}σ")
    k5 = max(abs(t_b - VALCIN_TU), abs(t_b_cmb - VALCIN_TU)) / VALCIN_ERR >= 3.0
    kill(kills, "K5", k5, f"t_□ {(t_b-VALCIN_TU)/VALCIN_ERR:+.2f}σ  t_∞ {(t_b_cmb-VALCIN_TU)/VALCIN_ERR:+.2f}σ")
    # clusters older than universe: t_b < t_GC at 2σ (t_b + 2σ_tgc < t_GC)
    older = (VALCIN_TGC - t_b) / VALCIN_TGC_ERR
    kill(kills, "K6", older >= 2.0, f"t_GC − t_b = {VALCIN_TGC-t_b:+.3f} Gyr ({older:+.2f}σ of t_GC)")

    # --- K7 SN lock cost ---
    print("\n--- K7  SN lock vs free Ω ---")
    dchi_sn = sn_lock["chi2"] - sn_free["chi2"]
    print(f"  lock χ²={sn_lock['chi2']:.2f}  free Ω={sn_free['omega_de']:.4f} χ²={sn_free['chi2']:.2f}  "
          f"Δ={dchi_sn:.2f} ({math.sqrt(max(dchi_sn,0)):.2f}σ)")
    kill(kills, "K7", dchi_sn >= 9.0, f"SN lock Δχ²={dchi_sn:.2f} (kill if ≥9)")

    # --- K8 time delays: ONE continuum ⇒ H0_TD = H0_L. Score PAST delay-only. ---
    print("\n--- K8  time delays (past delay-only; H0_TD = H0_L) ---")
    H0_td = H0_SN  # same metric, this octave
    td_samples = [
        ("H0LiCOW Wong+2020 6 lenses", 73.3, 1.75, True),
        ("TDCOSMO-1 Millon+2020 7 lenses", 74.2, 1.6, True),
        ("TDCOSMO-IV Birrer+2020 TD-only", 74.5, 5.85, True),
        ("TDCOSMO-2025 8 lenses", 72.1, 3.85, True),
        ("TDCOSMO-IV + SLACS (NOT delay-only)", 67.4, 3.65, False),
    ]
    td_kill = False
    td_rows = []
    for name, h, s, delay_only in td_samples:
        sig = tens(H0_td, h, 0.0, s)
        row = {"name": name, "H0": h, "sig_data": s, "pull": sig, "delay_only": delay_only}
        td_rows.append(row)
        tag = "delay-only" if delay_only else "mass-prior combo — not K8"
        print(f"  predict {H0_td:.2f} vs {name}: {h}±{s}  {sig:.2f}σ  [{tag}]")
        if delay_only and sig >= 3.0:
            td_kill = True
    kill(kills, "K8", td_kill, "any delay-only past sample ≥3σ from H0_L")

    # --- extra knives, not instant K ---
    print("\n--- extra (inform, do not auto-kill) ---")
    # DES is not SH0ES-calibrated; shape only
    des = SNLikelihood(
        load_des_dovekie(DATA / "DES-Dovekie_HD.csv", DATA / "DES-Dovekie_STAT+SYS.npz")
    )
    des_a = des.fit_h0(omega_de=OMEGA_DE_TODAY)
    des_c = des.fit_h0_omega_de(H0_guess=des_a["H0"], omega_guess=OMEGA_DE_TODAY)
    print(f"  DES lock Δχ² vs free Ω={des_a['chi2']-des_c['chi2']:.2f}  "
          f"Ω_free={des_c['omega_de']:.4f}")

    bao_planck_om = fit_h0_given_rd(
        PLANCK_RD_MPC, omega_de=PLANCK2018_OMEGA_DE, sample=DESI_DR2
    )
    sn_pl = like.fit_h0(omega_de=PLANCK2018_OMEGA_DE)
    print(
        f"  control Planck-Ω + f_□: SN {sn_pl['H0']:.3f} vs f_□×BAO {F_BAO*bao_planck_om['H0']:.3f}  "
        f"{tens(sn_pl['H0'], F_BAO*bao_planck_om['H0'], SIG_SN, SIG_BAO):.2f}σ"
    )

    R_lock_planckH = shift_parameter_R(PLANCK2018_H0, omega_de=OMEGA_DE_TODAY)
    print(f"  R lock @ Planck H0={R_lock_planckH:.4f}  {(R_lock_planckH-PLANCK_R)/PLANCK_R_SIGMA:+.2f}σ")

    n_fired = sum(1 for k in kills if k["fired"])
    print("\n" + "=" * 72)
    if n_fired:
        print(f"DISPROVED  {n_fired}/10 kill criteria fired")
    else:
        print("NOT DISPROVED  0/10 kill criteria fired")
        print("Survive ≠ proven. CAMB spectra, official DESI cov, and tighter")
        print("time-delay H0 remain outside this gauntlet.")
    print("=" * 72)

    out = {
        "kills": kills,
        "n_fired": n_fired,
        "splits": splits,
        "sn_lock": sn_lock,
        "sn_free": sn_free,
        "bao_early": bao_early,
        "chi_mixed": chi_mixed,
        "R_early": R_e,
        "omega_mh2_early": wm_data,
        "ages": {"t_gamma": t_L, "t_baryon_HS": t_b, "t_baryon_corr": t_b_corr, "t_Einstein": t_E},
        "H0_TD": H0_td,
        "td_history": td_rows,
        "survived": n_fired == 0,
    }
    path = RESULTS / "falsify.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"Wrote {path}")
    return 0 if n_fired == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

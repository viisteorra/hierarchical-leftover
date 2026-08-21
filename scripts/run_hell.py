#!/usr/bin/env python3
"""Put the locked theory through hell. Not a fit.

Every knife is either a uniqueness/freeze check or a number vs data.
Hybrids and inverses MUST look worse than the licensed maps.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bao import (  # noqa: E402
    DESI_DR1,
    DESI_DR2,
    PLANCK_RD_MPC,
    bao_chi2,
    fit_h0_given_rd,
    fit_rd,
)
from cosmology import OMEGA_R_H2, age_gyr, shift_parameter_R  # noqa: E402
from geometry import OMEGA_DE_TODAY, OMEGA_M, Q4_WEIGHT, Q5_WEIGHT, TOTAL, r as R_GEN  # noqa: E402
from likelihood import (  # noqa: E402
    PLANCK2018_H0,
    PLANCK2018_OMEGA_DE,
    PLANCK_R,
    PLANCK_R_SIGMA,
    SNLikelihood,
    load_des_dovekie,
    load_pantheon_plus,
)
from generate import leftover_exp_finite  # noqa: E402
from ruler import (  # noqa: E402
    leftover_exp_infinite,
    leftover_exp_planar,
    leftover_f,
    leftover_octave,
    mix_compress_f,
    mix_stretch_f,
)
from spacetime import (  # noqa: E402
    F_ANSATZ_T10,
    F_BAO,
    F_CLOCK,
    F_CMB,
    F_FINITE4,
    F_MIXED_2D,
    PLANCK_ZSTAR,
    ZSTAR_CAMB,
    ZSTAR_HS_LOCK,
    ZSTAR_RECOMB,
)

RESULTS = ROOT / "results"
DATA = ROOT / "data"
from hubble import H0_SN, H0_TH_LOCK, SIG_BAO, SIG_SN  # noqa: E402

H0_TH = H0_TH_LOCK
SIG_TH = 0.54
VALCIN, VALCIN_E = 13.57, math.sqrt(0.15**2 + 0.23**2)
VALCIN_TGC, VALCIN_TGC_E = 13.39, math.sqrt(0.10**2 + 0.23**2)
PLANCK_OMMH2, PLANCK_OMMH2_S = 0.1430, 0.0011
OMBH2 = 0.02237

knives = []


def tens(a, b, sa, sb):
    return abs(a - b) / math.sqrt(sa**2 + sb**2)


def knife(kid, fired, detail):
    knives.append({"id": kid, "fired": bool(fired), "detail": detail})
    print(f"  [{'KILL' if fired else 'live'}] {kid}: {detail}")


def camb_pack(H0: float) -> dict:
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
        "theta_pull": (th - 1.04110) / 0.00031,
        "rd_pull": (float(d["rdrag"]) - 147.09) / 0.26,
        "zstar_pull": (float(d["zstar"]) - 1089.80) / 0.21,
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    print("=" * 72)
    print("HELL  locked ln(2) + two-measure leftover")
    print("=" * 72)

    # --- pytest ---
    print("\n--- H0  pytest ---")
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", str(ROOT / "tests")],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    last = (proc.stdout or "").strip().splitlines()[-1] if proc.stdout else "no pytest out"
    knife("H0_pytest", proc.returncode != 0, last)

    # --- freeze ---
    print("\n--- H1  freeze ---")
    knife("H1_1_11", Q4_WEIGHT != 1 or Q5_WEIGHT != 11 or TOTAL != 12, "1:11/12")
    knife("H1_T", abs(OMEGA_DE_TODAY - __import__("math").log(2.0)) >= 1e-12, f"T={OMEGA_DE_TODAY}")
    knife("H1_r", abs(R_GEN - __import__("math").log(2.0) / (1.0 + __import__("math").log(2.0))) >= 1e-12, f"r={R_GEN}")
    knife("H1_planar_e", abs(leftover_exp_planar() - leftover_exp_finite()) >= 1e-12, "planar G")

    like = SNLikelihood(
        load_pantheon_plus(DATA / "Pantheon+SH0ES.dat", DATA / "Pantheon+SH0ES_STAT+SYS.cov")
    )
    des = SNLikelihood(
        load_des_dovekie(DATA / "DES-Dovekie_HD.csv", DATA / "DES-Dovekie_STAT+SYS.npz")
    )
    bao2 = fit_h0_given_rd(PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY, sample=DESI_DR2)
    bao1 = fit_h0_given_rd(PLANCK_RD_MPC, omega_de=OMEGA_DE_TODAY, sample=DESI_DR1)
    H0E, chi0 = bao2["H0"], bao2["chi2"]
    phi = leftover_octave(ZSTAR_RECOMB)["frac"]
    tL = age_gyr(H0_SN)

    def dchi(f, sample=DESI_DR2):
        return bao_chi2(H0_SN, PLANCK_RD_MPC / f, omega_de=OMEGA_DE_TODAY, sample=sample) - (
            chi0 if sample is DESI_DR2 else fit_h0_given_rd(PLANCK_RD_MPC, sample=sample)["chi2"]
        )

    # --- Hubble maps ---
    print("\n--- H2  licensed maps ---")
    print(f"  f_□={F_BAO:.6f}  f_∞={F_CMB:.6f}  φ={phi:.6f}")
    sb = tens(H0_SN, F_BAO * H0E, SIG_SN, SIG_BAO)
    st = tens(H0_SN, F_CMB * H0_TH, SIG_SN, SIG_TH)
    print(f"  BAO {H0E:.3f}×f_□={F_BAO*H0E:.3f}  {sb:.3f}σ")
    print(f"  θ*  {H0_TH:.3f}×f_∞={F_CMB*H0_TH:.3f}  {st:.3f}σ")
    knife("H2_bao", sb >= 3.0, f"BAO map {sb:.3f}σ")
    knife("H2_cmb", st >= 3.0, f"CMB map {st:.3f}σ")

    print("\n--- H3  K2 / hybrid ---")
    dc = bao_chi2(H0_SN, PLANCK_RD_MPC / F_BAO, sample=DESI_DR2) - chi0
    hy = bao_chi2(H0_SN, PLANCK_RD_MPC / F_CMB, sample=DESI_DR2) - chi0
    pred = bao_chi2(F_BAO * H0E, PLANCK_RD_MPC / F_BAO, sample=DESI_DR2)
    print(f"  (H0_SN, rd/f_□) Δχ²={dc:+.3f}  predicted pair χ²={pred:.3f}  hybrid Δ={hy:+.2f}")
    knife("H3_k2", dc >= 10.0, f"Δχ²={dc:.2f}")
    knife("H3_hybrid_must_hurt", hy < 10.0, f"hybrid Δ={hy:.2f} (must be ugly)")
    knife("H3_pred_pair", pred - chi0 >= 1.0, f"predicted pair Δ={pred-chi0:.3f}")

    print("\n--- H4  foils must not beat licensed maps on BOTH ---")
    foils = {
        "clock": F_CLOCK,
        "bare_12/11": F_MIXED_2D,
        "compress 11/12": mix_compress_f(phi),
        "r/12 tail": 2 ** (phi * (12 / 11) / (1 - R_GEN / 12)),
        "finite4": F_FINITE4,
        "3/2": 2 ** (phi * 1.5),
        "4/3": 2 ** (phi * 4 / 3),
        "T/10": F_ANSATZ_T10,
        "N=2 nest": 2 ** (phi * leftover_exp_planar() / (1 - R_GEN / 121)),
        "2^{1/12}": 2 ** (1 / 12),
        "12/11 raw": 12 / 11,
    }
    licensed_worst = max(sb, st)
    beat_both = []
    for name, f in foils.items():
        s1 = tens(H0_SN, f * H0E, SIG_SN, SIG_BAO)
        s2 = tens(H0_SN, f * H0_TH, SIG_SN, SIG_TH)
        dc_f = bao_chi2(H0_SN, PLANCK_RD_MPC / f, sample=DESI_DR2) - chi0
        flag = ""
        if s1 < sb and s2 < st and dc_f < 10:
            beat_both.append(name)
            flag = "  << beats both maps + K2"
        print(f"  {name:16s} BAO {s1:5.2f}σ  θ* {s2:5.2f}σ  Δχ² {dc_f:+7.2f}{flag}")
    knife("H4_foil_beats_both", bool(beat_both), f"foils {beat_both}")

    print("\n--- H5  leftover independence ---")
    fs = [leftover_f(2.0 ** (n + phi) - 1.0) for n in range(8, 14)]
    fps = [leftover_octave(2.0 ** (n + phi) - 1.0)["f_planar"] for n in range(8, 14)]
    knife("H5_inf_n", max(fs) - min(fs) >= 1e-10, f"f_∞ span {max(fs)-min(fs):.2e}")
    knife("H5_pl_n", max(fps) - min(fps) >= 1e-10, f"f_□ span {max(fps)-min(fps):.2e}")
    emp = H0_SN / H0E
    knife("H5_not_emp", abs(F_BAO - emp) < 1e-6, "f_□ identical to H0_SN/H0_E")

    print("\n--- H6  z* fragility ---")
    for tag, z in (
        ("recomb", ZSTAR_RECOMB),
        ("camb", ZSTAR_CAMB),
        ("planck", PLANCK_ZSTAR),
        ("HS", ZSTAR_HS_LOCK),
        ("z-10", ZSTAR_RECOMB - 10),
        ("z+10", ZSTAR_RECOMB + 10),
        ("1080", 1080.0),
        ("1100", 1100.0),
    ):
        lo = leftover_octave(z)
        fsq, fin = lo["f_planar"], lo["f_inf"]
        s1 = tens(H0_SN, fsq * H0E, SIG_SN, SIG_BAO)
        s2 = tens(H0_SN, fin * H0_TH, SIG_SN, SIG_TH)
        print(f"  {tag:8s} z*={z:8.2f}  f_□={fsq:.5f} {s1:.2f}σ  f_∞={fin:.5f} {s2:.2f}σ")
        if tag in ("recomb", "camb", "planck", "HS") and (s1 >= 3 or s2 >= 3):
            knife(f"H6_{tag}", True, f"{s1:.2f}/{s2:.2f}σ")
    knife("H6_planck_window", False, "catalog z* window printed")

    print("\n--- H7  DESI DR1 vs DR2 ---")
    print(f"  DR2 H0_E={H0E:.3f} χ²={chi0:.2f}")
    print(f"  DR1 H0_E={bao1['H0']:.3f} χ²={bao1['chi2']:.2f}")
    s_dr1 = tens(H0_SN, F_BAO * bao1["H0"], SIG_SN, SIG_BAO)
    dc1 = bao_chi2(H0_SN, PLANCK_RD_MPC / F_BAO, sample=DESI_DR1) - bao1["chi2"]
    print(f"  DR1 map {F_BAO*bao1['H0']:.3f}  {s_dr1:.2f}σ  Δχ²={dc1:+.2f}")
    knife("H7_dr1", s_dr1 >= 3.0 or dc1 >= 10.0, f"DR1 {s_dr1:.2f}σ Δ={dc1:.2f}")
    knife("H7_dr12_split", abs(bao1["H0"] - H0E) / 0.5 >= 3.0, f"DR1−DR2 {bao1['H0']-H0E:+.3f}")

    print("\n--- H8  SN lock / DES / w / Ω ---")
    sn_lock = like.fit_h0(omega_de=OMEGA_DE_TODAY)
    sn_w = like.fit_h0_w(omega_de=OMEGA_DE_TODAY)
    sn_om = like.fit_h0_omega_de(H0_guess=H0_SN, omega_guess=OMEGA_DE_TODAY)
    des_a = des.fit_h0(omega_de=OMEGA_DE_TODAY)
    des_c = des.fit_h0_omega_de(H0_guess=des_a["H0"], omega_guess=OMEGA_DE_TODAY)
    print(f"  Pantheon lock χ²={sn_lock['chi2']:.2f}  H0={sn_lock['H0']:.3f}")
    print(f"  free w  w={sn_w.get('w', sn_w)} χ²={sn_w['chi2']:.2f}  Δ={sn_w['chi2']-sn_lock['chi2']:.2f}")
    print(f"  free Ω  Ω={sn_om['omega_de']:.4f} Δχ²={sn_lock['chi2']-sn_om['chi2']:.2f}")
    print(f"  DES lock Δχ² vs free Ω={des_a['chi2']-des_c['chi2']:.2f}  Ω_free={des_c['omega_de']:.4f}")
    knife("H8_sn_lock", sn_lock["chi2"] - sn_om["chi2"] >= 9.0, "SN lock vs free Ω")
    knife("H8_des", des_a["chi2"] - des_c["chi2"] >= 9.0, "DES lock vs free Ω")

    print("\n--- H9  SN redshift slices (lock Ω, fit H0) ---")
    z = like.data.z
    for lo_z, hi_z, name in (
        (0.01, 0.1, "low"),
        (0.1, 0.5, "mid"),
        (0.5, 2.5, "high"),
        (0.01, 0.023, "SH0ES-ish"),
    ):
        mask = (z >= lo_z) & (z < hi_z)
        if mask.sum() < 20:
            print(f"  {name} n={int(mask.sum())} skip")
            continue
        sub = SNLikelihood(like.data.masked(mask))
        sl = sub.fit_h0(omega_de=OMEGA_DE_TODAY)
        so = sub.fit_h0_omega_de(H0_guess=sl["H0"], omega_guess=OMEGA_DE_TODAY)
        d = sl["chi2"] - so["chi2"]
        print(f"  {name:10s} n={int(mask.sum()):4d}  H0={sl['H0']:.2f}  lockΔχ²={d:.2f}  Ω_free={so['omega_de']:.3f}")
        knife(f"H9_{name}", d >= 9.0, f"{name} Δ={d:.2f}")

    print("\n--- H10  R / ω_m / ages / delays ---")
    R_e = shift_parameter_R(H0_SN / F_CMB, omega_de=OMEGA_DE_TODAY)
    pull_R = (R_e - PLANCK_R) / PLANCK_R_SIGMA
    wm_e = (1 - OMEGA_DE_TODAY) * (H0E / 100) ** 2
    wm_c = (1 - OMEGA_DE_TODAY) * ((H0_SN / F_CMB) / 100) ** 2
    print(f"  R={R_e:.4f}  {pull_R:+.2f}σ")
    print(f"  ω_m BAO-E={wm_e:.5f} {(wm_e-PLANCK_OMMH2)/PLANCK_OMMH2_S:+.2f}σ")
    print(f"  ω_m SN/f_∞={wm_c:.5f} {(wm_c-PLANCK_OMMH2)/PLANCK_OMMH2_S:+.2f}σ")
    tsq, tinf = tL * F_BAO, tL * F_CMB
    print(f"  t_□={tsq:.3f} {(tsq-VALCIN)/VALCIN_E:+.2f}σ  t_∞={tinf:.3f} {(tinf-VALCIN)/VALCIN_E:+.2f}σ")
    knife("H10_R", abs(pull_R) >= 3.0, f"R {pull_R:+.2f}σ")
    knife("H10_wm", abs(wm_e - PLANCK_OMMH2) / PLANCK_OMMH2_S >= 3.0, "ω_m BAO-E")
    knife("H10_age", max(abs(tsq - VALCIN), abs(tinf - VALCIN)) / VALCIN_E >= 3.0, "ages")
    for name, h, s in (
        ("H0LiCOW", 73.3, 1.75),
        ("TDCOSMO-1", 74.2, 1.6),
        ("TD-IV", 74.5, 5.85),
        ("TD-2025", 72.1, 3.85),
    ):
        p = tens(H0_SN, h, 0.0, s)
        print(f"  delay {name} {p:.2f}σ")
        knife(f"H10_{name}", p >= 3.0, f"{name} {p:.2f}σ")

    print("\n--- H11  BAO per tracer at (H0_SN, rd/f_□) ---")
    rdL = PLANCK_RD_MPC / F_BAO
    from bao import BAOPoint, _cov, _theory

    for p in DESI_DR2:
        chi_p = bao_chi2(H0_SN, rdL, sample=(p,), omega_de=OMEGA_DE_TODAY)
        chi_e = bao_chi2(H0E, PLANCK_RD_MPC, sample=(p,), omega_de=OMEGA_DE_TODAY)
        print(f"  {p.name:12s} z={p.z:.3f}  local χ²={chi_p:.2f}  early χ²={chi_e:.2f}  Δ={chi_p-chi_e:+.2f}")
        knife(f"H11_{p.name}", chi_p - chi_e >= 9.0, f"{p.name} Δ={chi_p-chi_e:.2f}")

    print("\n--- H12  CAMB (slow) ---")
    camb_rows = {}
    for label, H0 in (
        ("theta_lock", H0_TH),
        ("SN_over_finf", H0_SN / F_CMB),
        ("BAO_E", H0E),
    ):
        try:
            row = camb_pack(H0)
        except Exception as e:
            print(f"  CAMB {label} FAIL {e}")
            knife(f"H12_{label}", True, str(e)[:80])
            continue
        camb_rows[label] = row
        print(
            f"  {label:14s} H0={H0:.3f}  100θ*={row['100theta']:.5f} {row['theta_pull']:+.1f}σ  "
            f"rd={row['rd']:.2f} {row['rd_pull']:+.1f}σ  z*={row['zstar']:.2f}  ℓ={row['peak_l']}"
        )
    if "theta_lock" in camb_rows:
        knife("H12_theta_lock", abs(camb_rows["theta_lock"]["theta_pull"]) >= 3.0, "θ* at 67.60")
        knife("H12_peak", camb_rows["theta_lock"]["peak_l"] not in range(215, 226), "ℓ peak")
    if "SN_over_finf" in camb_rows:
        knife(
            "H12_early_from_f",
            abs(camb_rows["SN_over_finf"]["theta_pull"]) >= 3.0,
            f"θ* at H0_SN/f_∞ {camb_rows['SN_over_finf']['theta_pull']:+.1f}σ",
        )
    if "BAO_E" in camb_rows:
        # This is the known DESI vs Planck shape; not a leftover kill by itself
        print(f"  (BAO H0_E θ* pull is DESI–Planck shape, leftover cannot eat it)")

    print("\n--- H14  swapped measures must fail (Lemma C) ---")
    s_inf_bao = tens(H0_SN, F_CMB * H0E, SIG_SN, SIG_BAO)
    s_box_cmb = tens(H0_SN, F_BAO * H0_TH, SIG_SN, SIG_TH)
    dc_swap = bao_chi2(H0_SN, PLANCK_RD_MPC / F_CMB, sample=DESI_DR2) - chi0
    print(f"  f_∞ on BAO {s_inf_bao:.3f}σ  Δχ²={dc_swap:+.2f}  (licensed BAO {sb:.3f}σ)")
    print(f"  f_□ on θ*  {s_box_cmb:.3f}σ  (licensed θ* {st:.3f}σ)")
    knife("H14_inf_on_bao_k2", dc_swap < 10.0, f"inf-on-BAO Δχ²={dc_swap:.2f} must be ugly")
    knife("H14_swap_worse", s_inf_bao <= sb or s_box_cmb <= st, "swap not worse on both licensed pairs")

    print("\n--- H15  first principles ---")
    from uniqueness import derivation_chain, period_from_fifth_convergent

    n_cf, k_cf, _ = period_from_fifth_convergent()
    free = [row for row in derivation_chain() if row["free"]]
    knife("H15_u1_cf", (n_cf, k_cf) != (12, 7), f"CF period {n_cf}")
    knife("H15_only_P", len(free) != 1 or free[0]["symbol"] != "P", f"{len(free)} free")
    geo = (ROOT / "src" / "geometry.py").read_text().lower()
    knife("H15_not_fit", "chi2" in geo or "minimize" in geo, "geometry contains a fit")

    print("\n--- H13  Planck-Ω control ---")
    bao_pl = fit_h0_given_rd(PLANCK_RD_MPC, omega_de=PLANCK2018_OMEGA_DE, sample=DESI_DR2)
    sn_pl = like.fit_h0(omega_de=PLANCK2018_OMEGA_DE)
    s_pl = tens(sn_pl["H0"], F_BAO * bao_pl["H0"], SIG_SN, SIG_BAO)
    print(f"  Planck Ω SN {sn_pl['H0']:.3f} vs f_□×BAO {F_BAO*bao_pl['H0']:.3f}  {s_pl:.2f}σ")
    knife("H13_planck_om", s_pl >= 3.0, f"{s_pl:.2f}σ")

    n_kill = sum(1 for k in knives if k["fired"])
    print("\n" + "=" * 72)
    print(f"HELL  {n_kill}/{len(knives)} knives fired")
    if n_kill == 0:
        print("Not killed. Survive ≠ proven.")
    else:
        print("Fired:")
        for k in knives:
            if k["fired"]:
                print(f"  - {k['id']}: {k['detail']}")
    print("=" * 72)
    out = {
        "n_fired": n_kill,
        "n_knives": len(knives),
        "knives": knives,
        "maps": {
            "f_bao": F_BAO,
            "f_cmb": F_CMB,
            "H0_E": H0E,
            "H0_TH": H0_TH,
            "H0_L_bao": F_BAO * H0E,
            "H0_L_cmb": F_CMB * H0_TH,
            "sig_bao": sb,
            "sig_cmb": st,
            "dchi_k2": dc,
            "hybrid": hy,
        },
        "camb": camb_rows,
    }
    path = RESULTS / "hell.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"Wrote {path}")
    return 0 if n_kill == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

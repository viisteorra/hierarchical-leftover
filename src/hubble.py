"""Hubble scale. Leftover f is not a function of H0.

Theorem 3 (docs/THEOREM.md): two readings of one FRW, leftover maps
the scale. SN tests both. The early–early H0 ratio is the measure ratio.

    H0_L = f_□ H0_E,BAO = f_∞ H0_E,θ*
    H0_E,BAO / H0_E,θ* = f_∞ / f_□

No χ² in f. No T/10. Density is ln(2).
"""

from __future__ import annotations

import math

from spacetime import F_BAO, F_CMB

# Hubble-tension error budget: SH0ES calibration, not intercept scatter.
SIG_SN = 1.04
SIG_BAO = 0.190  # official DESI DR2 ALL_GCcomb, Δχ²=1 at lock
SIG_PLANCK_H0 = 0.54
PLANCK2018_H0 = 67.36
PLANCK_100THETA = 1.04110
PLANCK_100THETA_SIG = 0.00031
PLANCK_RD = 147.09
PLANCK_RD_SIG = 0.26
PLANCK_OMMH2 = 0.1430
PLANCK_OMMH2_SIG = 0.0011

# Frozen intercepts of the continuous lock. Scripts recompute; tests use
# these so pytest does not refit Pantheon+ or CAMB.
H0_SN = 73.50430785786638
H0_BAO = 68.53650184574981  # official DESI DR2 cov + Planck rd
H0_TH_LOCK = 67.79577217102053  # CAMB bisection at ln(2) → 100θ* = 1.04110


def tension(a: float, b: float, sa: float, sb: float) -> float:
    return abs(a - b) / math.sqrt(sa**2 + sb**2)


def solution(
    h0_sn: float = H0_SN,
    h0_bao: float = H0_BAO,
    h0_th: float = H0_TH_LOCK,
) -> dict:
    """Scale accounting. f comes from leftover, not from these H0s."""
    f_box = float(F_BAO)
    f_inf = float(F_CMB)
    pred_sn_from_bao = f_box * h0_bao
    pred_sn_from_th = f_inf * h0_th
    pred_bao_from_sn = h0_sn / f_box
    pred_th_from_sn = h0_sn / f_inf
    ratio_f = f_inf / f_box
    ratio_h = h0_bao / h0_th
    pred_bao_from_th = h0_th * ratio_f
    pred_th_from_bao = h0_bao / ratio_f
    return {
        "f_box": f_box,
        "f_inf": f_inf,
        "h0_sn": h0_sn,
        "h0_bao": h0_bao,
        "h0_th": h0_th,
        "pred_sn_from_bao": pred_sn_from_bao,
        "pred_sn_from_th": pred_sn_from_th,
        "pred_bao_from_sn": pred_bao_from_sn,
        "pred_th_from_sn": pred_th_from_sn,
        "ratio_f": ratio_f,
        "ratio_h": ratio_h,
        "pred_bao_from_th": pred_bao_from_th,
        "pred_th_from_bao": pred_th_from_bao,
        "raw_bao_sigma": tension(h0_sn, h0_bao, SIG_SN, SIG_BAO),
        "raw_th_sigma": tension(h0_sn, h0_th, SIG_SN, SIG_PLANCK_H0),
        "raw_planck_sigma": tension(h0_sn, PLANCK2018_H0, SIG_SN, SIG_PLANCK_H0),
        "map_bao_sigma": tension(h0_sn, pred_sn_from_bao, SIG_SN, SIG_BAO),
        "map_th_sigma": tension(h0_sn, pred_sn_from_th, SIG_SN, SIG_PLANCK_H0),
        "fwd_bao_sigma": tension(h0_bao, pred_bao_from_sn, SIG_BAO, 0.0),
        "early_early_sigma": tension(h0_bao, pred_bao_from_th, SIG_BAO, SIG_PLANCK_H0),
        "split_bao_kms": h0_sn - pred_sn_from_bao,
        "split_th_kms": h0_sn - pred_sn_from_th,
        "raw_split_bao_kms": h0_sn - h0_bao,
        "raw_split_th_kms": h0_sn - h0_th,
    }


def omega_m_h2(h0: float, omega_m: float) -> float:
    return omega_m * (h0 / 100.0) ** 2

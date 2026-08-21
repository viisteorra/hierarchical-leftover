#!/usr/bin/env python3
"""CAMB TT peaks at the early reading.

Early H0 is the θ* lock 67.5991 on Ω_m=22/71, not BAO 68.38 (mixed frame).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from geometry import OMEGA_M  # noqa: E402
from spacetime import F_CMB  # noqa: E402

RESULTS = ROOT / "results"
H0_TH = 67.5991
H0_PL = 67.36
OMBH2 = 0.02237


def peaks(tt, n=4):
    ell = np.arange(tt.size)
    out = []
    start = 80
    for _ in range(n):
        if start + 50 >= 900:
            break
        sl = tt[start : start + 220]
        i = int(np.argmax(sl))
        L = start + i
        out.append({"ell": L, "CTT": float(tt[L])})
        start = L + 80
    return out


def run(H0, om_m):
    import camb

    h = H0 / 100.0
    omch2 = om_m * h**2 - OMBH2
    p = camb.CAMBparams()
    p.set_cosmology(H0=H0, ombh2=OMBH2, omch2=omch2, omk=0.0, tau=0.0544)
    p.InitPower.set_params(As=1e-10 * np.exp(3.044), ns=0.9649)
    p.set_for_lmax(2500, lens_potential_accuracy=1)
    r = camb.get_results(p)
    d = r.get_derived_params()
    cls = r.get_lensed_scalar_cls(lmax=2500)
    tt = cls[:2501, 0]
    ts = float(d["thetastar"])
    th = 100.0 * ts if ts < 1 else ts
    return {
        "H0": H0,
        "om_m": om_m,
        "100theta": th,
        "rd": float(d["rdrag"]),
        "zstar": float(d["zstar"]),
        "age": float(d["age"]),
        "peaks": peaks(tt),
        "C220": float(tt[220]),
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    print("=" * 72)
    print("CAMB spectra  early reading")
    print("=" * 72)
    print(f"  lock Ω_m={OMEGA_M:.6f}=22/71  H0_θ={H0_TH}  f_∞={F_CMB:.6f}")
    lock = run(H0_TH, OMEGA_M)
    fid = run(H0_PL, 0.3153)
    print(f"  lock  100θ*={lock['100theta']:.5f}  rd={lock['rd']:.2f}  z*={lock['zstar']:.2f}")
    print(f"  fid   100θ*={fid['100theta']:.5f}  rd={fid['rd']:.2f}  z*={fid['zstar']:.2f}")
    print("  TT peaks (ℓ):")
    for a, b in zip(lock["peaks"], fid["peaks"]):
        dC = (a["CTT"] - b["CTT"]) / b["CTT"] if b["CTT"] else 0.0
        print(f"    lock ℓ={a['ell']}  fid ℓ={b['ell']}  dC/C={dC:+.3%}")
    print(f"  C_220 lock/fid {(lock['C220']/fid['C220']-1):+.3%}")
    path = RESULTS / "spectra.json"
    path.write_text(json.dumps({"lock": lock, "fiducial_planck_om": fid}, indent=2, default=float))
    print(f"Wrote {path}")
    # knives: first peak 220±5, |dC/C| at ℓ=220 not huge
    p0 = lock["peaks"][0]["ell"]
    ok_peak = 215 <= p0 <= 225
    print(f"  first peak {'OK' if ok_peak else 'OFF'} ℓ={p0}")
    return 0 if ok_peak else 1


if __name__ == "__main__":
    raise SystemExit(main())

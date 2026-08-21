#!/usr/bin/env python3
"""Print the frozen solution."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from geometry import LN2, OMEGA_DE_TODAY, T_DISCRETE, W  # noqa: E402
from hubble import solution  # noqa: E402
from spacetime import ZSTAR_RECOMB  # noqa: E402
from uniqueness import PERIOD_N, SEED, HYPER_STEPS  # noqa: E402


def main() -> int:
    s = solution()
    print("SOLUTION  given P  continuous")
    print(
        f"  Ω_DE=ln2={OMEGA_DE_TODAY:.6f}  discrete 12-fold {float(T_DISCRETE):.6f}=49/71  "
        f"period {PERIOD_N} seed {SEED}:{HYPER_STEPS}  w={W}"
    )
    print(f"  z*={ZSTAR_RECOMB}  f_□={s['f_box']:.6f}  f_∞={s['f_inf']:.6f}")
    print(
        f"  BAO  {s['h0_bao']:.3f} × f_□ = {s['pred_sn_from_bao']:.3f}   "
        f"SN {s['h0_sn']:.3f}  {s['map_bao_sigma']:.3f}σ  (raw {s['raw_bao_sigma']:.2f}σ)"
    )
    print(
        f"  θ*   {s['h0_th']:.3f} × f_∞ = {s['pred_sn_from_th']:.3f}   "
        f"SN {s['h0_sn']:.3f}  {s['map_th_sigma']:.3f}σ  (raw {s['raw_th_sigma']:.2f}σ)"
    )
    print(
        f"  early–early  f_∞/f_□={s['ratio_f']:.6f}  DESI/θ*={s['ratio_h']:.6f}  "
        f"{s['early_early_sigma']:.3f}σ"
    )
    print(f"  forward  SN/f_□={s['pred_bao_from_sn']:.3f}  SN/f_∞={s['pred_th_from_sn']:.4f}")
    print("  see docs/THEORY.md  docs/CLAIM.md  docs/COMPLETE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

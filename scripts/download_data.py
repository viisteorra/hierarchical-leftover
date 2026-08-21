#!/usr/bin/env python3
"""Download Pantheon+ SH0ES Hubble-diagram distances and STAT+SYS covariance."""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"

DAT_NAME = "Pantheon+SH0ES.dat"
COV_NAME = "Pantheon+SH0ES_STAT+SYS.cov"

# Official data release (Brout+2022 / Scolnic+2022 / Riess+2022).
BASE = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/"
    "main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/"
)
DAT_URL = BASE + "Pantheon%2BSH0ES.dat"
COV_URL = BASE + "Pantheon%2BSH0ES_STAT%2BSYS.cov"

DES_HD_NAME = "DES-Dovekie_HD.csv"
DES_NPZ_NAME = "DES-Dovekie_STAT+SYS.npz"
DES_HD_URL = (
    "https://raw.githubusercontent.com/des-science/DES-SN5YR/main/"
    "4_DISTANCES_COVMAT/DES-Dovekie_HD.csv"
)
DES_NPZ_URL = (
    "https://raw.githubusercontent.com/des-science/DES-SN5YR/main/"
    "4_DISTANCES_COVMAT/STAT%2BSYS.npz"
)

USER_AGENT = "hierarchical-de/0.1 (Pantheon+ cosmology smoke test)"


def _download(url: str, dest: Path, min_bytes: int) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"already have {dest} ({dest.stat().st_size} bytes)")
        return
    print(f"downloading {url}")
    print(f"         -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=120) as resp, tmp.open("wb") as fh:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    size = tmp.stat().st_size
    if size < min_bytes:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{dest.name} too small ({size} bytes); expected >= {min_bytes}")
    tmp.replace(dest)
    print(f"saved {dest} ({size} bytes)")


def main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    _download(DAT_URL, DATA_DIR / DAT_NAME, min_bytes=100_000)
    _download(COV_URL, DATA_DIR / COV_NAME, min_bytes=1_000_000)
    _download(DES_HD_URL, DATA_DIR / DES_HD_NAME, min_bytes=50_000)
    _download(DES_NPZ_URL, DATA_DIR / DES_NPZ_NAME, min_bytes=1_000_000)
    print("Pantheon+ and DES-SN5YR data ready.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

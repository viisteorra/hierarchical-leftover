"""Pantheon+ Hubble-diagram χ² with geometrically locked Ω_DE.

Observed distances are the SH0ES-calibrated moduli MU_SH0ES, so a fit
for H0 is a real intercept fit on a calibrated ladder rather than a
pure shape-only SN test. The covariance is the official STAT+SYS
matrix; diagonal-only errors are never used for cosmological χ².
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize, minimize_scalar

from cosmology import distance_modulus
from geometry import OMEGA_DE_TODAY
from stretch import stretch_dmu

# Pantheon+ cosmology papers cut z < 0.01 (peculiar-velocity dominated).
DEFAULT_ZMIN = 0.01

# Named comparison values. Not geometry, never written back into geometry.py.
PLANCK_OMEGA_DE = 0.689  # design-doc Planck-era figure (≈ Planck+BAO)
PLANCK2018_OMEGA_M = 0.3153  # TT,TE,EE+lowE+lensing
PLANCK2018_OMEGA_M_SIGMA = 0.0073
PLANCK2018_OMEGA_DE = 1.0 - PLANCK2018_OMEGA_M  # 0.6847
PLANCK2018_OMEGA_DE_SIGMA = PLANCK2018_OMEGA_M_SIGMA
PLANCK2018_H0 = 67.36
PLANCK_R = 1.7502  # Chen, Huang, Wang 2019; TT,TE,EE+lowE
PLANCK_R_SIGMA = 0.0046


@dataclass
class SNData:
    z: np.ndarray
    mu: np.ndarray
    cov: np.ndarray
    zmin: float
    n_raw: int
    path_dat: Path
    path_cov: Path

    @property
    def n(self) -> int:
        return int(self.z.size)

    def masked(self, mask: np.ndarray) -> "SNData":
        idx = np.flatnonzero(np.asarray(mask, dtype=bool))
        return SNData(
            z=self.z[idx],
            mu=self.mu[idx],
            cov=self.cov[np.ix_(idx, idx)],
            zmin=self.zmin,
            n_raw=self.n_raw,
            path_dat=self.path_dat,
            path_cov=self.path_cov,
        )

    def with_mu(self, mu: np.ndarray) -> "SNData":
        mu = np.asarray(mu, dtype=float)
        if mu.shape != self.mu.shape:
            raise ValueError("mu shape mismatch")
        return SNData(
            z=self.z,
            mu=mu,
            cov=self.cov,
            zmin=self.zmin,
            n_raw=self.n_raw,
            path_dat=self.path_dat,
            path_cov=self.path_cov,
        )


def load_pantheon_plus(
    dat_path: Path,
    cov_path: Path,
    zmin: float = DEFAULT_ZMIN,
) -> SNData:
    """Load Pantheon+SH0ES Hubble-diagram distances and STAT+SYS covariance."""
    dat_path = Path(dat_path)
    cov_path = Path(cov_path)
    df = pd.read_csv(dat_path, sep=r"\s+")
    for col in ("zHD", "MU_SH0ES"):
        if col not in df.columns:
            raise ValueError(f"{dat_path} missing column {col}; got {list(df.columns)}")

    n_raw = len(df)
    cov = _read_cov(cov_path, n_raw)

    z_all = df["zHD"].to_numpy(dtype=float)
    mu_all = df["MU_SH0ES"].to_numpy(dtype=float)
    mask = z_all > zmin
    if not np.any(mask):
        raise ValueError(f"no supernovae with zHD > {zmin}")
    idx = np.flatnonzero(mask)

    return SNData(
        z=z_all[idx],
        mu=mu_all[idx],
        cov=cov[np.ix_(idx, idx)],
        zmin=zmin,
        n_raw=n_raw,
        path_dat=dat_path,
        path_cov=cov_path,
    )


def load_des_dovekie(hd_path: Path, npz_path: Path) -> SNData:
    """Load DES-SN5YR / Dovekie Hubble diagram (MU at fiducial H0=70).

    Covariance is stored as a packed upper-triangular *inverse* matrix
    (DES-SN5YR README). We invert once so SNLikelihood can use Cholesky.
    """
    hd_path = Path(hd_path)
    npz_path = Path(npz_path)
    z = []
    mu = []
    for line in hd_path.read_text().splitlines():
        if not line.startswith("SN:"):
            continue
        parts = line.split()
        z.append(float(parts[3]))
        mu.append(float(parts[5]))
    z = np.asarray(z, dtype=float)
    mu = np.asarray(mu, dtype=float)
    packed = np.load(npz_path)
    n = int(packed["nsn"][0])
    if n != z.size:
        raise ValueError(f"DES nsn={n} != Hubble-diagram rows {z.size}")
    inv = np.zeros((n, n), dtype=float)
    inv[np.triu_indices(n)] = packed["cov"]
    inv = inv + np.tril(inv.T, -1)
    cov = np.linalg.inv(inv)
    return SNData(
        z=z,
        mu=mu,
        cov=cov,
        zmin=float(z.min()),
        n_raw=n,
        path_dat=hd_path,
        path_cov=npz_path,
    )


def _read_cov(path: Path, n: int) -> np.ndarray:
    raw = np.loadtxt(path)
    if raw.size == n * n + 1:
        n_cov = int(raw.flat[0])
        if n_cov != n:
            raise ValueError(f"covariance N={n_cov} does not match catalogue N={n}")
        cov = raw.flat[1:].reshape(n, n)
    elif raw.size == n * n:
        cov = raw.reshape(n, n)
    else:
        raise ValueError(
            f"covariance size {raw.size} is neither {n*n} nor {n*n + 1}"
        )
    if cov.shape != (n, n):
        raise ValueError(f"covariance shape {cov.shape} != {(n, n)}")
    return cov


class SNLikelihood:
    """Gaussian SN Ia likelihood with a pre-factorized covariance."""

    def __init__(self, data: SNData):
        self.data = data
        self._cho = cho_factor(data.cov, lower=True, check_finite=False)

    def chi2(
        self,
        H0: float,
        w: float = -1.0,
        omega_de: float = OMEGA_DE_TODAY,
        delta0: float = 0.0,
    ) -> float:
        mu_model = distance_modulus(self.data.z, H0, w=w, omega_de=omega_de)
        if delta0 != 0.0:
            mu_model = mu_model + stretch_dmu(self.data.z, delta0)
        residual = self.data.mu - mu_model
        solved = cho_solve(self._cho, residual, check_finite=False)
        return float(residual @ solved)

    def fit_h0(
        self,
        w: float = -1.0,
        omega_de: float = OMEGA_DE_TODAY,
        delta0: float = 0.0,
        H0_bounds: tuple[float, float] = (50.0, 100.0),
    ) -> dict:
        """Lock w and Ω_DE, fit only H0. Optional SN stretch δ0 is held fixed."""

        def objective(H0: float) -> float:
            return self.chi2(float(H0), w=w, omega_de=omega_de, delta0=delta0)

        result = minimize_scalar(objective, bounds=H0_bounds, method="bounded")
        H0 = float(result.x)
        chi2 = float(result.fun)
        n = self.data.n
        npar = 1
        return {
            "H0": H0,
            "w": w,
            "omega_de": omega_de,
            "delta0": delta0,
            "chi2": chi2,
            "n": n,
            "npar": npar,
            "dof": n - npar,
            "chi2_dof": chi2 / (n - npar),
            "success": bool(result.success),
        }

    def fit_h0_delta0(
        self,
        omega_de: float = OMEGA_DE_TODAY,
        H0_bounds: tuple[float, float] = (50.0, 100.0),
        delta_bounds: tuple[float, float] = (-0.2, 0.2),
        H0_guess: float = 70.0,
        delta_guess: float = 0.0,
    ) -> dict:
        """Stage B: lock Ω_DE and w=-1, stretch SN distances only."""

        def objective(params) -> float:
            H0, delta0 = float(params[0]), float(params[1])
            return self.chi2(H0, w=-1.0, omega_de=omega_de, delta0=delta0)

        result = minimize(
            objective,
            x0=np.array([H0_guess, delta_guess]),
            method="Nelder-Mead",
            bounds=[H0_bounds, delta_bounds],
            options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 600},
        )
        H0, delta0 = (float(v) for v in result.x)
        chi2 = float(result.fun)
        n = self.data.n
        npar = 2
        return {
            "H0": H0,
            "w": -1.0,
            "omega_de": omega_de,
            "delta0": delta0,
            "chi2": chi2,
            "n": n,
            "npar": npar,
            "dof": n - npar,
            "chi2_dof": chi2 / (n - npar),
            "success": bool(np.isfinite(chi2)),
            "sn_only_stretch": True,
        }

    def fit_h0_w(
        self,
        omega_de: float = OMEGA_DE_TODAY,
        H0_bounds: tuple[float, float] = (50.0, 100.0),
        w_bounds: tuple[float, float] = (-3.0, 0.5),
        H0_guess: float = 70.0,
        w_guess: float = -1.0,
    ) -> dict:
        """Experiment 2: lock present-day Ω_DE, fit constant w and H0."""

        def objective(params) -> float:
            H0, w = float(params[0]), float(params[1])
            return self.chi2(H0, w=w, omega_de=omega_de)

        result = minimize(
            objective,
            x0=np.array([H0_guess, w_guess]),
            method="Nelder-Mead",
            bounds=[H0_bounds, w_bounds],
            options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 500},
        )
        H0, w = (float(v) for v in result.x)
        chi2 = float(result.fun)
        n = self.data.n
        npar = 2
        in_bounds = (
            H0_bounds[0] <= H0 <= H0_bounds[1]
            and w_bounds[0] <= w <= w_bounds[1]
        )
        return {
            "H0": H0,
            "w": w,
            "omega_de": omega_de,
            "chi2": chi2,
            "n": n,
            "npar": npar,
            "dof": n - npar,
            "chi2_dof": chi2 / (n - npar),
            "success": bool(np.isfinite(chi2) and in_bounds),
        }

    def fit_h0_omega_de(
        self,
        w: float = -1.0,
        H0_bounds: tuple[float, float] = (50.0, 100.0),
        omega_bounds: tuple[float, float] = (0.20, 0.90),
        H0_guess: float = 70.0,
        omega_guess: float = 0.70,
    ) -> dict:
        """COMPARISON ONLY.

        Float present-day Ω_DE against the supernovae. This is how we
        measure the χ² *cost of the lock*. It is not a license to replace
        the 1:11 derivation.
        """

        def objective(params) -> float:
            H0, omega_de = float(params[0]), float(params[1])
            return self.chi2(H0, w=w, omega_de=omega_de)

        result = minimize(
            objective,
            x0=np.array([H0_guess, omega_guess]),
            method="Nelder-Mead",
            bounds=[H0_bounds, omega_bounds],
            options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 600},
        )
        H0, omega_de = (float(v) for v in result.x)
        chi2 = float(result.fun)
        n = self.data.n
        npar = 2
        return {
            "H0": H0,
            "w": w,
            "omega_de": omega_de,
            "chi2": chi2,
            "n": n,
            "npar": npar,
            "dof": n - npar,
            "chi2_dof": chi2 / (n - npar),
            "success": bool(np.isfinite(chi2)),
            "comparison_only": True,
        }

    def profile(
        self,
        values: np.ndarray,
        param: str,
        w: float = -1.0,
        omega_de: float = OMEGA_DE_TODAY,
        H0: float | None = None,
        H0_bounds: tuple[float, float] = (50.0, 100.0),
    ) -> np.ndarray:
        """χ² along `param` in {H0, w, omega_de}. Other params as given.

        For w and omega_de, H0 is profiled (optimized) at each grid point.
        """
        chi2 = np.empty(len(values), dtype=float)
        for i, v in enumerate(values):
            if param == "H0":
                chi2[i] = self.chi2(float(v), w=w, omega_de=omega_de)
            elif param == "w":
                chi2[i] = self.fit_h0(w=float(v), omega_de=omega_de, H0_bounds=H0_bounds)["chi2"]
            elif param == "omega_de":
                chi2[i] = self.fit_h0(w=w, omega_de=float(v), H0_bounds=H0_bounds)["chi2"]
            else:
                raise ValueError(param)
        return chi2


def dchi2_interval(x: np.ndarray, chi2: np.ndarray, delta: float = 1.0) -> dict:
    """1-D interval where χ² ≤ min χ² + delta, linearly interpolated."""
    x = np.asarray(x, dtype=float)
    chi2 = np.asarray(chi2, dtype=float)
    imin = int(np.argmin(chi2))
    cmin = float(chi2[imin])
    y = chi2 - cmin
    lo = float(x[0])
    hi = float(x[-1])
    crossed_lo = False
    crossed_hi = False
    for i in range(imin, 0, -1):
        if y[i - 1] >= delta:
            lo = float(np.interp(delta, [y[i], y[i - 1]], [x[i], x[i - 1]]))
            crossed_lo = True
            break
    for i in range(imin, len(x) - 1):
        if y[i + 1] >= delta:
            hi = float(np.interp(delta, [y[i], y[i + 1]], [x[i], x[i + 1]]))
            crossed_hi = True
            break
    best = float(x[imin])
    return {
        "best": best,
        "chi2_min": cmin,
        "lo": lo,
        "hi": hi,
        "minus": best - lo,
        "plus": hi - best,
        "bounded_low": not crossed_lo,
        "bounded_high": not crossed_hi,
    }

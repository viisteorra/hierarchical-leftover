"""Synthetic-catalogue smoke tests: recover known H0, χ² behaves."""

import numpy as np

from cosmology import distance_modulus
from geometry import OMEGA_DE_TODAY
from likelihood import PLANCK_OMEGA_DE, SNData, SNLikelihood


def _synthetic(H0=70.0, w=-1.0, n=40, sigma=0.0, seed=0, zmin=0.02, zmax=1.2):
    rng = np.random.default_rng(seed)
    z = np.linspace(zmin, zmax, n)
    mu = distance_modulus(z, H0, w=w, omega_de=OMEGA_DE_TODAY)
    if sigma > 0.0:
        mu = mu + rng.normal(0.0, sigma, size=n)
        cov = np.diag(np.full(n, sigma**2))
    else:
        # Tiny but not pathological: 1e-12-mag² covariances make χ²
        # gradients explode and the 1-D bounded minimizer looks "failed".
        cov = np.diag(np.full(n, 1e-4**2))
    data = SNData(
        z=z,
        mu=mu,
        cov=cov,
        zmin=zmin,
        n_raw=n,
        path_dat="synthetic.dat",
        path_cov="synthetic.cov",
    )
    return SNLikelihood(data)


def test_chi2_zero_when_model_equals_data():
    like = _synthetic(H0=70.0, w=-1.0, sigma=0.0)
    chi2 = like.chi2(70.0, w=-1.0, omega_de=OMEGA_DE_TODAY)
    assert chi2 < 1e-6


def test_fit_h0_recovers_noiseless_truth():
    true_H0 = 73.2
    like = _synthetic(H0=true_H0, w=-1.0, sigma=0.0)
    fit = like.fit_h0(w=-1.0, omega_de=OMEGA_DE_TODAY)
    assert fit["success"]
    assert abs(fit["H0"] - true_H0) < 0.02
    assert fit["chi2"] < 1e-2
    assert fit["npar"] == 1


def test_fit_h0_w_recovers_noiseless_lambda():
    true_H0 = 71.5
    like = _synthetic(H0=true_H0, w=-1.0, sigma=0.0, n=50)
    fit = like.fit_h0_w(omega_de=OMEGA_DE_TODAY, H0_guess=70.0, w_guess=-0.8)
    assert fit["success"]
    assert abs(fit["H0"] - true_H0) < 0.05
    assert abs(fit["w"] + 1.0) < 0.02
    assert fit["npar"] == 2


def test_noisy_chi2_per_dof_is_order_unity():
    like = _synthetic(H0=70.0, w=-1.0, n=80, sigma=0.05, seed=7)
    fit = like.fit_h0(w=-1.0)
    # 80 points, one parameter; χ²/dof should land near 1, not 0 or 100.
    assert 0.3 < fit["chi2_dof"] < 3.0


def test_wrong_omega_de_increases_chi2_on_noiseless_data():
    like = _synthetic(H0=70.0, w=-1.0, sigma=0.0, n=60)
    chi2_true = like.chi2(70.0, w=-1.0, omega_de=OMEGA_DE_TODAY)
    chi2_wrong = like.chi2(70.0, w=-1.0, omega_de=PLANCK_OMEGA_DE)
    # Shape difference is small (0.69014 vs 0.689) but must be detectable
    # with noiseless data and a long lever arm in redshift.
    assert chi2_true < 1e-6
    assert chi2_wrong > chi2_true


def test_omega_de_is_not_a_free_parameter_of_the_fits():
    like = _synthetic(H0=70.0, sigma=0.0)
    fit1 = like.fit_h0()
    fit2 = like.fit_h0_w()
    assert fit1["omega_de"] == OMEGA_DE_TODAY
    assert fit2["omega_de"] == OMEGA_DE_TODAY
    assert "omega_de" in fit1 and fit1["npar"] == 1
    assert fit2["npar"] == 2  # H0, w — not Ω_DE

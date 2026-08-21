"""Profile-likelihood interval helper and free-Ω_DE comparison flag."""

import numpy as np

from cosmology import distance_modulus
from geometry import OMEGA_DE_TODAY
from likelihood import SNData, SNLikelihood, dchi2_interval


def _parabola_like():
    x = np.linspace(68.0, 76.0, 81)
    chi2 = ((x - 72.0) / 0.5) ** 2 + 10.0
    return x, chi2


def test_dchi2_interval_recovers_one_sigma_of_a_parabola():
    x, chi2 = _parabola_like()
    iv = dchi2_interval(x, chi2, delta=1.0)
    assert abs(iv["best"] - 72.0) < 0.05
    assert abs(iv["minus"] - 0.5) < 0.05
    assert abs(iv["plus"] - 0.5) < 0.05
    assert not iv["bounded_low"]
    assert not iv["bounded_high"]


def test_free_omega_de_is_flagged_comparison_only():
    z = np.linspace(0.05, 1.0, 20)
    mu = distance_modulus(z, 70.0, w=-1.0, omega_de=OMEGA_DE_TODAY)
    like = SNLikelihood(
        SNData(
            z=z,
            mu=mu,
            cov=np.diag(np.full(z.size, 0.02**2)),
            zmin=0.05,
            n_raw=z.size,
            path_dat="s.dat",
            path_cov="s.cov",
        )
    )
    fit = like.fit_h0_omega_de(H0_guess=70.0, omega_guess=0.65)
    assert fit["comparison_only"] is True
    assert abs(fit["omega_de"] - OMEGA_DE_TODAY) < 0.03

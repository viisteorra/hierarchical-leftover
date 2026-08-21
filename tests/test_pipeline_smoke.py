"""End-to-end smoke: both design experiments on a tiny synthetic catalogue."""

import numpy as np

from cosmology import distance_modulus
from geometry import OMEGA_DE_TODAY
from likelihood import SNData, SNLikelihood


def test_both_design_experiments_run_and_keep_omega_de_locked():
    z = np.linspace(0.02, 1.2, 25)
    mu = distance_modulus(z, 72.0, w=-1.0, omega_de=OMEGA_DE_TODAY)
    like = SNLikelihood(
        SNData(
            z=z,
            mu=mu,
            cov=np.diag(np.full(z.size, 1e-4**2)),
            zmin=0.02,
            n_raw=z.size,
            path_dat="synthetic.dat",
            path_cov="synthetic.cov",
        )
    )
    exp1 = like.fit_h0(w=-1.0, omega_de=OMEGA_DE_TODAY)
    exp2 = like.fit_h0_w(omega_de=OMEGA_DE_TODAY, H0_guess=exp1["H0"])

    assert exp1["success"] and exp2["success"]
    assert exp1["omega_de"] == OMEGA_DE_TODAY
    assert exp2["omega_de"] == OMEGA_DE_TODAY
    assert exp1["w"] == -1.0
    assert exp2["chi2"] <= exp1["chi2"] + 1e-6
    assert abs(exp1["H0"] - 72.0) < 0.05

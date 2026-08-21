"""Published flat-ΛCDM Ω_DE measurements. These are 'previous days'.

The 1:11 lock is *not* fitted to this table. The table answers whether
49/71 is a modern Planck-era number or something the data have always wanted.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HistPoint:
    year: float
    name: str
    omega_de: float
    sigma: float
    note: str


# Representative published values, not a complete meta-analysis.
# Ω_DE = 1 − Ω_m for the flat ΛCDM number quoted in each paper.
HISTORY = (
    HistPoint(2003, "WMAP1+other", 0.73, 0.04, "Spergel et al. 2003"),
    HistPoint(2013, "WMAP9", 0.721, 0.025, "Hinshaw et al. 2013, WMAP-only"),
    HistPoint(2014, "Planck13+WP", 0.686, 0.020, "Ade et al. 2014"),
    HistPoint(2016, "Planck15 TT+lensing", 0.692, 0.012, "Ade et al. 2016, Ωm=0.308±0.012"),
    HistPoint(2018, "Planck18 lensing", 0.6847, 0.0073, "Aghanim et al. 2020"),
    HistPoint(2018.5, "Planck18+BAO", 0.6889, 0.0056, "Aghanim et al. 2020 + BAO"),
    HistPoint(2022, "Pantheon+ SN-only", 0.666, 0.018, "Brout et al. 2022, Ωm=0.334±0.018"),
    HistPoint(2024, "DES-SN5YR", 0.648, 0.017, "Vincenzi et al. 2024, Ωm≈0.352"),
    HistPoint(2024.5, "DESI DR1 BAO+BBN", 0.705, 0.015, "Adame et al. 2025, Ωm≈0.295"),
    HistPoint(2025, "DESI DR2 BAO+BBN", 0.703, 0.009, "Abdul-Karim et al. 2025, Ωm≈0.297"),
)

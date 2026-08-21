# Hierarchical dark energy

A 0-parameter geometric lock of the present-day dark-energy density, and a leftover ruler that maps early and local Hubble readings of one flat FRW universe.

**Given primitive P** (the vacuum *is* this hierarchy; its residual tail is \(\Omega_{DE}\)):

\[
\Omega_{DE}=\frac{49}{71},\qquad w=-1,\qquad k=0.
\]

The Hubble *scale* split on this lock is leftover of two readings, not two universes. Leftover \(f\) does not contain \(H_0\).

P is an identification, not a theorem of GR or QFT. If P is false, the numbers are a coincidence.

Spine: [`docs/THEORY.md`](docs/THEORY.md). Hubble sentence: [`docs/CLAIM.md`](docs/CLAIM.md).

## Density

One Euclidean seed on a 12-fold period of square coordinations \(q=4,5\) (1:11). The 12 is the first continued-fraction convergent of \(\log_2(3/2)\) after \(3/5\). Generation multiplier \(r=2/q\). The hierarchical tail is the generating function \(G(\rho)=1/(1-\rho)\):

\[
r=\frac{49}{120},\qquad T=G(r)-1=\frac{49}{71}.
\]

That lock is hard-coded in `src/geometry.py` and is never a fit parameter. SN’s own \(\Omega\) minimum is \(\approx 0.669\), not \(49/71\).

Joint Pantheon+ SH0ES + DESI DR2 BAO + CMB \(R\): lock \(\Delta\chi^2\approx 0.05\).

```python
Q4_WEIGHT = 1
Q5_WEIGHT = 11
r = (Q4_WEIGHT * (2/4) + Q5_WEIGHT * (2/5)) / 12   # 49/120
OMEGA_DE_TODAY = r / (1 - r)                        # 49/71
W = -1.0
K_CURVATURE = 0
```

## Leftover and Hubble

Integer octaves of \(a(t)\) already sit in the metric. Leftover is \(\varphi=\{\log_2(1+z_*)\}\) at photon decoupling (\(z_*=1089.84\)). Same \(G\), two evaluations:

| Set | Exponent | Probe |
|---|---|---|
| Finite (11 generators) | \(e_\square=(12/11)G(r/11)=1440/1271\) | BAO spatial slices |
| Infinite (dimensions) | \(e_\infty=1+(1-r)\ln G(r)\) | CMB \(\theta_*\) |
| This octave | \(f=1\) | SN, time delays |

\[
f=2^{\varphi e},\qquad H_{0,L}=f\,H_{0,E}.
\]

| | Raw | After leftover |
|---|---|---|
| SN vs DESI+\(r_d\) | \(4.79\sigma\) | **\(0.01\sigma\)** |
| SN \(\to\) CAMB \(100\theta_*\) | \(5.01\sigma\) as \(H_0\) | **\(+0.69\sigma\)** |
| BAO pair \(\Delta\chi^2\) | — | **\(+0.00\)** |

Dumping DESI \(H_0=68.38\) into CAMB is the mixed frame (\(+18.6\sigma\) on \(\theta_*\)). Leftover does not eat that. There is no third leftover \(f\).

## Reproduce

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python scripts/download_data.py
python -m pytest tests
python scripts/run_final.py
```

CAMB (optional; Hubble forward prediction and hell H12):

```bash
pip install camb
python scripts/run_hubble.py
python scripts/run_falsify.py
python scripts/run_hell.py
```

## How to kill it

Rejecting P is refusing the premise. The *numbers* die if any of:

- leftover \(f\) depends on \(H_0\)
- 1:11 or 12 is chosen from \(\chi^2\)
- licensed BAO or \(\theta_*\) map \(\ge 3\sigma\)
- licensed BAO \(\Delta\chi^2\ge 10\)
- swapped leftover measures beat the licensed assignment
- CAMB \(100\theta_*\) at \(\mathrm{SN}/f_\infty\) \(\ge 3\sigma\)
- combined SN+BAO+CMB in flat \(\Lambda\) leaves \(49/71\) at \(\gtrsim 3\sigma\) (full Planck likelihood and official DESI covariance count)

## Limits

- CMB shape test here is CAMB \(100\theta_*\), first-peak \(\ell\), and the compressed shift \(R\), not a full Planck \(C_\ell\) likelihood.
- BAO uses published distance ratios, not the official cobaya likelihood.
- Mixed-frame CAMB at DESI \(H_0\) is predicted to fail; that is not a leftover defect.
- Weak lensing is classified (finite set \(\to e_\square\)) and not scored in this repo.

## Data

Pantheon+ SH0ES (Scolnic, Brout, Riess et al. 2022). DES-SN5YR / Dovekie (Vincenzi et al. 2024). DESI DR1/DR2 BAO (Adame et al. 2024; Abdul-Karim et al. 2025). Planck 2018. CMB \(R\): Chen, Huang & Wang 2019.

Catalogues are gitignored. `python scripts/download_data.py` fetches them.

## How this was worked out

This repo was developed in conversation with [Grok](https://x.ai/grok) (xAI): the hierarchy, leftover ruler, tests, and write-up were iterated together. The human side set the axioms and the freeze rules (do not retune 1:11 from \(\chi^2\), do not put \(H_0\) into leftover \(f\), do not lock \(T/10\)). The model side implemented, tried to kill the numbers, and wrote the proofs under those constraints.

That is disclosed for transparency. It is not a credential. The claim is the locked formulae and the tests, not who typed them.

## License

MIT.

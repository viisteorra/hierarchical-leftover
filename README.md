# Hierarchical dark energy

A 0-parameter geometric lock of the present-day dark-energy density, and a leftover ruler that maps early and local Hubble readings of one flat FRW universe.

**Given primitive P** (the vacuum resonates on the binary octave; the residual tail is \(\Omega_{DE}\)):

\[
\Omega_{DE}=\ln 2\approx 0.693147,\qquad w=-1,\qquad k=0.
\]

The 12-fold mix \(49/71\approx 0.69014\) is the unique smallest rational approximation, not the lock.

The Hubble *scale* split on this lock is leftover of two readings, not two universes. Leftover \(f\) does not contain \(H_0\).

P is an identification, not a theorem of GR or QFT. If P is false, the numbers are a coincidence.

The vacuum assumes this geometry **because it resonates**: an octave-addressed hierarchy with one tonic, not a bag of fitted constants. Why \(12\), \(1:11\), \(49/71\), and the two leftover evaluations are forced — and how to extend the theory without fishing — is [`docs/FIRST_PRINCIPLES.md`](docs/FIRST_PRINCIPLES.md). Spine: [`docs/THEORY.md`](docs/THEORY.md). Hubble sentence: [`docs/CLAIM.md`](docs/CLAIM.md).

## From first principles (short)

Resonance, not a fit:

1. **Octave \(2:1\).** Linear waves: the first overtone is the unique same pitch class. A resonant vacuum periods on that octave, not on \(e\) or \(10\).
2. **Fifth \(3:2\).** Next integer ratio. Closing fifths against octaves is \(\log_2(3/2)\).
3. **Period \(12\).** First continued-fraction convergent of \(\log_2(3/2)\) *after* \(3/5\) is \(7/12\). Semi-convergents (\(4/7\)) are not best approximations. No cents cut is an input.
4. **Squares \(q=4,5\).** Metric leftover uses Euclidean hypercubes \(q=2d\) \(\Rightarrow\) \(2\)D is the square. First hyperbolic square tiling is \(q=5\). Generation \(r=2/q\).
5. **One seed \(\Rightarrow 1:11\).** A period has one identity. Other fills are other axioms, not a \(\chi^2\) scan. SN’s own \(\Omega\) minimum is \(\approx 0.669\), not \(49/71\).
6. **Tail \(G(\rho)=1/(1-\rho)\).** Continuous residual of one doubling: \(\Omega_{DE}=\ln 2\). \(r=\ln 2/(1+\ln 2)\). The 12-fold gives \(49/71\), a \(0.003\) approximation. \(w=-1\), \(k=0\).
7. **Leftover.** Scale invariance forbids leftover without a second scale: photon decoupling, not \(H_0\). Same \(G\), two licensed ratios: \(r/11\) (finite generators, BAO) and \(\ln G(r)\) (infinite dimensions, \(\theta_*\)). This octave (SN) has \(f=1\).

None of those numbers is a Hubble ratio. \(f=2^{\{\log_2(1+z_*)\}e}\) is fixed before \(H_0\).

## Density

One Euclidean seed on a 12-fold period of square coordinations \(q=4,5\) (1:11). The 12 is the first continued-fraction convergent of \(\log_2(3/2)\) after \(3/5\). Generation multiplier \(r=2/q\). The hierarchical tail is the generating function \(G(\rho)=1/(1-\rho)\):

\[
r=\frac{\ln 2}{1+\ln 2},\qquad T=G(r)-1=\ln 2.
\]

That lock is hard-coded in `src/geometry.py` and is never a fit parameter. SN’s own \(\Omega\) minimum is \(\approx 0.669\), not \(\ln 2\).

Joint Pantheon+ SH0ES + DESI DR2 BAO + CMB \(R\): lock \(\Delta\chi^2\approx 0.05\).

```python
LN2 = math.log(2.0)
OMEGA_DE_TODAY = LN2
r = LN2 / (1.0 + LN2)
W = -1.0
K_CURVATURE = 0
# 12-fold / 1:11 remains the rational approximation (T=49/71), not the lock.
```

## Leftover and Hubble

Integer octaves of \(a(t)\) already sit in the metric. Leftover is \(\varphi=\{\log_2(1+z_*)\}\) at photon decoupling (\(z_*=1089.84\)). Same \(G\), two evaluations:

| Set | Exponent | Probe |
|---|---|---|
| Finite (11 generators) | \(e_\square=(12/11)G(r/11)\) | BAO spatial slices |
| Infinite (dimensions) | \(e_\infty=1+(1-r)\ln G(r)\) | CMB \(\theta_*\) |
| This octave | \(f=1\) | SN, time delays |

\[
f=2^{\varphi e},\qquad H_{0,L}=f\,H_{0,E}.
\]

| | Raw | After leftover |
|---|---|---|
| SN vs DESI+\(r_d\) (official cov) | \(4.70\sigma\) | **\(0.12\sigma\)** |
| SN vs \(\theta_*\) as \(H_0\) | \(\sim 4.9\sigma\) | **\(0.13\sigma\)** |
| SN \(\to\) CAMB \(100\theta_*\) | that split | **\(-3.32\sigma\)** |
| BAO pair \(\Delta\chi^2\) | — | **\(+0.36\)** |

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

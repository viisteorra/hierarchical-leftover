# The theory (one generating function)

Primitive **P**: the vacuum *resonates* on the binary octave \(2:1\). In continuous log-space the measure of one doubling is \(\ln 2\); that residual tail is \(\Omega_{DE}\). The rationals \(49/71\), \(49/120\), \(71/120\) were the 12-fold images of \(\ln 2\), \(r\), and \(1-r\). Leftover addressing is still 12-fold (\(s=12\varphi\), \(12/11\)). Everything else is lemmas + one measurement \(z_*\).

Why the values are not arbitrary, and how to take it further: `FIRST_PRINCIPLES.md`. Proofs: `THEOREM.md`. Hubble sentence: `CLAIM.md`.

## Hierarchy

Octave \(2:1\) (waves). Fifth \(3:2\) (next integer). Period \(12\) = first continued-fraction convergent of \(\log_2(3/2)\) after \(3/5\). Squares \(q=4=2d\), first hyperbolic \(q=5\). One Euclidean seed \(\Rightarrow 1:11\).

\[
\Omega_{DE}=\ln 2,\qquad r=\frac{\ln 2}{1+\ln 2},\qquad G(\rho)=\frac{1}{1-\rho}=\sum_{k=0}^{\infty}\rho^k.
\]

\(T=G(r)-1=\ln 2\). Then \(r=\ln 2/(1+\ln 2)\) and \(1-r=1/(1+\ln 2)\) (old \(49/120\) and \(71/120\)). FRW \(k=0\).

## One \(G\), three evaluations, two leftovers

\[
T=G(r)-1=\ln 2=\Omega_{DE},\qquad \Omega_m=1-\ln 2,\qquad w=-1.
\]

Leftover of two readings of a frozen ruler (photon decoupling \(z_*\)). Density is continuous; **the reading is 12-fold:**

\[
\log_2(1+z_*)=n+\varphi,\qquad s=12\varphi,\qquad f=2^{(s/12)e}.
\]

| Generating set | Exponent | Probe |
|---|---|---|
| Finite (11 hyperbolic steps of the 12-fold) | \(e_\square=(12/11)\,G(r/11)\) | spatial slices: BAO, weak lensing |
| Infinite (dimensions, continuous \(r\)) | \(e_\infty=1+(1-r)\ln G(r)\) | null metric: \(\theta_*\) |
| This octave | \(e=0\Rightarrow f=1\) | SN, time delays |

No third leftover. \(G(r/12)\) is the clock-unit foil (Lemma 3: leftover lives on the 11). Do not mix the two \(f\) on one pair.

\[
H_{0,L}=f\,H_{0,E}.
\]

\(H_0\) is not in \(G\). \(\chi^2\) is not in \(G\).

## Closed

Given P, the theorem is complete. There is no further \(f\), no further fill, no further period. New probes reuse the table above: classify the generating set, do not minimise \(\chi^2\). Analogous theory: swap one uniqueness lemma (e.g. \(53\)-fold) and copy \(G\).

Not remaining work (and not this theorem): proving P from GR/QFT; official DESI full covariance / Planck \(C_\ell\).

Code: `src/generate.py`.

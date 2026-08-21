# The theory (one generating function)

Primitive **P**: the vacuum *resonates*; it is this octave-addressed hierarchy, and the residual tail of that resonance is \(\Omega_{DE}\). The geometry is assumed *because of that resonance*, not because a catalogue preferred these integers. Everything below is lemmas + one measurement \(z_*\).

Why the values are not arbitrary, and how to take it further: `FIRST_PRINCIPLES.md`. Proofs: `THEOREM.md`. Hubble sentence: `CLAIM.md`.

## Hierarchy

Octave \(2:1\) (waves). Fifth \(3:2\) (next integer). Period \(12\) = first continued-fraction convergent of \(\log_2(3/2)\) after \(3/5\). Squares \(q=4=2d\), first hyperbolic \(q=5\). One Euclidean seed \(\Rightarrow 1:11\).

\[
r=\frac{49}{120},\qquad G(\rho)=\frac{1}{1-\rho}=\sum_{k=0}^{\infty}\rho^k.
\]

\(G\) is Theorem 1’s series, not a second object. FRW \(k=0\): metric leftover is Euclidean hypercubes; hyperbolic \(q=5\) lives in the 2D vacuum mix, not the 3-geometry.

## One \(G\), three evaluations, two leftovers

\[
T=G(r)-1=\frac{49}{71}=\Omega_{DE},\qquad \Omega_m=\frac{22}{71},\qquad w=-1.
\]

Leftover of two readings of a frozen ruler (photon decoupling \(z_*\)): \(\varphi=\{\log_2(1+z_*)\}\), \(f=2^{\varphi e}\).

| Generating set | Exponent | Probe |
|---|---|---|
| Finite (11 hyperbolic steps) | \(e_\square=(12/11)\,G(r/11)=1440/1271\) | spatial slices: BAO, weak lensing |
| Infinite (dimensions) | \(e_\infty=1+(1-r)\ln G(r)\) | null metric: \(\theta_*\) |
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

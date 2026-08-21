# Solution (given P)

Primitive **P**: the vacuum is this hierarchy; its tail is \(\Omega_{DE}\).
Uniqueness U1–U3 and leftover proofs: `AXIOMS.md`, `THEOREM.md`. No \(H_0\) in the leftover.

## Density

\[
\Omega_{DE}=\frac{49}{71},\qquad w=-1,\qquad \text{flat}.
\]

SN + DESI DR2 + CMB \(R\): lock \(\Delta\chi^2\approx 0.05\).

## Leftover (same \(\varphi\), two measures)

\(\varphi=\{\log_2(1+z_*)\}\), \(z_*=1089.84\) (photon decoupling / CAMB lock).

\[
f_\square=2^{\varphi\cdot 1440/1271}\approx 1.074268
\quad\text{(BAO, planar + generator tail)}
\]
\[
f_\infty=2^{\varphi[1-(1-r)\ln(1-r)]}\approx 1.086396
\quad\text{(CMB }\theta_*\text{, infinite Euclidean)}
\]

## Hubble scale (Theorem 3)

Raw lock split: SN \(73.47\) vs DESI+\(r_d\) \(68.38\) = \(5.09\) km/s (**4.79σ**); vs \(\theta_*\)-lock \(67.60\) = \(5.87\) km/s (**5.01σ**).

| Direction | Input | Map | Output | vs data |
|---|---|---|---|---|
| inverse | DESI+\(r_d\) \(68.38\) | \(\times f_\square\) | **73.46** | SN \(73.47\), **0.01σ** |
| inverse | CAMB \(\theta_*\) \(67.60\) | \(\times f_\infty\) | **73.44** | SN \(73.47\), **0.03σ** |
| forward | SN \(73.47\) | \(/f_\square\) | **68.39** | DESI+\(r_d\) \(68.38\), **0.04σ** |
| forward | SN \(73.47\) | \(/f_\infty\) | **67.63** | CAMB \(100\theta_*=1.04132\), **+0.69σ** |

Corollary: \(f_\infty/f_\square\approx 1.01129\) vs DESI/\(\theta_*\approx 1.01160\) (**0.04σ**). K2 \(\Delta\chi^2=+0.00\). Hybrid unused. First peak \(\ell=220\).

Hell 0/46. Falsify 0/10. Nails 0/20. CMB-frame \(\omega_m\) scored in the \(\theta_*\) frame. Spine: `THEORY.md` (\(G(\rho)=1/(1-\rho)\)).

## What is solved

- Dark-energy *density* as a locked fraction, given P.
- Hubble *scale* split: SH0ES vs inverse-ladder and SH0ES vs Planck \(\theta_*\), given P and two leftover measures.
- DESI vs \(\theta_*\) *Hubble readings* as the leftover-measure ratio (not a third \(f\)).
- Raw 4.8σ on this lock → **0.01σ** (BAO) / **0.69σ** (forward CAMB \(\theta_*\)). Inverse \(\theta_*\) map in SH0ES \(H_0\) units is **0.03σ**.

## What is not solved (and is not supposed to be)

- P itself (not in GR/QFT).
- One leftover \(f\) for every probe.
- Dumping DESI \(H_0=68.38\) into CAMB (\(+18.6\sigma\) on \(\theta_*\)): mixed frame.
- Official DESI full covariance / full Planck \(C_\ell\) likelihood.

The sentence: `CLAIM.md`. Spine: `THEORY.md`. U1 is the fifth’s continued fraction (not a 2-cent cut). Lemma C is finite vs infinite generating set of the same \(G\). \(k=0\) is Euclidean metric leftover, not a fit of \(\Omega_k\).

If P is true, this is the theorem. If P is false, the numbers are a coincidence. There is no further \(f\) to tune.

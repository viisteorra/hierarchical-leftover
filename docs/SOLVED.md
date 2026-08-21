# Solution (given P)

Primitive **P**: the vacuum resonates on the binary octave; the residual tail is \(\Omega_{DE}\).
Uniqueness U1–U3: `AXIOMS.md`. Continuous lock and leftover: `THEOREM.md`, `THEORY.md`, `FIRST_PRINCIPLES.md`. No \(H_0\) in leftover.

## Density

\[
\Omega_{DE}=\ln 2\approx 0.693147,\qquad w=-1,\qquad k=0,\qquad \Omega_m=1-\ln 2.
\]

\(r=\ln 2/(1+\ln 2)\), \(1-r=1/(1+\ln 2)\) (old \(49/120\), \(71/120\)). Leftover addressing is still 12-fold. Not a \(\chi^2\) retune.

Official DESI DR2 13×13 covariance. SN lock \(\Delta\chi^2=1.77\) vs free \(\Omega\). BAO \(\chi^2=11.42\) on 12 dof.

## Leftover (same \(\varphi\), two measures)

\(\varphi=\{\log_2(1+z_*)\}\), \(z_*=1089.84\), \(r=\ln 2/(1+\ln 2)\).

\[
f_\square=2^{\varphi\cdot(12/11)G(r/11)}\approx 1.074275
\quad\text{(BAO; }1440/1271\text{ was the }49/120\text{ image)}
\quad\text{(BAO, finite generating set)}
\]
\[
f_\infty=2^{\varphi[1-(1-r)\ln(1-r)]}\approx 1.086430
\quad\text{(CMB }\theta_*\text{, infinite Euclidean)}
\]

## Hubble scale (Theorem 3)

Raw lock split: SN \(73.50\) vs DESI+\(r_d\) \(68.54\) = \(4.97\) km/s (**4.70σ**).

| Direction | Input | Map | Output | vs data |
|---|---|---|---|---|
| inverse | DESI+\(r_d\) \(68.54\) | \(\times f_\square\) | **73.63** | SN \(73.50\), **0.12σ** |
| inverse | CAMB \(\theta_*\) \(67.80\) | \(\times f_\infty\) | **73.66** | SN \(73.50\), **0.13σ** |
| forward | SN \(73.50\) | \(/f_\square\) | **68.42** | DESI+\(r_d\) \(68.54\) |
| forward | SN \(73.50\) | \(/f_\infty\) | **67.66** | CAMB \(100\theta_*=1.04007\), **−3.32σ** |

K2 \(\Delta\chi^2=+0.36\). Hybrid unused. First peak \(\ell=220\). Inverse maps close the SH0ES-scale split. Forward CAMB \(\theta_*\) at leftover-predicted \(H_0\) is **−3.32σ** because \(\Omega_m=1-\ln 2\); leftover \(f\) does not contain \(H_0\) and does not eat that shape.

## What is solved

- Dark-energy *density* as \(\ln 2\), given P. 12-fold is the approximation.
- Hubble *scale* inverse maps: 4.7σ → 0.12σ (BAO) and 0.13σ (\(\theta_*\) as \(H_0\)).

## What is not solved

- P itself.
- Forward CAMB \(100\theta_*\) at \(\mathrm{SN}/f_\infty\) (−3.32σ).
- One leftover \(f\) for every probe.
- Mixed-frame CAMB at DESI \(H_0\).

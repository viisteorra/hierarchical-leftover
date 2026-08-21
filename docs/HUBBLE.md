# Hubble tension, given P

Leftover \(f\) is Theorem 2; it does not contain \(H_0\). Accounting: `src/hubble.py`. Live CAMB: `scripts/run_hubble.py`.

## What the tension was

On this lock (\(\Omega_{DE}=49/71\), \(w=-1\)):

| | \(H_0\) | vs SN \(73.47\pm 1.04\) |
|---|---|---|
| Pantheon+ SH0ES (lock \(\Omega\)) | \(73.47\) | — |
| DESI DR2 + Planck \(r_d\) | \(68.38\) | **4.79σ** (\(5.09\) km/s) |
| CAMB \(\theta_*\) lock | \(67.60\) | **5.01σ** |
| Planck 2018 catalog | \(67.36\) | **5.22σ** |

That is the ordinary Hubble scale split. Stretch, \(w(a)\), stacks, and one leftover \(f\) on every probe do not eat it. Leftover preserves \(H_0 r_d\), so mixed-frame CAMB at DESI \(H_0\) stays \(+18.6\sigma\) on \(\theta_*\).

## What leftover does

Same \(\varphi=\{\log_2(1+z_*)\}\). BAO reads planar \(f_\square=2^{\varphi\cdot 1440/1271}\). \(\theta_*\) reads infinite Euclidean \(f_\infty=2^{\varphi[1-(1-r)\ln(1-r)]}\).

**Inverse (early → local).** DESI \(68.38\times f_\square=73.46\) (**0.01σ**). \(\theta_*\)-lock \(67.60\times f_\infty=73.44\) (**0.03σ** in SH0ES \(H_0\) units).

**Forward (local → early).** SN \(/f_\square=68.39\) vs DESI \(68.38\) (**0.04σ**). SN \(/f_\infty=67.63\); CAMB at that \(H_0\) gives \(100\theta_*=1.04132\) vs Planck \(1.04110\pm 0.00031\) (**+0.69σ**), \(r_d=147.28\) (**+0.71σ**), \(z_*=1089.85\) (**+0.24σ**), first peak \(\ell=220\).

**Corollary.** \(H_{0,\mathrm{BAO}}/H_{0,\theta_*}=f_\infty/f_\square\). Predicted \(1.01129\) vs measured \(1.01160\) (**0.04σ**). The DESI vs \(\theta_*\) *scale* split is the two-measure ratio. It is not a third leftover.

## Residuals that are not Hubble scale

- Mixed-frame CAMB at DESI \(H_0=68.38\): \(+18.6\sigma\) on \(\theta_*\). Wrong pair, wrong measure.
- BAO-frame \(\omega_m\) at \(68.38\): \(+1.73\sigma\). CMB-frame \(\omega_m\) at SN\(/f_\infty\): \(\sim-1.2\sigma\).
- Full Planck \(C_\ell\) likelihood and official DESI covariance: not run here.
- P is unproved.

## Claim

Given P and two leftover measures, the Hubble *scale* split on this lock is closed: **4.8σ → 0.01σ (BAO) and 0.69σ (Planck \(\theta_*\) predicted from SN)**. That is not a fit of \(f\) to \(H_0\), and it is not a proof of P.

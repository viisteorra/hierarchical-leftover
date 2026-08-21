# Analogous axioms and uniqueness lemmas

No \(H_0\). No \(\chi^2\). Code: `src/uniqueness.py`.

You do not prove axioms **in the same system**. You prove them *elsewhere* (waves, tilings, continued fractions, operational meaning of rulers) until **one** physical identification remains. That remainder is a known primitive, not an unknown unknown.

## The one primitive (cannot be proved)

**P.** The vacuum *is* this hierarchy, and its residual tail is the present dark-energy density.

Maths cannot extract P from GR or from QFT (the cosmological-constant problem is a different number). P is the dictionary between a combinatorial tail and \(T_{\mu\nu}\). Everything below is from P plus *elsewhere*.

## Proved elsewhere (no longer axioms)

**Octave \(2:1\).** Elsewhere: linear waves. The first overtone \(2f\) is the unique interval that is the same pitch class. If the vacuum is harmonic (P), its period is that octave. Not \(e\), not \(10\).

**U1 — period 12.** Elsewhere: continued fraction of \(\log_2(3/2)\). Convergents: \(1/1\), \(1/2\), \(3/5\), \(7/12\), \(24/41\), \(31/53,\ldots\). The period is the denominator of the first convergent *after* \(3/5\), uniquely \(7/12\). Semi-convergents (\(4/7\)) are not best approximations and are not used. No cents threshold is an input. (Check: every naive cut in \((1.96, 16.2)\) cents also yields \(n=12\); \(2\) sits in that basin, it is not fitted. \(16.2\) is \(4/7\), which CF already discarded.)

**U2 — \(q=4,5\).** Elsewhere: regular tilings \(\{4,q\}\). Euclidean \(\Leftrightarrow q=4\); first hyperbolic \(q=5\). Unique boundary pair for squares. (Squares: unique translationally invariant 2D lattice \(\mathbb{Z}^2\) has degree 4.)

**U3 — one seed \(\Rightarrow 1:11\).** Elsewhere: a period has one identity (the seed / tonic). One identity in 12 parts \(\Rightarrow n_4=1\), \(n_5=11\). Unique.

**\(w=-1\).** Lemma from P + U1: if the period is a constant of the hierarchy, it does not drift, so the residual density does not dilute. That is \(w=-1\), not a fit.

**B — last scattering.** Scale invariance of the octave forbids leftover without *some* second scale. Elsewhere: leftover is the mismatch of two readings of the **same standard ruler**. BAO and \(\theta_*\) use the sound horizon, which is frozen at photon/baryon decoupling. So the second scale is that freeze-out, not \(H_0\), not \(z_{\mathrm{eq}}\). (\(z_*\) vs \(z_d\): leftover uses **photon** decoupling, when the acoustic *pattern* is last scattered. \(r_d\) is the drag horizon; using \(z_d\) as leftover’s second scale is a foil, scored in `scripts/run_close_gaps.py`, not a retune of B.)

**D — one metric.** Elsewhere: “two readings” *means* one line element, two frames. A second metric would be a different theory.

**C — which \(f\) which probe (lemma).** Given P+U2: density lives in the 2D vacuum mix; leftover of that mix is planar \(f_\square\). The tail that produces \(T=r/(1-r)\) is already an infinite generating function; leftover of that tail on Euclidean hypercubes \(q=2d\) is \(f_\infty\). BAO \(D_M,D_H,D_V\) are spatial rulers on the vacuum-filled slices \(\Rightarrow f_\square\). \(\theta_*\) is a null angle through the metric \(\Rightarrow f_\infty\). One pair, one measure: mixed \(r_d/f_\infty\) on BAO is two measures on one pair. Not chosen from \(\chi^2\). Code: `spacetime.measure_for`.

Under P + U1–U3:

\[
r=\frac{49}{120},\qquad T=\frac{49}{71}=\Omega_{DE}.
\]

That is Theorem 1. Leftover formulas are Theorem 2. No \(\chi^2\) in either.

## What is *not* an unknown unknown

An unknown unknown would be an extra parameter we did not know we were using. Leftover \(f\) is fixed before \(H_0\); 1:11 is not a fit; hell 0/41. The remaining risk is **P is false** — a known, single identification. If P is false, \(49/71\) is a coincidence. If P is true, the maths is the theorem.

## Uniqueness lemmas (math)

### U1 — period is 12

**Lemma.** The period is the denominator of the first continued-fraction convergent of \(\log_2(3/2)\) after \(3/5\). That convergent is uniquely \(7/12\).

**Proof.** \(\log_2(3/2)\) has continued-fraction convergents \(1/1\), \(1/2\), \(3/5\), \(7/12\), \(24/41\), \(31/53,\ldots\). Convergents are the unique best rational approximations. Semi-convergents such as \(4/7\) are discarded. After \(3/5\), the next convergent is \(7/12\). No cents cut enters.

**Check (not the proof).** Errors of \(2^{k/n}\) vs \(3:2\): \(7/12\) \(-1.96\) cents, \(4/7\) \(-16.2\) cents (semi-convergent, discarded), \(3/5\) \(+18.0\) cents. Every naive cut in \((1.96, 16.2)\) selects \(n=12\). The figure \(2\) cents sits in the basin; it is not an input.

The fourth \(4:3\) is \(5/12\) at \(+1.96\) cents (octave complement of the fifth). That is the same lemma, not a second period.

### U2 — square coordinations \(q=4\) and \(q=5\)

**Lemma.** The unique Euclidean square tiling of the plane is \(\{4,4\}\) (\(q=4\) squares per vertex). The first hyperbolic square tiling is \(\{4,5\}\) (\(q=5\)).

**Proof.** Regular tilings \(\{p,q\}\): \((p-2)(q-2)=4\) Euclidean, \(<4\) spherical, \(>4\) hyperbolic. For squares \(p=4\): \(2(q-2)=4\Rightarrow q=4\) Euclidean; \(q=5\) is the first integer with \(2(q-2)>4\). No other integer \(q\) sits on that Euclidean/hyperbolic boundary.

Generation multiplier \(r=2/q\) (octave ratio over coordination) is the definition of the hierarchy, not a fit.

### U3 — one Euclidean seed ⇒ 1:11

**Lemma.** A 12-period that is mixed (both coordinations) and has **exactly one** Euclidean seed has weights \(1:11\).

**Proof.** \(n_4+n_5=12\), \(n_4,n_5\ge 1\). “Exactly one seed” is \(n_4=1\), hence \(n_5=11\). Unique. Other mixed fills \(2:10,\ldots,11:1\) exist; they are other axioms, not this one. 1:11 is not chosen from \(\chi^2\).

Then
\[
r=\frac{1\cdot(2/4)+11\cdot(2/5)}{12}=\frac{49}{120},\qquad
T=\frac{r}{1-r}=\frac{49}{71}.
\]

That is Theorem 1 under P + U1–U3.

## Primitive vs lemma (summary)

| Claim | Status |
|---|---|
| Vacuum = this hierarchy, tail = \(\Omega_{DE}\) | **P** — only primitive |
| Octave \(2:1\) | elsewhere (waves) given P |
| 12, \(q=4,5\), \(1:11\), \(T=49/71\) | U1–U3 + algebra |
| \(w=-1\) | lemma: period does not drift |
| Second scale = last scattering | elsewhere: sound-horizon freeze |
| One metric | definition of two readings |
| \(G(\rho)=1/(1-\rho)\) | Theorem 1’s series \(\sum\rho^k\), not a new knob |
| \(k=0\) | lemma: metric leftover is Euclidean hypercubes |
| BAO \(\to f_\square\), \(\theta_*\to f_\infty\) | lemma C: finite vs infinite generating set of \(G\) |

Analogous theory: keep P, swap a uniqueness lemma (e.g. 53-fold). Tail and leftover proofs copy.

## Empirical (not in the proof)

`scripts/run_hell.py`, `scripts/run_falsify.py`, `scripts/run_hubble.py`, `scripts/run_nails.py`. Density sits (\(\Delta\chi^2\approx 0.05\)). SN-free \(\Omega\approx 0.669\), not \(49/71\). Leftover maps: BAO \(0.01\sigma\), inverse \(\theta_*\) \(0.03\sigma\), forward CAMB \(100\theta_*\) \(+0.69\sigma\). Early–early ratio \(f_\infty/f_\square\) \(0.04\sigma\). Swapped measures fail (inf-on-BAO \(\Delta\chi^2=+18\)). Hybrid forbidden. Dumping BAO \(H_0=68.38\) into CAMB is mixed-frame (\(+18.6\sigma\) on \(\theta_*\)); leftover does not eat it. Hell **0/46**. Nails **0/20**.

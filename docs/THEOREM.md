# The leftover theorem

No \(H_0\) in the proof. No \(\chi^2\) in the proof.

Analogous axioms **A–D** and uniqueness **U1–U3** (12 from the fifth, \(\{4,5\}\), one seed \(\Rightarrow 1:11\)): `AXIOMS.md`, `src/uniqueness.py`.

## Axioms

**A (vacuum).** Octave-addressed hierarchy; residual tail \(T\equiv\Omega_{DE}\) today, flat, \(w=-1\).

**B (second scale).** Photon decoupling. Leftover is not \(H_0\).

**C (measures, lemma).** BAO reads planar leftover \(f_\square\) (spatial slices of the 2D vacuum). \(\theta_*\) reads infinite-D leftover \(f_\infty\) (null angle in the metric). One pair, one measure. Not a fit.

**D (one metric).** \(ds^2=-c^2\,dt^2+a(t)^2\,d\chi^2\). Early and local are two readings.

**U1–U3 (math).** Period \(=12\) (first CF convergent of \(\log_2(3/2)\) after \(3/5\)). Squares \(q=4=2d\) (2D hypercube), \(q=5\) first hyperbolic. Exactly one Euclidean seed \(\Rightarrow 1:11\).

Under A and U1–U3 the old G1–G3 are lemmas: period 12, fill 1:11, \(r=2/q\), \(T=r/(1-r)\). G4 is D. G5 is B. C is the two-measure rule.

## Theorem 1 (density)

\[
\Omega_{DE}=T=\frac{49}{71}.
\]

**Proof.** G2:

\[
r=\frac{1\cdot(2/4)+11\cdot(2/5)}{12}=\frac{4.9}{12}=\frac{49}{120}.
\]

G3:

\[
T=\frac{r}{1-r}=\frac{49/120}{71/120}=\frac{49}{71}.
\]

Flatness closes \(\Omega_m=22/71\). No likelihood is used. \(49=7^2\) is not a primality claim and is not used.

## One generating function (Theorems 1–2 are evaluations of \(G\))

Let \(G(\rho)=1/(1-\rho)\). Density and both leftovers are the same tail:

\[
T=G(r)-1,
\qquad
e_\square=\frac{12}{11}\,G\!\left(\frac{r}{11}\right),
\qquad
e_\infty=1+(1-r)\ln G(r).
\]

**Identities.** \(G(r)-1=r/(1-r)\). Finite: \((12/11)/(1-r/11)=(12/11)G(r/11)=1440/1271\). Infinite: \(\ln G(r)=-\ln(1-r)\), hence \(1+(1-r)\ln G(r)=1-(1-r)\ln(1-r)\).

**No third leftover.** U1–U3 put two ratios into \(G\): \(r\) (density) and \(r/11\) (leftover on the 11 generators). \(G(r/12)\) is the clock-unit foil (Lemma 3). \(T/N\) is not an evaluation of \(G\). Spine: `THEORY.md`. Code: `src/generate.py`.

## Lemma (no integer \(N\))

There is no theorem \(f=1/(1-T/N)\) for an integer \(N\) forced by G1–G4.

**Proof.** The internal integers of G1–G2 are \(\{1,2,4,5,11,12\}\). \(10\notin\) that set. Arithmetic coincidences \(11-1=10\) and \(12-2=10\) are not generative rules of the fill. Even among internal integers, \(N=12\), \(N=11\), and \(N=1\) are equally licensed as partitions of \(T\); G1–G4 do not pick one. The geometric series that produces \(T=r/(1-r)\) does not select a partition of \(T\) at all.

Separately: G1 is homogeneous in logarithmic address. Every octave is equivalent. A scale-invariant hierarchy has no intrinsic number of octaves between “then” and “now.” An integer \(\lfloor\log_2(1+z_*)\rfloor\) requires G5, and even with G5 that integer is the *completed* count, not a license to divide \(T\).

So \(N=10\) and \(f=1/(1-T/10)=710/661\) remain an **ansatz**. They are not this theorem.

## Lemma 2 (clock leftover)

Under G1, G4, and G5 only — omitting G2 — the leftover of the *metric expansion* is unique:

\[
f_{\mathrm{clock}}=2^{\{\log_2(1+z_*)\}}=\frac{1+z_*}{2^{\lfloor\log_2(1+z_*)\rfloor}}.
\]

**Proof.** Let \(a_0=1\) today and \(a_*=(1+z_*)^{-1}\) at photon decoupling (G5). G1 addresses those scales by \(\nu(a)=\log_2(a^{-1})\). Split \(\log_2(1+z_*)=n+\varphi\) with \(\varphi=\{\log_2(1+z_*)\}\in[0,1)\). The integer \(n\) is \(n\) completed octaves of expansion: those doublings already sit in \(a(t)\) (G4) and cannot be a mismatch between two readings. The remainder as a scale-factor ratio (G1, base 2) is \(2^{\varphi}=(1+z_*)/2^n\).

This lemma does not use G2. It is the leftover of equal-log FRW, not of the mixed continuum.

## Lemma 3 (direction of the mix)

G1+G5 produce \(\varphi\) already in octave units. G1 identifies one octave with one 12-part period, so \(\varphi\) is a fraction of the **clock**. G2's leftover-generating set is the 11 hyperbolic steps; the Euclidean seed is the identity of the period (the reference reading), not a mismatch.

A fraction \(\varphi\) of a 12-part clock, referred to an 11-part generating set, is \(\varphi\cdot 12/11\) of that set.

The inverse \(11/12\) would take a fraction of the generating set and refer it to the clock. We do not have leftover in generator units from G1. Direction is forced: **\(12/11\), not \(11/12\)**.

G1 is logarithmic, so a change of measure rescales the exponent, not the scale factor: \(2^{\varphi\cdot 12/11}\), not \(2^{\varphi}\times 12/11\).

## Lemma 4 (planar leftover)

Under G1, G2, G4, and G5, leftover referred to the 2D vacuum mix is

\[
f_2=2^{\varphi\cdot 12/11}.
\]

Proof: Lemma 2 plus Lemma 3. This omits the ambient 3-space and 4-spacetime of the metric.

## Lemma 5 (Euclidean hypercubes)

In \(d\) Euclidean dimensions a hypercube has coordination \(q=2d\). Consecutive ratios:

\[
\frac{q_d}{q_{d-1}}=\frac{d}{d-1}.
\]

So \(d=3\) (space) gives \(6/4=3/2\), and \(d=4\) (spacetime) gives \(8/6=4/3\). Same direction as Lemma 3: refer the lower measure to the higher. Inverses \(2/3\) and \(3/4\) are forbidden. \(d=2\) Euclidean-only would be \(2/1=2\), but the 2D *vacuum* is the mixed fill G2, not a pure square lattice, so \(d=2\) uses \(12/11\) not \(2\).

FRW has no fifth Euclidean spatial dimension. Hyperbolic \(q=5\) is already inside the 2D mix.

## Theorem 2 (leftover)

Spacetime is infinite-dimensional. Euclidean hypercubes have \(q=2d\), leftover exponents \(e_1=1\) and \(e_d=d/(d-1)\) for \(d\ge 2\). A **flat** mean over \(d=1,\ldots,N\) goes to \(1\) as \(N\to\infty\) (the clock). The unique regulator already in the theory is the hierarchical tail, weights \(r^{d-1}\):

\[
\langle e\rangle=\frac{\sum_{d=1}^{\infty} r^{d-1} e_d}{\sum_{d=1}^{\infty} r^{d-1}}=1-(1-r)\ln(1-r),
\qquad
f=2^{\{\log_2(1+z_*)\}\cdot\langle e\rangle}.
\]

\(r=49/120\) is the mix generator (A1), not a partition of \(T\). Integer octaves already live in \(a(t)\). \(H_0\) does not enter. There is no cutoff \(N\).

**Proof.** Lemma 2: clock leftover \(\varphi\). Lemma 5: \(e_d=d/(d-1)\). G1 is logarithmic, so the mean of identifications is \(2^{\varphi\langle e\rangle}\). The generating function with weights \(r^{d-1}\) is the same tail that produces \(T=r/(1-r)\). Sum: \(e_1=1\), and for \(d\ge 2\) set \(j=d-1\),

\[
\sum_{d=2}^{\infty} r^{d-1}\frac{d}{d-1}=\frac{r}{1-r}-\ln(1-r),
\]

so the weighted mean is \(1-(1-r)\ln(1-r)\). A 4D cutoff (\(325/264\)) is a fake wall. Putting the 2D *mix* \(12/11\) in place of Euclidean \(e_2=2\) and sending \(d\to\infty\) recovers planar leftover — that is Lemma 4, not infinite spacetime.

**Uniqueness.**

| Rival | Why it is not the theorem |
|---|---|
| Clock \(2^{\varphi}\) only | Flat \(N\to\infty\) of \(d/(d-1)\), omits the tail weights |
| Planar \(2^{\varphi\cdot 12/11}\) only | Lemma 4: 2D vacuum, not infinite spacetime |
| Finite mean \(325/264\) | 4D cutoff. Not the metric |
| Inverse \(11/12\), \(2/3\), \(3/4\) | Direction lemmas |
| \(f=1/(1-T/N)\) | \(N\) is not forced; partitions \(T\) |
| Solve \(f=H_{0,\mathrm{SN}}/H_{0,E}\) | \(H_0\) is a reading. Circular |
| Higher-\(D\) **density** \(T=1/(d-1)\) | Density stays the 2D mix \(49/71\) |

**Independence.** Adding a completed octave leaves \(\varphi\) invariant. The formula uses \(z_*\) and \(r\). It does not use \(H_0\) or a \(\chi^2\).

## Numerical lock (not a proof input)

G5 takes recombination from CAMB at the locked \(\Omega_m=22/71\), which agrees with the Planck catalog:

\[
z_*=1089.84\qquad\text{(CAMB lock; Planck \(1089.80\pm 0.21\))}.
\]

Then \(\varphi=\log_2(1090.84)-10\), \(\langle e\rangle\approx 1.310514\), and

\[
f=2^{\varphi\langle e\rangle}\approx 1.086396.
\]

Clock \(1090.84/1024=1.0652734375\). Planar \(2^{\varphi\cdot 12/11}\approx 1.071415\). 4D cutoff \(2^{\varphi\cdot 325/264}\approx 1.080952\).

## Predictions (after the theorem)

Same continuum, two readings (G4):

\[
H_{0,L}=f\,H_{0,E},\qquad r_{d,L}=r_{d,E}/f,\qquad t_U=t_{\mathrm{FRW}}(H_{0,L})\times f=t_{\mathrm{FRW}}(H_{0,E}).
\]

Time delays in this octave read the local frame: \(H_{0,\mathrm{TD}}=H_{0,L}\).

These are tests of the identification, not steps in the proof.

**Two measures, one \(\varphi\).** Mixing \(f_\infty\) onto BAO (or \(f_\square\) onto \(\theta_*\)) is the hybrid that dies on \(\chi^2\). The observables do not share a leftover measure:

- BAO \(D_M,D_H,D_V\) live in spatial slices of the 2D vacuum hierarchy. Leftover already sits on the 11 generators, so the leftover-period tail is the same series as \(T\) with ratio \(r/11\):

\[
f_\square=2^{\varphi\cdot (12/11)/(1-r/11)}=2^{\varphi\cdot 1440/1271}.
\]

Bare \(12/11\) is the lemma without that tail. Rival \(r/12\) would be clock-unit tail (wrong direction after Lemma 3). Do not tail \(f_\infty\) again (already an infinite sum).

- CMB \(\theta_*\) is last-scattering in infinite-D spacetime: \(f_\infty=2^{\varphi[1-(1-r)\ln(1-r)]}\).

SN tests both mapped locals. Density stays \(49/71\).

## Theorem 3 (Hubble scale)

Under Theorem 2, two readings (D), and the measure rule (C):

\[
H_{0,L}=f_\square\,H_{0,E}^{\mathrm{BAO}}=f_\infty\,H_{0,E}^{\theta_*}.
\]

SN tests both locals. \(f_\square\) and \(f_\infty\) are Theorem 2; they do not contain \(H_0\).

**Corollary (early–early).** The two early Hubble readings differ by the leftover-measure ratio, not by a third scale:

\[
\frac{H_{0,E}^{\mathrm{BAO}}}{H_{0,E}^{\theta_*}}=\frac{f_\infty}{f_\square}.
\]

Dumping DESI \(H_0\) into CAMB is two measures on one pair (the hybrid). Leftover does not eat that mixed frame. Official \(C_\ell\) and DESI full covariance are tests, not this proof.

Numerical scores: `src/hubble.py`, `scripts/run_hubble.py`. Inverse maps and the forward CAMB prediction at \(H_0=H_{0,\mathrm{SN}}/f_\infty\) are tests of the identification.

## What the theorem is not

- Not a derivation of \(z_*\). Recombination is G5.
- Not a derivation of \(N\), and not \(710/661\).
- Not the claim that leftover equals the empirical Hubble ratio. That is an empirical question, scored separately (`scripts/run_falsify.py`).

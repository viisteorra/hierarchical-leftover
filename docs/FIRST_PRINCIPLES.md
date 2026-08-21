# From first principles

No \(H_0\). No \(\chi^2\). This is the path from a resonant vacuum to every locked number, and the rules for taking it further.

Proofs: `THEOREM.md`. Uniqueness lemmas: `AXIOMS.md`. Spine: `THEORY.md`. Code: `src/uniqueness.py`, `src/generate.py`.

## Why the vacuum has this geometry

The one primitive **P** is not “dark energy equals 0.69.” It is:

> The vacuum *resonates*. It is an octave-addressed hierarchy. The residual tail of that hierarchy *is* the present dark-energy density.

Resonance is the reason the vacuum assumes this geometry. A resonator is not a bag of free parameters. It has:

- a **period** (what comes back as the same class),
- a **tonic** (one identity per period),
- **overtones** (integer ratios),
- a **tail** (what is left after the period repeats).

If that is the vacuum, then GR’s \(T_{\mu\nu}\) on the dark-energy slot is that tail. Maths cannot extract P from Einstein’s equation or from QFT (the cosmological-constant problem is a different number). P is the dictionary. Everything below is forced *given* P, or it is a measurement (recombination), or it is a different theory.

## The walk

### 1. Octave \(2:1\)

Elsewhere: linear waves. The first overtone \(2f\) is the unique interval that is the same pitch class. If the vacuum is harmonic, its period is that octave.

Not \(e\) (that is a calculus of growth, not a resonator). Not \(10\) (that is a counting base). Base \(2\) is the resonance.

### 2. Fifth \(3:2\)

The next integer after \(2\) is \(3\). The just fifth \(3:2\) is the next independent small-integer ratio. A resonator that knows octaves and knows the next harmonic knows fifths.

### 3. Period \(12\) — not a cents fit

Closing fifths against octaves is \(\log_2(3/2)\). Continued-fraction convergents:

\[
\frac{1}{1},\;\frac{1}{2},\;\frac{3}{5},\;\frac{7}{12},\;\frac{24}{41},\;\frac{31}{53},\;\ldots
\]

Convergents are the unique best rational approximations. Semi-convergents such as \(4/7\) are discarded. The period is the denominator of the first convergent *after* \(3/5\): uniquely **\(7/12\)**.

No cents threshold is an input. (Check, not the proof: every naive cut in \((1.96, 16.2)\) cents also yields \(n=12\). The figure \(2\) sits in that basin; it is not fitted.)

A \(53\)-fold is a *better* fifth. That is an analogous theory (keep P, swap U1). It is not this one. \(12\) is the smallest equal fold that carries a near-just fifth.

### 4. Squares \(q=4,5\) — not triangles, not a fit

Leftover of the metric uses Euclidean hypercubes, coordination \(q=2d\). In \(2\)D that is the square: \(q=4\). Regular tilings \(\{4,q\}\): Euclidean iff \(q=4\); the first hyperbolic square tiling is \(q=5\). Unique boundary pair.

Hyperbolic \(q=5\) lives in the **2D vacuum mix** (density). It is not FRW spatial curvature. Metric leftover is Euclidean \(\Rightarrow k=0\).

Generation multiplier \(r=2/q\): one octave (\(2\)) per coordination \(q\). That is the definition of the hierarchy, not a fit.

### 5. One seed \(\Rightarrow 1:11\) — not \(6:6\)

A period has one identity (the tonic / Euclidean seed). One identity in \(12\) parts \(\Rightarrow n_4=1\), \(n_5=11\). Unique.

Other mixed fills \(2:10,\ldots,11:1\) exist. They are **other axioms**. They are not selected from \(\chi^2\). SN’s own \(\Omega\) minimum is \(\approx 0.669\), not \(49/71\); the lock is not a supernova fit.

### 6. The tail is \(G\), not a second idea

The *continuous* residual of binary scaling is the Haar measure of one doubling:

\[
\Omega_{DE}=T=\ln 2,\qquad r=\frac{\ln 2}{1+\ln 2},\qquad
G(\rho)=\frac{1}{1-\rho}=\sum_{k=0}^{\infty}\rho^k,\qquad T=G(r)-1.
\]

The 12-fold / 1:11 mix is the unique smallest rational approximation of that same resonance (U1–U3):

\[
r_{12}=\frac{49}{120},\qquad T_{12}=\frac{49}{71}\approx 0.69014=\ln 2-0.0030.
\]

\(49=7^2\) is not a reason and is not used. Flatness (\(k=0\)) closes \(\Omega_m=1-\ln 2\). The lock is \(\ln 2\), not \(49/71\). Neither is chosen from \(\chi^2\).

If the period is a constant of the hierarchy, it does not drift, so the residual density does not dilute: \(w=-1\).

### 7. Leftover — mismatch of two readings, not \(H_0\)

A scale-invariant octave has no intrinsic count of octaves between “then” and “now.” Leftover needs a *second scale*. BAO and \(\theta_*\) share the sound-horizon ruler, frozen at photon decoupling. That freeze-out is the second scale. Not \(H_0\). Not \(z_{\mathrm{eq}}\). Drag \(z_d\) is a foil.

Integer octaves already live in \(a(t)\). Leftover is the fractional part

\[
\varphi=\{\log_2(1+z_*)\},\qquad f=2^{\varphi e}.
\]

Same \(G\), two licensed ratios only:

| Generating set | Why | Exponent |
|---|---|---|
| Finite: the \(11\) hyperbolic generators | Leftover of the 2D vacuum mix lives on the generators, not the clock (Lemma 3) | \(e_\square=(12/11)\,G(r/11)\) |
| Infinite: dimensions of the metric | The same tail, unbounded, on Euclidean hypercubes \(q=2d\) | \(e_\infty=1+(1-r)\ln G(r)\) |

\(G(r/12)\) is the clock-unit foil (wrong direction after Lemma 3). \(T/N\) (including \(T/10\)) is not an evaluation of \(G\). \(10\notin\{1,2,4,5,11,12\}\).

Which probe reads which set is the observable, not a fit: spatial slices (BAO, lensing) \(\to\) finite; null angle (\(\theta_*\)) \(\to\) infinite; this octave (SN, time delays) \(\to f=1\).

Two readings of **one** FRW metric:

\[
H_{0,L}=f\,H_{0,E}.
\]

\(H_0\) is a reading. Putting \(H_0\) into \(f\) is circular.

## Why these values are not arbitrary

| Temptation | Why it is not this theory |
|---|---|
| Pick \(12\) because it fitted \(\Omega_{DE}\) | \(12\) is U1 from \(\log_2(3/2)\), before any catalogue |
| Pick \(1:11\) from \(\chi^2\) | One seed. SN prefers \(\Omega\approx 0.669\), not \(49/71\) |
| Pick \(f\) so that \(H_{0,\mathrm{SN}}/H_{0,E}\) matches | \(f=2^{\varphi e}\) from \(z_*\) and \(r\) only |
| Use \(T/10\) because it is close | \(N=10\) is not in the hierarchy; \(G\) does not partition \(T\) |
| Nest \(11^2\), scan \(n_4:n_5\), average dimensions of density | Other axioms, or they worsen DESI vs Planck |
| One leftover \(f\) for every probe | Finite set and infinite set are different evaluations of \(G\) |
| Dump DESI \(H_0\) into CAMB | Two measures on one pair |

An analogous theory is allowed: keep P, swap a uniqueness lemma (e.g. period \(53\)). Then \(r\), \(T\), and leftover **copy** with the new integers. That is how you take it further without fishing.

## How to take it further

Do **not** add a parameter and minimise \(\chi^2\). Do this:

1. **New probe.** Classify its generating set: finite (spatial slice), infinite (null metric), none (this octave), or early microphysics (Einstein frame). Then \(e\) is already fixed. Weak lensing \(\to e_\square\). BBN stays early. Full \(C_\ell\) tests the \(\theta_*\) reading; it is not a third leftover.
2. **Analogous resonator.** Change one uniqueness lemma (period \(41\) or \(53\); a different seed rule). Recompute \(r\), \(T\), \(G\), leftover. Score the same kill tests. Do not mix lemmas to eat a residual.
3. **Kill tests.** Licensed BAO and \(\theta_*\) maps must stay \(<3\sigma\); licensed BAO \(\Delta\chi^2<10\); swapped measures must be worse; \(f\) must not depend on \(H_0\). `scripts/run_nails.py`, `run_falsify.py`, `run_hell.py`.

If a new number is not in \(\{1,2,4,5,11,12\}\), \(G\), \(\varphi\), or a measured second scale (\(z_*\), \(\omega_b\), \(T_{\mathrm{CMB}}\)), it is not this hierarchy.

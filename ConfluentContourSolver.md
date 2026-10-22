# The Confluent Contour Solver (CCS)

The **Confluent Contour Solver (CCS)** is a Lambert solver built on the contour-integral formulation of Negrete and Abdelkhalik:

> A. Negrete and O. Abdelkhalik, *Exact Solution to Lambert's Problem Using Contour Integrals*, Journal of Guidance, Control, and Dynamics 48(4), 2025, DOI 10.2514/1.G008499.

It keeps the paper's central idea: the time-of-flight root is a pole of $1/(T-t)$, located by a ratio of contour integrals. It changes nearly everything about how that idea is evaluated:
- the unknown it solves for;
- the form of the time function;
- where and how large the contours are;
- how many nodes they use;
- how the multi-revolution and hyperbolic cases are handled.

Most importantly, wherever it can, CCS replaces the sampled contour with its exact **zero-radius (confluent) limit**, computed from derivatives at one real point.

This document describes the paper's method (Section 1), lists every difference (Section 2), explains the mathematics behind CCS (Section 3), walks through the algorithm in full (Section 4), and reports measured accuracy and speed (Section 5).

- Implementation: `fortran/confluentcontour.f90`, module `confluentcontour`, entry points `contour_one` and `contour_batch`. This is the project's solver ID 0 ("Contour").
- Derivations and experiments: `research/REPORT.md`.

---

## Contents

1. [The method in the paper](#1-the-method-in-the-paper)
2. [What CCS does differently](#2-what-ccs-does-differently)
3. [Mathematical foundations of CCS](#3-mathematical-foundations-of-ccs)
4. [The CCS algorithm in detail](#4-the-ccs-algorithm-in-detail)
5. [Accuracy and performance](#5-accuracy-and-performance)
6. [Design alternatives that were tested and rejected](#6-design-alternatives-that-were-tested-and-rejected)
7. [Implementation map](#7-implementation-map)

---

## 1. The method in the paper

### 1.1 Universal-variable formulation (paper Eqs. 1–7)

The paper uses Bate's universal variable $z$:
- $z=\Delta E^2$ (the square of the eccentric-anomaly change) on elliptic arcs;
- $z<0$ on hyperbolic arcs;
- $0<z<4\pi^2$ for zero revolutions, and $(2n\pi)^2\le z<(2(n+1)\pi)^2$ for $n$ revolutions.

$$
\cos\phi=\frac{\mathbf r_1\cdot\mathbf r_2}{r_1r_2},\qquad
A=\frac{\sin\phi\,\sqrt{r_1r_2}}{\sqrt{1-\cos\phi}},
$$

$$
C(z)=\frac{1-\cos\sqrt z}{z},\qquad S(z)=\frac{\sqrt z-\sin\sqrt z}{\sqrt z^{\,3}},
$$

$$
y(z)=r_1+r_2+A\,\frac{zS(z)-1}{\sqrt{C(z)}},\qquad \chi(z)=\sqrt{\frac{y(z)}{C(z)}},
$$

$$
\sqrt\mu\,\Delta t(z)=\chi^3S+A\sqrt{y}.
$$

Once $z$ is known, the velocities come from the Lagrange coefficients:

$$
f=1-\frac{y}{r_1},\quad g=A\sqrt{\frac{y}{\mu}},\quad \dot g=1-\frac{y}{r_2},\qquad
\mathbf v_1=\frac{\mathbf r_2-f\mathbf r_1}{g},\quad \mathbf v_2=\frac{\dot g\,\mathbf r_2-\mathbf r_1}{g}.
$$

### 1.2 The root as a pole (Eqs. 8–11)

Define

$$
f(z)=\frac{1}{\chi^3S+A\sqrt y-\sqrt\mu\,\Delta t}.
$$

A simple root $z_*$ of the time equation is a simple pole of $f$. If $C$ encloses $z_*$ and no other singularity, and $f$ is analytic on and inside $C$ apart from $z_*$, then $(z-z_*)f(z)$ is analytic and Cauchy's theorem gives

$$
z_*=\frac{\oint_C z\,f(z)\,dz}{\oint_C f(z)\,dz}. \tag{Eq. 9/11}
$$

### 1.3 Quadrature (Eqs. 12–14)

On the circle $C(x)=z_0+C_re^{2\pi ix}$ the ratio becomes

$$
z_*=z_0+C_r\,\frac{\int_0^1e^{4\pi ix}f_C(x)\,dx}{\int_0^1e^{2\pi ix}f_C(x)\,dx}.
$$

- Both integrals are evaluated with the composite trapezoidal rule, which converges exponentially for periodic analytic integrands.
- Conjugate symmetry about the real axis halves the work (Eq. 14).
- The paper recommends $N=6$ points.

### 1.4 Contours (§III.B, §III.C, §IV)

| Case | Contour in the paper |
|---|---|
| Zero-revolution elliptic ($\Delta t>t_p$) | One large circle: center $2\pi^2$, radius $2\pi^2-\varepsilon$ ($\varepsilon=0.5$), covering the whole elliptic interval |
| Accuracy refinement (optional) | A second circle of radius 2 centered on the large-contour answer. This doubles the cost; it gave about $10^{-13}$ relative error at $N=6$ |
| Hyperbolic ($\Delta t<t_p$) | A lower bound $z_l$ on the hyperbolic root from a 5-term series reversion of $\Delta t(z)$ about $z=0$; circle centered at $z_l/2$, radius slightly less than $\lvert z_l\rvert/2$ |
| $n$ revolutions | A circle covering $(2n\pi)^2\dots(2(n+1)\pi)^2$; four moments $I_k=\oint z^kf\,dz$, $k=0..3$ |

### 1.5 Multi-revolution roots (Eqs. 19–25)

Both $n$-revolution roots $z_0,z_1$ lie inside one contour. Requiring $\oint(z-z_0)(z-z_1)f\,dz=0$ and $\oint z(z-z_0)(z-z_1)f\,dz=0$ gives

$$
(I_1^2-I_0I_2)\,z^2+(I_3I_0-I_1I_2)\,z+(I_2^2-I_1I_3)=0. \tag{Eq. 24}
$$

- Real roots mean a feasible transfer.
- A complex-conjugate pair means $\Delta t<t_{\min}$; the paper uses this as its feasibility test.

### 1.6 What the paper does not do

- **No root refinement:** the contour value is the answer.
- **No error estimate.**
- **No closed-form initial estimate:** the large contour plays that role.
- **No explicit map of where $f$ fails to be analytic,** beyond noting the pole at $z=4\pi^2$ and the non-analytic point at $\Delta t=0$.
- **Reported cost:** about 2.9 µs for one contour and 5.9 µs for large plus small (C++, Intel Xeon W-1370, $N=6$).

---

## 2. What CCS does differently

| # | Aspect | Paper | CCS |
|---|---|---|---|
| 1 | **Unknown** | $z$ | The half-anomaly angle $v=\sqrt z/2-m\pi\in(0,\pi)$ on every elliptic revolution interval. Near singular limits it switches to $w=\pi-v$, $w=\sqrt{-z}/2$ (hyperbolic) or $Y=\sqrt y$ (very short hyperbolic arcs) |
| 2 | **Time function** | Stumpff functions $C,S$; $\chi=\sqrt{y/C}$, $\sqrt{C}$, $\sqrt y$ | $T=\sqrt y\,(yF+a)$ with $y=y_0+b(1-\cos v)$ and $F=(v+m\pi-\sin v\cos v)/(\sqrt2\sin^3v)$. No Stumpff functions and one square root |
| 3 | **Analytic structure** | Assumed, with a contour shrunk by $\varepsilon$ | Mapped exactly: closed-form branch points of $\sqrt y$; the poles of $T$ are *zeros* of $1/(T-t)$, so they are harmless; nearby competing roots are modeled explicitly |
| 4 | **Starting information** | None; a contour covers the whole admissible interval | Closed-form initial estimates (inverse Hermite models), refined by about 66 KB of compiled-in tables in the Lancaster–Blanchard parameter (§4.4a) |
| 5 | **Root extraction (main path)** | Sampled contour ratio, $N=6$, complex arithmetic | **Confluent contour step:** the exact $r\to0$ limit of the paper's $N=4$ ratio, computed from $T,T',T'',T'''$ at one real point in closed form (Householder order 3), for elliptic, hyperbolic and near-parabolic arcs. Real arithmetic only; usually one step |
| 6 | **Root extraction (other regimes and fallback)** | — | The paper's ratio on a *small, sign-verified* contour around the estimate, with adaptive $N\in\{4,6,8,12,16\}$ and a built-in error estimate |
| 7 | **Need to enclose the root** | Required | Not required (Theorem 1). Enclosure is verified only to *size* the fallback contour |
| 8 | **Hyperbolic lower bound** | 5-term series reversion (approximate) | Closed form: the branch point $z_b=-16\,\operatorname{asinh}^2\!\sqrt{y_0/2b}$ |
| 9 | **Multi-revolution** | One large contour, four moments, quadratic Eq. 24 | Each branch solved separately on its own half-interval, seeded from a table of both branches. A rigorous closed-form lower bound on $t_{\min}$, a tabulated $t_{\min}$, and near $t_{\min}$ the exact minimum from the zero of an entire function $K(v)$ |
| 10 | **Feasibility test** | Complex roots of Eq. 24 | $t$ below the bound or clearly below the tabulated $t_{\min}$ → infeasible; within $10^{-4}$ of it, compare with the exact $t_{\min}$. Merged roots are reported explicitly |
| 11 | **Error control** | None | Order-2/order-3 agreement at the rounding level (confluent path; counts as verified), $N/2$-subsample estimate plus real residual check $\lvert T-t\rvert\le10^{-10}t$ (contour path), and a certified interval safeguard |
| 12 | **Near-parabolic** | Same contour | Exact $z=0$ at $t=t_p$; the linear term $z=(t-t_p)/T'(0)$ within rounding of $t_p$; parabolic series in $F$ for $\lvert v\rvert<0.3$ |
| 13 | **Velocity recovery** | $f,g,\dot g$ | Radial/transverse components from the half-angle pair, written without cancellation near $0$ and $180^\circ$ |
| 14 | **Geometry** | Raw $r_1,r_2,\cos\phi$ | Normalized by $\max(r_1,r_2)$; $\lvert\cos\frac\phi2\rvert$ and $\lvert\sin\frac\phi2\rvert$ from $\lVert\hat e_1\pm\hat e_2\rVert$. Collinear transfers are detected and flagged |
| 15 | **Status reporting** | Not specified | 0 success, 1 merged roots, 5 infeasible revolution count, −1 invalid, −2 undefined plane, −4 no accurate root |
| 16 | **Cost** | 2.9–5.9 µs (C++) | 0.05–0.27 µs per solve depending on regime, about 0.195 µs on average (Fortran, Section 5) |

The following sections explain each change.

---

## 3. Mathematical foundations of CCS

### 3.1 The half-anomaly form of the time equation

Put $u=\sqrt z/2$ (half the eccentric-anomaly change). The Stumpff functions collapse into elementary trigonometry:

$$
C(z)=\frac{\sin^2u}{2u^2},\qquad zS(z)-1=-\frac{\sin 2u}{2u},\qquad S(z)=\frac{2u-\sin2u}{8u^3}.
$$

Taking the analytic square root $\sqrt{C}=\sin u/(\sqrt2\,u)$, which is entire in $z$:

$$
y=r_1+r_2-\sqrt2A\cos u,\qquad \chi^3S=y^{3/2}\,\frac{u-\sin u\cos u}{\sqrt2\,\sin^3u}.
$$

These identities were checked to 40 digits at complex $z$.

On the $m$-revolution interval set $v=u-m\pi\in(0,\pi)$. Normalize lengths by $s=\max(r_1,r_2)$ and time by $\sqrt{s^3/\mu}$, and write $a$ for the normalized $A$ and $b=\sqrt2a$. Then

$$
\boxed{T(v)=\sqrt{y(v)}\,\big(y(v)\,F(v)+a\big),\qquad y(v)=y_0+b(1-\cos v),\qquad F(v)=\frac{v+m\pi-\sin v\cos v}{\sqrt2\,\sin^3v}}
$$

with $y_0=y(0)=r_1+r_2-b$. The same expression covers every revolution count, and each interval has the same width $\pi$. The sign of $\sin u$ on odd revolutions is absorbed by working in $v$.

The universal-variable identities themselves are classical (Lagrange/Lancaster–Blanchard structure). CCS's use of them is computational: one sine/cosine pair and one square root per evaluation instead of nested Stumpff functions.

### 3.2 Where $1/(T-t)$ stops being analytic

Trapezoidal contour accuracy, and the reach of the confluent step, depend on the nearest point where $1/(T-t)$ is not analytic. In the $v$ plane there are exactly three kinds:

1. **Branch points of $\sqrt y$.** $y=0$ when $\cos v=1+y_0/b$. In Lancaster–Blanchard form these are $v=\pm i\eta$ (short way, $b>0$) or $v=\pi\pm i\eta$ (long way), with $\eta=-\ln\lvert\lambda\rvert$. They never lie on the real elliptic interval, and they lie outside the disk $\lvert v-\pi/2\rvert<\pi/2$ for every geometry. On the imaginary axis this is exactly the paper's "$\Delta t=0$" non-analytic point: CCS has it in closed form, $z_b=-16\,\operatorname{asinh}^2\sqrt{y_0/(2b)}$, whereas the paper approximates it by series reversion.
2. **The poles of $T$** at $v=k\pi$ (removable at $v=0$ when $m=0$). These are *zeros* of $1/(T-t)$, not singularities. This is why the paper's reciprocal choice is good, and why a log-derivative $T'/(T-t)$ is bad: it turns them into poles.
3. **Other solutions of $T(v)=t$** are additional poles of $1/(T-t)$. These are what actually limit convergence:
   - the twin real root on multi-revolution intervals;
   - the mirror root $-v_*$ for $m=0$ ($T$ is even in $v$);
   - near $v=\pi$, where $T\sim K/(\pi-v)^3$, the two complex companions of the real root at $120^\circ$, a distance $\sqrt3\,(\pi-v_*)$ away.

### 3.3 Theorem 1: what the trapezoidal ratio really computes

Write $1/(T-t)=R/(x-p)+h(x)$ near the root $p$, where:
- $R=1/T'(p)$ is the residue;
- $h$ is analytic in $\lvert x-c\rvert<\rho$, with $\rho$ the distance from $c$ to the nearest point of §3.2;
- $\delta=p-c$.

Take $N$ trapezoidal nodes $c+r\omega_j$, $\omega_j=e^{2\pi ij/N}$, with $0<r<\rho$ and no node at $p$. The paper's estimate is

$$
\hat p_N(r)=c+r\,\frac{\hat f_{-2}}{\hat f_{-1}},\qquad \hat f_{-k}=\frac1N\sum_{j}\omega_j^{\,k}\,\frac{1}{T(c+r\omega_j)-t}.
$$

**Theorem 1.** *Exactly,*

$$
\hat p_N(r)-p=\frac{(r^N-\delta^N)\,E(r)}{R+(r^N-\delta^N)\,E'(r)},
\qquad
E=\sum_{l\ge1}\big(h_{lN-2}-\delta\,h_{lN-1}\big)r^{(l-1)N},\quad
E'=\sum_{l\ge1}h_{lN-1}\,r^{(l-1)N},
$$

*where $h_k$ are the Taylor coefficients of $h$ at $c$.*

*Proof sketch.*
- **Pole part.** Its discrete coefficients are $\hat f_{-k}=R\,\delta^{k-1}r^{-k}/(1-q)$ with $q=(\delta/r)^N$. This one formula holds whether the pole is inside ($\lvert\delta\rvert<r$) or outside ($\lvert\delta\rvert>r$) the circle; it is the same rational function of $\delta$. Its contribution to the ratio is therefore exactly $\delta$.
- **Regular part.** It aliases as $\hat h_{-k}=\sum_l h_{lN-k}r^{lN-k}$.
- Substituting both into the ratio and simplifying gives the formula. $\square$

The formula was verified against 50-digit computations to 6 significant digits, for real and complex $\delta$.

**Consequences for the paper's design:**

| Consequence | Meaning |
|---|---|
| (a) Exact for any pure pole, at any radius | Only the regular part $h$ causes error |
| (b) **Enclosure is not needed** | A circle that misses the root is just as valid. The paper's test "the denominator evaluates to zero if the pole is outside" holds only for exact integrals: the discrete $\hat f_{-1}$ is small but nonzero, and the ratio stays accurate |
| (c) Big contours are expensive | For fixed $N$, choosing $r\gg\lvert\delta\rvert$ inflates the error by about $(r/\lvert\delta\rvert)^N$. This is why the paper's large first contour is inaccurate and its radius-2 refinement helps |
| (d) **Confluent limit** | As $r\to0$, $\hat p_N\to c+a_{N-2}/a_{N-1}$, where $a_k$ are the Taylor coefficients of $1/(T-t)$ at $c$. That is König's method, i.e. **Householder's method of order $N-1$**, with error exactly $-\delta^N(h_{N-2}-\delta h_{N-1})/(R-\delta^Nh_{N-1})$ |
| (e) Why not shrink a *sampled* circle | In floating point, a sampled circle with $r<\lvert\delta\rvert$ computes $\hat f_{-1}$ by cancellation, losing about $(\lvert\delta\rvert/r)^{N-1}$ in relative accuracy. Computing the Taylor coefficients analytically has no such cancellation |

So the paper's ratio and high-order Householder iteration are two members of one family. CCS uses the confluent member where it can compute derivatives cheaply, and the sampled member on small verified contours elsewhere.

**Two-root analog.** The paper's two-root quadratic (Eq. 24) built from trapezoidal moments is likewise exact for two pure poles at any radius. Each pole's moments stay geometric, with weights $R_i/(1-q_i)$. Its $r\to0$ limit is the Hadamard (Padé $[n/2]$) denominator built from Taylor coefficients.

### 3.4 An exact ODE for the time function

Since $\frac{d}{dv}\big(v+m\pi-\sin v\cos v\big)=2\sin^2v$:

$$
\frac{d}{dv}\big(\sqrt2\sin^3v\,F\big)=2\sin^2v
\quad\Longrightarrow\quad
\boxed{\sin v\,F'(v)+3\cos v\,F(v)=\sqrt2}
$$

This holds for every revolution count $m$; $m$ enters only through $F$ itself. Differentiating repeatedly gives all derivatives of $F$ from $s=\sin v$, $c=\cos v$ and a single reciprocal $1/s$:

$$
F'=\frac{\sqrt2-3cF}{s},\qquad
F''=\frac{3sF-4cF'}{s},\qquad
F'''=\frac{3cF+7sF'-5cF''}{s}.
$$

The other factors are elementary: $y'=bs$, $y''=bc$, $y'''=-bs$. So $T$ and its first three derivatives cost one sine/cosine pair (carried by rotation, §4.5), one square root and about 60 floating-point operations, all in real arithmetic.

---

## 4. The CCS algorithm in detail

### 4.1 Call sequence

```
contour_one(r1, r2, tof, direction, v1, v2, status, half, [nrev], [branch])
│
├─ 1. prepare_geometry                        normalize, half-angle data, plane check
├─ 2. t = tof / scale^(3/2)                   (mu = 1 units)
├─ 3. zero revolutions, a > 0, t < 0.1 tp  →  very short hyperbola in Y = sqrt(y)   (4.8)
└─ 4. solve_z
      ├─ m = 0, t = tp to rounding              →  z = 0 (parabola)
      ├─ m = 0, |t − tp| at rounding level      →  z = (t − tp)/T'(0)
      ├─ m = 0, t > tp  →  elliptic estimate    ─┐
      ├─ m = 0, t < tp  →  hyperbolic estimate   ├─→ root finder (4.5 confluent step,
      ├─ m > 0          →  feasibility + branch  ─┘      4.6 verified small contour)
      │                    estimate (4.7)
      ├─ residual check |T(z) − t| ≤ 1e−10 t   (contour paths only)
      └─ safeguards: interval contours, certified interval contour (with two-root moments)
└─ 5. reconstruct                              velocities (4.9)
```

### 4.2 Geometry

- $s=\max(\lvert\mathbf r_1\rvert,\lvert\mathbf r_2\rvert)$; $\hat r_i=r_i/s$; unit vectors $\hat e_i$; $\hat h=\pm(\hat e_1\times\hat e_2)/\lVert\cdot\rVert$, with the sign from the transfer direction $d=\pm1$; transverse unit vectors $\hat t_i=\hat h\times\hat e_i$.
- If $\lVert\hat e_1\times\hat e_2\rVert<32\varepsilon$ the plane is undefined: status −2. Below $10^{-4}$ the `half` flag reports near-collinear geometry.
- $c_s=\lvert\cos\frac\phi2\rvert=\frac12\lVert\hat e_1+\hat e_2\rVert$ and $s_s=\lvert\sin\frac\phi2\rvert=\frac12\lVert\hat e_2-\hat e_1\rVert$. Using norms avoids the cancellation in $1\pm\cos\phi$.
- $a=d\sqrt2\sqrt{\hat r_1\hat r_2}\,c_s$ (the paper's $A$, normalized) and $b=\sqrt2a$.
- $y_0=q^2+2\sqrt{\hat r_1\hat r_2}\,s_s^2/(1+c_s)$ for $d>0$, or $q^2+2\sqrt{\hat r_1\hat r_2}(1+c_s)$ for $d<0$, with $q=(\hat r_1-\hat r_2)/(\sqrt{\hat r_1}+\sqrt{\hat r_2})$. Also $y_{\text{end}}=y_0+2b$, written without cancellation.
- In terms of half angles, $y=y_0+2b\sin^2\frac v2$, or the equivalent $y_{\text{end}}-2b\cos^2\frac v2$ when $b<0$. Whichever form has no cancellation is used.
- Parabolic time and slope: $t_p=\sqrt{y_0}\big(\tfrac{\sqrt2}{3}y_0+a\big)$ and $T'(0)=\frac{y_0^{3/2}}{20\sqrt2}+\frac b8\big(\frac{\sqrt{y_0}}{\sqrt2}+\frac{a}{2\sqrt{y_0}}\big)$ (derivative with respect to $z$).

### 4.3 Coordinates

Each regime is solved in the coordinate where its formula has no cancellation and its singularities are far away:

| Kind | Coordinate | Used for | Why |
|---|---|---|---|
| 2 | $v=\sqrt z/2-m\pi$ | Zero-rev elliptic with $0.49<z\le2\pi^2$; the low-$z$ multi-rev branch | Avoids a complex $\sqrt z$; sin/cos of $v/2$ feed everything |
| 1 | $w=\pi-v$ | Zero-rev elliptic with $z>2\pi^2$; the high-$z$ multi-rev branch | Near the $v=\pi$ pole the small quantity $w$ is carried exactly. Velocity recovery uses $q=4\pi^2-z=4w(2\pi-w)$ instead of the rounded $z$ |
| 0 | $z$ | Zero-rev elliptic near the parabola ($z\le0.49$); mildly hyperbolic | Parabolic series region |
| 4 | $w=\sqrt{-z}/2$ | Hyperbolic with $z<-1$ | Hyperbolic sine/cosine from one `expm1` |
| 5 | $Y=\sqrt y$ | Very short hyperbolic arcs, $t<0.1\,t_p$ | Never forms $y$ by cancellation (4.8) |
| 3 | $v$, root of $K(v)$ | Locating the multi-rev minimum time | $K$ is entire (4.7) |

Near $v=0$ on zero-revolution arcs, $F$ uses a 10-term parabolic series in $v^2$ (truncation below $3\times10^{-19}$ for $\lvert v\rvert<0.3$) instead of the cancelling closed form.

### 4.4 Closed-form initial estimates

**Zero-revolution elliptic.**
- Let $\varphi=T^{-1/3}$ and $q=4\pi^2-z$. Near the upper boundary $T\sim K_{\text{top}}/q^3$, so $q\propto\varphi$ there; at the parabola $T=t_p$ with known slope.
- An inverse cubic Hermite model of $q(\varphi)$ (or of $z$ in the lower half) matches the exact pole coefficient $K_{\text{top}}=8\pi^{4/3}2^{-1/6}\sqrt{y_0+2b}$ at $q=0$, and $t_p$ and $T'(0)$ at $z=0$.
- It costs two cube roots. It also gives the local slope $\lvert dx/d\varphi\rvert$, which predicts how far the estimate is from the root.

**Hyperbolic.**
- A cubic Hermite model of $z(T^2)$ between the branch point ($T=0$ at $z_b$, with the exact slope there) and the parabola ($t_p$, $T'(0)$).
- For $b<0$ (long-way hyperbolic) there is no branch point: $w^2=t_p\ln(t_p/t)/(4T'(0))$ near the parabola, and for $t<t_p/e$ the larger of that and the asymptote $w=2\ln\!\big(\sqrt{-b}\,(y_0+b)/t\big)$ from $T\sim\sqrt{-b}\,(y_0+b)\,e^{-w/2}$. Before this estimate (a linear model capped at $z\ge-4$), long-way arcs with $t<10^{-2}t_p$ mostly fell through to the safeguards at about 14 µs per solve; they now take about 210 ns.

**Very short hyperbolic** ($t<0.1t_p$): the same Hermite construction in $Y=\sqrt y$.

**Multi-revolution.**
- Write $T=H/(\sqrt2\sin^3v)$ with $H=\sqrt y\,(y(v+m\pi-\sin v\cos v)+\sqrt2a\sin^3v)$. Then $\varphi=T^{-1/3}=\sin v\cdot g(v)$, where $g=(\sqrt2/H)^{1/3}$ is smooth and positive.
- $\log g$ is modeled by a quadratic through $v=0$ and $v=\pi$ (exact end values) and one exact value at the model's minimum-time point $v_m$.
- The branch estimate solves $\sin v\,e^{\ell(v)}=\varphi_t$ by two fixed-point iterations, $v\leftarrow\arcsin\!\big(\varphi_t\,e^{-\ell(v)}\big)$.
- Near a located minimum $v_*$, the estimate instead comes from the quadratic with zero slope at $v_*$ through the current point: $v=v_*+(v_0-v_*)\sqrt{(t-t_{\min})/(T(v_0)-t_{\min})}$.

These estimates typically land within $10^{-2}$–$10^{-3}$ of the root in the working coordinate. Inside the tables of §4.4a they are replaced by table values that are usually within $10^{-5}$.

### 4.4a Table-refined estimates (compiled in, about 66 KB)

The closed-form estimates of §4.4 usually land $10^{-2}$–$10^{-3}$ from the root, which costs a second confluent step. CCS therefore refines them with small tables compiled into the source (no data file). They sit in `fortran/contour_table_zero_rev.inc` and `fortran/contour_table_multi_rev.inc`, generated by `tools/gen_table.py` and `tools/gen_mrev.py`.

**Why two variables suffice.** In Lancaster–Blanchard form the root depends only on
- the geometry parameter $\lambda=d\,\sqrt{r_1r_2}\cos\frac\phi2/s_L$, where $s_L=(r_1+r_2+c)/2$ is the semiperimeter and $c$ the chord ($1-\lambda^2=c/s_L$ exactly), and
- one normalized time.

Times scale with $s_L^{3/2}$. Tables use the stretched coordinate $\tau=\lambda/(1+\sqrt{1-\lambda^2})$, which costs one division. It is $\tan(\tfrac12\arcsin\lambda)$, so it packs grid points toward the degenerate limits $|\lambda|\to1$, where the roots change fastest.

**Zero revolutions.** The time coordinate is $s=(t_p/t)^{1/3}\in(0,1)$, already computed for the closed-form estimate. Two tables store functions that stay smooth through both ends:

| Table | Range | Stored | Behavior it absorbs |
|---|---|---|---|
| A | $s\le\tfrac12$ | $(\pi-v)/s$ | $\pi-v\propto s$ as $t\to\infty$ |
| B | $s>\tfrac12$ | $z/(1-s)$ | $z\propto 1-s$ at the parabola |

The size is $32\times32$ per half (16 KB).

**Multiple revolutions.**
- **Minimum time.** $t_{\min}^{1/3}$ (stored as $\sqrt[3]{t_{\min}/m}/\sqrt{s_L}$) and the location $v_*$ of the minimum are tabulated on a fine 128-point $\tau$ grid. They are accurate to about $4\times10^{-7}$ relative in $t_{\min}^{1/3}$; a 256-point grid gives $2\times10^{-8}$.
- **Branch roots.** Both branches lie on one smooth curve in $\sigma=\text{branch}\cdot\sqrt{1-(t_{\min}/t)^{1/3}}\in(-1,1)$. Near $t_{\min}$ the roots behave like $v_*\pm\kappa\sqrt{t-t_{\min}}$, which is analytic in $\sigma$; as $t\to\infty$, $v$ and $\pi-v$ are proportional to $1-\sigma^2$. The table is $16\times24$ per slice.
- **Revolution count.** The roots vary smoothly in $\mu=1/m$, so ten slices ($m=1,2,3,4,6,8,12,16,32,\infty$) with 4-point Lagrange interpolation in $\mu$ cover every $m$. All multi-revolution tables together are 50 KB.

**Use at run time.**
- Every lookup is a local $4\times4$ (or $4\times4\times4$) cubic Lagrange interpolation of stored values: about 30–100 flops, with no transcendental functions beyond the $\sqrt{\ }$ and $\sqrt[3]{\ }$ the estimate needs anyway.
- For multiple revolutions, $r=(t_{\min}/t)^{1/3}$ is formed from the table, and:
  - $r>1+10^{-4}$: status 5 (infeasible) immediately;
  - $|r-1|\le10^{-4}$ (within about $3\times10^{-4}$ of $t_{\min}$): the exact logic of §4.7 runs, including the minimum from $K(v)$;
  - otherwise the tabulated branch root seeds the confluent step, with $v_*$ as the branch boundary.
- The margin is about 250 times the measured table error. `tests/test_contour_tables.f90` checks the classification and velocities against the exact minimum time at $t=t_{\min}(1\pm\delta)$, $\delta=10^{-1}\ldots10^{-6}$ (against an independent golden-section minimum of the exact time function), for 150 random geometries and $m\in\{1,2,3,5,7,12,20,50,200\}$.
- Outside the tables ($|\tau|>0.995$, or $s>0.99$) or on any failure, the closed-form path of §4.4 runs unchanged.

**Effect.**
- Zero-revolution elliptic solves take about **one** confluent step (previously about two).
- Multi-revolution solves take about one step and rarely need the minimum-time search.
- Infeasible revolution counts are usually decided by one table lookup.
- Larger tables (up to 3.3 MB) were tested and give no further speedup (Section 6).

### 4.5 The confluent contour step (main root extraction)

It is used for coordinates 1 and 2 (every zero-revolution elliptic case with $z>0.49$ and both multi-revolution branches), and in the two forms described at the end of this section: hyperbolic, for coordinate 4 ($z<-1$), and near-parabolic, for coordinate 0 ($-1\le z\le0.49$). It runs before any bracket or contour set-up; that set-up happens only if the step declines.

**Inputs:**
- the estimate $x_0$;
- the half-angle pair $(\sin\frac{x_0}2,\cos\frac{x_0}2)$;
- the admissible interval $(lo,hi)$: $(0,\pi)$ for zero revolutions, or the branch's side of the minimum for multi-revolution.

**Per step:**

1. From the half-angle pair: $s=2\sin\frac x2\cos\frac x2$, $c=(\cos\frac x2-\sin\frac x2)(\cos\frac x2+\sin\frac x2)$. In coordinate $w=\pi-v$, $\cos v=-c$ and $u=(m+1)\pi-w$; in coordinate $v$, $u=v+m\pi$. $y$ comes from the cancellation-free half-angle form.
2. Derivatives with respect to $v$ (§3.4), using $\iota=1/s$:

   $$F=\frac{u-sc}{\sqrt2}\iota^3,\quad F_1=(\sqrt2-3cF)\iota,\quad F_2=(3sF-4cF_1)\iota,\quad F_3=(3cF+7sF_1-5cF_2)\iota$$

   $$y_1=bs,\quad y_2=bc,\quad y_3=-bs$$

   $$Q=\sqrt y,\quad Q_1=\tfrac{y_1}{2Q},\quad Q_2=\tfrac{y_2/2-Q_1^2}{Q},\quad Q_3=\tfrac{y_3/2-3Q_1Q_2}{Q}$$

   $$P=yF+a,\quad P_1=y_1F+yF_1,\quad P_2=y_2F+2y_1F_1+yF_2,\quad P_3=y_3F+3y_2F_1+3y_1F_2+yF_3$$
3. Residual and Taylor coefficients of $T=QP$ in the working coordinate, with $\sigma=+1$ for $v$ and $-1$ for $w$:

   $$g_0=QP-t,\quad t_1=\sigma(Q_1P+QP_1),\quad t_2=\tfrac12(Q_2P+2Q_1P_1+QP_2),\quad t_3=\tfrac\sigma6(Q_3P+3Q_2P_1+3Q_1P_2+QP_3).$$
4. Scaled Taylor coefficients $b_n=a_ng_0^{\,n+1}$ of $1/(T-t)$:

   $$b_1=-t_1,\qquad b_2=-(t_2g_0+t_1b_1),\qquad b_3=-(t_3g_0^2+t_2g_0b_1+t_1b_2).$$
5. Order-2 (Halley) and order-3 confluent values, i.e. the $r\to0$ limits of the 3- and 4-node trapezoidal ratios:

   $$x_2=x_0+g_0\frac{b_1}{b_2},\qquad x_3=x_0+g_0\frac{b_2}{b_3}.$$
6. Carry the half-angle pair to $x_3$ by an exact rotation: short Taylor polynomials in $(x_3-x_0)/2$, or one `sincos` call for large steps or near $0$ and $\pi$. No libm sine/cosine is needed on the ordinary path.

**Acceptance.**
- Let $D=\lvert x_3-x_2\rvert$ and $\Delta=\lvert x_3-x_0\rvert$. With geometric (König) convergence, $D$ estimates the order-2 error, and the order-3 error is about $D\,(D/\Delta)^{1/3}$.
- Accept $x_3$ when

  $$D\le4\varepsilon\lvert x_3\rvert\qquad\text{or}\qquad D<\Delta\ \text{ and }\ D^4\le\Delta\,(\varepsilon\lvert x_3\rvert)^3 .$$
- Otherwise repeat from $x_3$, up to **4 steps**. From a table estimate (error about $10^{-5}$ or less) one step usually reaches machine precision and passes this test. From a closed-form estimate (error $10^{-2}$–$10^{-3}$) two steps are usually needed.

**Decline**, and hand the same estimate to the verified contour (4.6), if any of these happen:
- $x_3$ leaves $(lo,hi)$;
- $y\le0$;
- the iteration is not contracting ($D\ge\Delta$);
- the point is in the zero-revolution parabolic-series region $\lvert v\vert<0.3$.

An accepted step counts as verified, so the separate real residual check is skipped. Over 220 000 benchmark cases that check never rejected a converged confluent root: the worst $\lvert T-t\rvert/t$ was $1.4\times10^{-14}$, against a threshold of $10^{-10}$. Contour fallbacks always run it.

**Why it is faster than the paper's contour:**
- One sine/cosine (by rotation) and one square root per step, in real arithmetic.
- A complex node needs a sine/cosine *and* an exponential, a complex square root and a complex division.
- By Theorem 1(d), one order-3 confluent step extracts as much as a 4-node contour around the estimate, with no sampling cancellation.

**Order 4 on elliptic arcs.** For coordinates 1 and 2 the step is order 4, the zero-radius limit of the 5-node ratio. Differentiating the ODE once more gives
$$sF''''=-3sF+10cF'+12sF''-6cF''',$$
with $y''''=-bc$, $Q''''=(y''''/2-4Q'Q'''-3Q''^2)/Q$ and the Leibniz rule for $P''''$ and $T''''$. Then $b_4=-(t_4g_0^3+t_3g_0^2b_1+t_2g_0b_2+t_1b_3)$ and $x_4=x_0+g_0b_3/b_4$. The acceptance test becomes $D=|x_4-x_3|$, $D^5\le\Delta(\varepsilon|x_4|)^4$. No extra transcendental calls are needed, and starting errors up to about $10^{-4}$ (the table estimates) finish in one step. Zero-revolution solves take exactly one step; multi-revolution solves take about 1.04.

**Hyperbolic form** (coordinate $w=\sqrt{-z}/2$). With $S=\sinh w$ and $C=\cosh w$, $F=(SC-w)/(\sqrt2S^3)$ satisfies the same ODE, $S\,F'+3C\,F=\sqrt2$. Hence
$$F''=-\frac{4CF'+3SF}{S},\qquad F'''=-\frac{5CF''+7SF'+3CF}{S},\qquad y=y_0-2b\sinh^2\tfrac w2,\quad y'=-bS,\ y''=-bC,\ y'''=-bS.$$
The step and acceptance test are unchanged. $\sinh$ and $\cosh$ come from one `expm1` per step, and the step declines for $w<0.3$, where the parabolic series applies.

**Near-parabolic form** (coordinate 0, $x=z/4=v^2$, $-\tfrac14\le x\le0.1225$). Here $T$ is a polynomial and one square root:
- $F(x)=\tfrac1{\sqrt2}\sum_{k=0}^{13}f_kx^k$ uses the exact series coefficients ($f_0=\tfrac23$, $f_1=\tfrac15$, $f_2=\tfrac{17}{420}$, …), with truncation below $3\times10^{-19}$ on this range;
- $1-\cos v=x\sum_{k\ge0}(-x)^k/(2k+2)!$, so $y=y_0+b\,(1-\cos v)$;
- $T,T',T'',T'''$ in $x$ come from one Horner pass with derivatives, followed by the same $Q$, $P$, Leibniz and step formulas;
- the same expression covers both elliptic ($x>0$) and hyperbolic ($x<0$) arcs, with no change of formula at the parabola.

If an iterate leaves the series range below $z=-1$, the hyperbolic form continues from it in $w=\sqrt{-z}/2$.

**Very short hyperbolic arcs** (coordinate 5, $Y=\sqrt y$, $t<0.1\,t_p$). With $S=\sinh^2\frac w2=(y_0-Y^2)/(2b)$, the function $G(S)=F(w(S))$ satisfies an ODE with polynomial coefficients,
$$2S(1+S)\,G'+3(1+2S)\,G=\sqrt2,\qquad G''=-\frac{5(1+2S)G'+6G}{2S(1+S)},\qquad G'''=-\frac{7(1+2S)G''+16G'}{2S(1+S)},$$
and $T=Y^3G(S(Y))+aY$, with $S'=-Y/b$ and $S''=-1/b$.
- For $S<0.0225$, where the closed form of $G$ cancels and the ODE divides by a small $S$, $G=\tfrac{\sqrt2}{3}\,{}_2F_1(1,3;\tfrac52;-S)$ is summed directly. That is 12 terms, coefficients from $c_{k+1}=-c_k(2k+6)/(2k+5)$, with derivatives by Horner.
- The root is $Y$, then $z=-16\operatorname{asinh}^2\sqrt{(y_0-Y^2)/2b}$ as before.

The routines keep no state between calls (all variables local, iteration limits compile-time constants), so batch calls under OpenMP are safe.

### 4.6 The verified small contour (other regimes and fallback)

It is used for coordinate 3 (the minimum-time function near $t_{\min}$) and whenever a confluent step declines.

1. **Sign bracket.** Step from the estimate by $1.25\times$ the predicted distance $\lvert T(x_0)-t\rvert\cdot\lvert dx/dT\rvert$, in the direction indicated by the residual's sign. If the sign has not changed, the first retry aims just past the secant root; after that the step grows by $4\times$ (at most 6 expansions). A sign change *proves* a root lies in the bracket. Because the endpoints are also the real nodes of the contour, this check costs nothing extra.
2. **Contour.** A circle centered on the bracket with radius $r$ equal to half its width.
3. **Node count.** With $\rho$ the modeled distance to the nearest point of §3.2, the smallest $N$ with $(r/\rho)^{2N-2}$ below about $10^{-14}$:

   | $\rho/r$ | $\ge10^{7/3}$ | $\ge10^{7/5}$ | $\ge10$ | $\ge10^{7/11}$ | otherwise |
   |---|---|---|---|---|---|
   | $N$ (semicircle intervals) | 4 | 6 | 8 | 12 | 16 |

   The circle has $2N$ points; conjugate symmetry leaves $N-1$ complex evaluations, plus the two real endpoint values already known.
4. **Node evaluation.** A branch-free loop that the compiler vectorizes. $T-t$ is formed as one fraction, $\big[\sqrt y\,(y\,n+\sqrt2a s^3)-\sqrt2t\,s^3\big]/(\sqrt2s^3)$ with $n=u-sc$, so each node needs one complex division. Sine/cosine of node angles come from the center's half-angle pair plus short Taylor polynomials in the small node offsets.
5. **Estimate and error.** The paper's ratio from all nodes gives $x$; the ratio from every other node (the $N/2$ rule) gives $x_h$. The error estimate is $e=\lvert x-x_h\rvert\cdot\max\!\big(\lvert x-x_h\rvert/r,\ \min(1,r/\rho)^N\big)$. Accept when $e\le2\varepsilon\lvert x\rvert$. Otherwise repeat from $x$ with a bracket of size $\min(4e,r)$, up to 4 passes.
6. **By-product.** The contour residue $r\cdot\text{den}/(2N)=1/T'(p)$ is the derivative of the root with respect to time of flight. CCS uses it to size later brackets and to obtain $T''$ at the multi-revolution minimum.

### 4.7 Multi-revolution logic

1. **Rigorous infeasibility test.** Exactly, $T=y^{3/2}h_m(v)/\sqrt2+a\sqrt y$ with $h_m(v)=(m\pi+v-\sin v\cos v)/\sin^3v$. A table of $\min_v h_m$ ($m\le64$, extended by $+\pi$ per revolution beyond) and the extreme values of $y$ give a lower bound on $t_{\min}$. If $t$ is below it: **status 5**, with no root finding at all.
2. **Table classification** (§4.4a). If $t$ is clearly below the tabulated $t_{\min}$: **status 5**. If it is clearly above, the tabulated branch root seeds the confluent step on $(0,v_*)$ or $(0,\pi-v_*)$, and the solve is usually done. Within $10^{-4}$ of $t_{\min}$, or outside the tables, continue with steps 3–5.
3. **Model minimum.** The $\log g$ model of §4.4 gives $v_m$ and an exact $T(v_m)$. If $t>T(v_m)$, one root lies on each side of $v_m$ and the branches are solved directly.
4. **Undecided case** ($t\le T(v_m)$).
   - Find the exact minimum $v_*$ as the root of the entire function $K=s\,(y'P+2yP')-6c\,yP$, where $P=y\,n+\sqrt2as^3$. This has $\operatorname{sign}K=\operatorname{sign}dT/dv$ and contains no square root. It is solved with the verified contour (kind 3).
   - Then $t_{\min}=T(v_*)$:
     - $t<t_{\min}$: **status 5**;
     - $t=t_{\min}$ to rounding: **status 1** (merged roots, returned at $v_*$);
     - otherwise $v_*$ splits the branches, and $T''(v_*)$ (from the contour residue of $K$) seeds the near-minimum quadratic estimate.
5. **Branches.** `branch = −1` solves in $v\in(0,v_*)$; `branch = +1` solves in $w=\pi-v\in(0,\pi-v_*)$. Each goes through the confluent step first (4.5), then the verified contour if needed (4.6). Solving each branch on its own half-interval removes the twin root from that solve; near $t_{\min}$ the twin sits just beyond the interval end.

Compared with the paper's single large contour and Eq. 24, this:
- never computes four weighted moments whose residues diverge with opposite signs near $t_{\min}$;
- decides feasibility from a bound or an exact minimum, not from the sign of a numerically computed discriminant;
- reports merged roots explicitly.

### 4.8 Hyperbolic and near-parabolic transfers

- **Mild hyperbolic** ($-1<z<0$): near-parabolic confluent steps in $x=z/4$ (§4.5), handing off to the hyperbolic form if the root lies below $z=-1$, with the verified contour in $z$ (kind 0, bracketed below by the branch point $z_b$) as fallback.
- **Strong hyperbolic** ($z<-1$): confluent steps in $w=\sqrt{-z}/2$ (kind 4, hyperbolic form of §4.5), with the verified contour as fallback.
- **Very short arcs** ($t<0.1t_p$, $a>0$): confluent steps in $Y=\sqrt y$ (§4.5), with the contour in $Y$ as fallback. With $S=\sinh^2\frac w2=(y_0-Y^2)/(2b)$:
  - $\sinh w=2\sqrt{S(1+S)}$, $\cosh w=1+2S$ and $w=2\operatorname{asinh}\sqrt S$, so $y$ is never formed by cancellation;
  - small $S$ uses series;
  - $z=-16\operatorname{asinh}^2\sqrt{(y_0-Y^2)/2b}$.
- **Parabolic:**
  - $t=t_p$ to within $8\varepsilon t_p$: $z=0$ exactly;
  - within about $10^3\varepsilon$ of $t_p$: $z=(t-t_p)/T'(0)$, exact to $O(z^2)$, where no contour could resolve $T-t$.

### 4.9 Velocity reconstruction

From the half-angle pair at the root, form $\cos v$ and $\sin^2v$. They are recomputed by libm when $\lvert\cos v\rvert<\frac14$, so that near-antipodal transfers keep full relative accuracy. For upper-boundary roots, the small distance $q=4\pi^2-z$ is used directly. Then, with $d=\operatorname{sign}a$ and $\Delta c=d\,c_s-\cos v$ (replaced by $(\sin^2v-s_s^2)/(d\,c_s+\cos v)$ when that form avoids cancellation):

$$
v_{r1}=\frac{\sqrt2\Big(\Delta c+d\,c_s\,\dfrac{\hat r_2-\hat r_1}{\hat r_1(\rho_r+1)}\Big)}{\sqrt y},\qquad
v_{r2}=-\frac{\sqrt2\Big(\Delta c+d\,c_s\,\dfrac{\hat r_1-\hat r_2}{\hat r_2(1/\rho_r+1)}\Big)}{\sqrt y},
$$

$$
v_{t1}=\frac{\sqrt2\,\rho_r\,s_s}{\sqrt y},\qquad v_{t2}=\frac{\sqrt2\,s_s}{\rho_r\sqrt y},\qquad \rho_r=\sqrt{\hat r_2/\hat r_1},
$$

$$
\mathbf v_i=\frac{v_{ri}\,\hat e_i+v_{ti}\,\hat t_i}{\sqrt{s}}\quad(\mu=1\ \text{units}).
$$

This is algebraically the paper's $f,g,\dot g$ reconstruction (Eqs. 6–7), rearranged so that no component is a difference of nearly equal vectors.

### 4.10 Safeguards

These are rarely used, and every one is still a contour method:
- **Interval contours** (`elliptic_root`) on the full zero-revolution elliptic interval, used if the fast path fails.
- **Negative-axis contours** (`hyperbolic_root`) for zero-revolution hyperbolic arcs the fast path declines. In practice these are long-way arcs with $t\lesssim10^{-5}t_p$ ($w\approx30$), where $yF+a$ cancels to about $10^{-10}$ and the confluent acceptance test cannot pass.
- **A certified interval contour** (`safeguarded_z`) as the last resort for every regime: a sign bracket shrunk by contour updates with bisection every third pass. For multiple revolutions it is seeded by the paper's Eq. 24 two-root moments (`pair_moments`) on centered, scaled moments.

If the very-short-hyperbola path (§4.8) declines, the case continues through the normal hyperbolic route above.
- **Optional diagnostics:** the accepted $z$, $q$, the relative residual, and whether a safeguard was used.

### 4.11 Batch solves (`contour_batch`)

A single solve is one long chain of dependent operations: geometry, table, trigonometry, step, trigonometry, velocities. In a loop of single solves, the processor cannot overlap consecutive cases.

`contour_batch` instead processes chunks of 256 cases in stages:
1. geometry, classification and table estimates (including table-based infeasibility);
2. half-angle pairs;
3. confluent steps;
4. velocities;
5. everything else (hyperbolic, near-parabolic, near-$t_{\min}$, declined steps), solved from the already prepared geometry.

Results are identical to `contour_one` bit for bit (`tests/test_contour_batch.f90`, serial and OpenMP). The MATLAB gateway's batch style uses this path.

- **Parallelism.** Chunks run in parallel under OpenMP with dynamic scheduling. Chunk times vary with the case mix, and static scheduling left 8 threads slower than 6 on a CPU with efficiency cores. Dynamic scheduling cut 8-thread time from 37 to 28 ns per case.
- **Small batches.** A batch of at most one chunk (256 cases) skips the OpenMP region entirely. A serialized region costs about 180 ns, and waking a thread team costs 10–50 µs, which earlier dominated batches of 1–100 cases at 2–8 threads.
- **Work arrays.** The chunk work arrays have fixed size: gfortran heap-allocates automatic arrays on every call, about 200 ns.
- **Result.** A batch of one case costs about the same as `contour_one`, and batches of 10 or more are faster than the scalar loop.

### 4.12 Status codes and outputs

| Status | Meaning |
|---:|---|
| 0 | Solution found |
| 1 | Coincident multi-revolution roots at minimum time |
| 5 | Flight time below the multi-revolution minimum |
| −1 | Invalid input |
| −2 | Nearly collinear geometry: transfer plane undefined |
| −4 | No sufficiently accurate root |

`half`: −1 degenerate, 0 ordinary, +1 near-collinear (the transfer plane is poorly conditioned). Inputs are in consistent units with $\mu=1$; for other $\mu$, scale time by $\sqrt\mu$.

---

## 5. Accuracy and performance

**Setup:**
- Apple M1 Pro, macOS, GNU Fortran 16.2, `-O3 -fno-fast-math`, single thread.
- Native benchmark with 2000 random cases per class.
- Truth from closed-form quad-precision Kepler propagation, independent of any Lambert solver.
- Pass criterion: endpoint velocity error $\le10^{-10}(1+\lVert\mathbf v\rVert)$.
- Comparison solver: Russell's ivLam (interpolated initial guess plus Householder iterations), a current production standard.

### 5.1 Time per solve (ns)

| Class | CCS | Russell (ivLam) |
|---|---:|---:|
| Circular, short way | ≈209 | ≈210 |
| Elliptic, short way | ≈218 | ≈215 |
| Elliptic, long way | ≈197 | ≈198 |
| High eccentricity | ≈213 | ≈205 |
| Hyperbolic | ≈241 | ≈205 |
| Near-parabolic | ≈152 | ≈184 |
| Multi-rev, 1 revolution | ≈228 | ≈269 |
| Multi-rev, 2–3 revolutions | ≈235 | ≈281 |
| Multi-rev, 1–10 revolutions | ≈274 | ≈276 |
| Infeasible revolution count (status 5) | ≈53 | ≈57 |
| Strong hyperbola | ≈190 | ≈216 |
| **Overall mean** (3 seeds) | **193–198** | **200–206** |

- In this benchmark CCS averages **0.96–0.99×** Russell's time. It is faster on near-parabolic, strongly hyperbolic and multi-revolution transfers, level on elliptic ones, and about 1.2× on mildly-to-moderately hyperbolic ones.
- On a 30 000-case "hard" mixed dataset exported from the MATLAB generator (`research/speed/ds_bench.f90`), CCS averages **168 ns against Russell's 209**:

  | Category | CCS (ns) | Russell (ns) |
  |---|---:|---:|
  | Elliptic, 0 revolutions | 211 | 202 |
  | Elliptic, 1–4 revolutions | 227–234 | 267–276 |
  | Hard hyperbolic (40% very short arcs) | 244 | 199 |
  | Parabolic | 36 | 171 |

  Through the MATLAB gateway (3009 cases) CCS is about 15% faster than Russell and passes the same cases as the earlier native port of the paper's method, which it replaced.
- The paper reports 2.9 µs (one contour) and 5.9 µs (two contours) on its own hardware in C++. That is not a controlled comparison, but it is roughly an order of magnitude slower.

### 5.2 Accuracy (20 000 cases per class)

| Class | CCS max relative velocity error | CCS failures | Russell max error | Russell failures |
|---|---:|---:|---:|---:|
| Circular | 2.1e−15 | 0 | 4.2e−15 | 0 |
| Elliptic, short way | 4.4e−14 | 0 | 6.6e−13 | 0 |
| Elliptic, long way | 7.4e−14 | 0 | 9.1e−13 | 0 |
| High eccentricity | 2.2e−14 | 0 | 7.3e−12 | 0 |
| Hyperbolic | 6.1e−14 | 0 | 3.6e−13 | 0 |
| Strong hyperbola | 9.5e−15 | 0 | 2.6e−15 | 0 |
| Near-parabolic | 2.6e−15 | 0 | 2.1e−15 | 0 |
| Multi-rev, 1 revolution | 5.1e−11 | 0 | 8.2e−12 | 5 |
| Multi-rev, 2–3 revolutions | 2.3e−12 | 1 | 1.5e−12 | 3 |
| Multi-rev, 1–10 revolutions | 9.5e−12 | 0 | 1.1e−12 | 1 |
| Infeasible (must return status 5) | — | 0 | — | 0 |

Max errors are over passing cases. Totals: CCS 1 failure, Russell 9, out of 220 000 cases (seed 99).

The largest CCS errors are multi-revolution roots very close to $t_{\min}$. There the root itself is ill-conditioned ($dz/dt\propto(t-t_{\min})^{-1/2}$), so this is a property of the problem, not of the method. A stress audit of singular geometries (near-$180^\circ$, near-collinear, extreme radius ratios) fails 11 cases for CCS and 322 for Russell.

**Order 4 and batching** (`benchmark_mixed`, interleaved A/B, M1 Pro):

| | Scalar, order 3 | Scalar, order 4 | Batch, order 4 | Russell |
|---|---:|---:|---:|---:|
| Elliptic, zero revolutions (ns) | 193–205 | 186–196 | 131–155 | 189–215 |
| Multi-revolution (ns) | 224–253 | 203–232 | 133–164 | 253–281 |
| Overall mean (ns) | 191 | 182 | 147 | ≈200 |

Through the MATLAB gateway (3009 hard mixed cases): 0.51 ms scalar and 0.47 ms batch, against 0.63–0.64 ms for Russell.

### 5.3 Where the time goes (one elliptic solve)

| Component | ns |
|---|---:|
| Geometry and velocity reconstruction | ≈22 |
| Table estimate (one cube root, lookup) | ≈25 |
| Confluent step (order 4, one step) | ≈60 |
| Call overhead, dispatch | rest |

- Each step is limited by latency: about six dependent divides and square roots. Folding its two divisions into one gave no gain, because they overlap.
- The remaining contour paths are locating the minimum time near $t_{\min}$ and fallbacks. The hyperbolic gap to Russell comes from the strongly hyperbolic step needing about 2.4 steps from its closed-form estimate; a hyperbolic table was tried and did not pay for its own lookup.

### 5.4 Limitations

- Exactly antipodal ($180^\circ$) transfers have no unique plane; CCS reports status −2, like the paper's formulation.
- Accuracy near $t_{\min}$ is limited by the problem's conditioning.
- The confluent step covers every zero-revolution regime and both multi-revolution branches. Only the minimum-time search near $t_{\min}$ and fallbacks use the verified small contour.
- The tables add about 66 KB to the compiled code (no data file). Classification of infeasible multi-revolution requests from the table is validated empirically (§4.4a), with a margin of about 250 times the measured table error. Requests near $t_{\min}$ always use the exact minimum.

---

## 6. Design alternatives that were tested and rejected

Each of these was tested in 40–50-digit arithmetic and, where relevant, timed in Fortran (`research/REPORT.md`).

| Alternative | Finding |
|---|---|
| The paper's global contours (large circle, then radius 2) | Theorem 1: the error grows like $(r/\lvert\delta\rvert)^N$. Replaced by estimate plus confluent step or small contour |
| Single log-derivative integral $\oint z\,T'/(T-t)$ | Needs $T'$ at every node; does not cancel the pole alias; turns the poles of $T$ into poles of the integrand. Worse at equal cost |
| Pole-deflated residual $R=\sqrt y(y\,n+\sqrt2as^3)-\sqrt2t\,s^3$ with log-derivative moments | Fixes the poles of $T$, but is never more accurate than the reciprocal $1/(T-t)$ at equal node count |
| Unweighted power sums $p_1,p_2$ for the two multi-rev roots | Worse conditioned near $t_{\min}$ than Eq. 24 (the discriminant $2p_2-p_1^2$ cancels). Per-branch solving is used instead |
| Thin ellipse or rectangle (Gauss–Legendre sides) contours | 4–8 digits better than the paper's global circle at equal cost, but only about 1% faster than the small circle as local contours. Not used |
| Contour collapsed onto the real axis (Chebyshev nodes) | Exact for a pure pole, but half the convergence rate per evaluation (one real number per sample) |
| Higher-order generic Taylor kernels; Newton for later steps | Slower than the fixed order-3 closed form |
| Squared residual (removes $\sqrt y$) | Introduces roots from the wrong sign sheet on multi-revolution intervals |
| Larger estimate tables (136 KB, 224 KB, 3.3 MB) | 136 KB: about 4% faster than the 66 KB default; larger: no further gain, since all already make nearly every solve one step |
| Tabulating a correction to the closed-form estimate | The correction varies by O(1) near $\lvert\lambda\rvert\to1$; tabulating the root itself in endpoint-adapted variables needs far fewer points |

---

## 7. Implementation map

All routines are in `fortran/confluentcontour.f90`.

| Routine | Role | Section |
|---|---|---|
| `contour_one`, `contour_batch` | Public entry points (single case, batch) | 4.1 |
| `prepare_geometry` | Normalization, half-angle geometry, plane check | 4.2 |
| `solve_z` | Regime dispatch, residual check, safeguards | 4.1 |
| `elliptic_estimate`, `elliptic_fast` | Zero-rev elliptic estimate and coordinate choice | 4.4 |
| `table_estimate`, `mrev_estimate`, `lb_coord`, `lagrange4` | Table refinement of the estimates; Lancaster–Blanchard coordinate | 4.4a |
| `contour_table_zero_rev.inc`, `contour_table_multi_rev.inc` | Generated tables (`tools/gen_table.py`, `tools/gen_mrev.py`) | 4.4a |
| `hyperbolic_fast`, `hyperbolic_short`, `sqrty_time` | Hyperbolic paths | 4.8 |
| `multirev_fast`, `critical_core` | Feasibility bound, minimum-time function $K$, branch estimates | 4.7 |
| `fast_root` | Root extraction driver: confluent step first, then bracket and verified contour | 4.5, 4.6 |
| `confluent_step`, `confluent_step_h`, `confluent_step_z`, `confluent_step_y` | The confluent contour step: order 4 elliptic; order 3 hyperbolic, near-parabolic, very short hyperbolic | 4.5 |
| `contour_nodes`, `bracket_pole`, `coordinate_rho` | Vectorized node evaluation, the paper's ratio with error estimate, singularity distance | 4.6 |
| `frac_core`, `real_core*`, `rotate_half` | Time-function kernels and exact half-angle rotation | 3.1, 4.5 |
| `reconstruct` | Velocities | 4.9 |
| `elliptic_root`, `hyperbolic_root`, `pair_moments`, `safeguarded_z` | Safeguard contours | 4.10 |

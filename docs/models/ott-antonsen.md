# Ott–Antonsen Dimensionality Reduction

The Ott–Antonsen (OA) ansatz collapses the infinite-dimensional continuum
Kuramoto model onto a finite set of ODEs for the order parameter(s). This page
gives the derivation; the implementation notes live in
[`src/MM/models/OA-Ansatz.md`](../../src/MM/models/OA-Ansatz.md).

---

## 1. Mean field and the continuity equation

For all-to-all coupling, the Kuramoto model reduces (see
[kuramoto.md](kuramoto.md)) to the mean-field form

$$
\frac{\mathrm{d}\theta}{\mathrm{d}t}
= \omega + K\,\mathrm{Im}\!\left[r(t)\,e^{-i\theta}\right],
\qquad
r(t) = \rho\,e^{i\varphi}.
$$

In the limit $N\to\infty$ the state is a density $f(\theta,\omega;t)$, with
$\int f\,\mathrm{d}\theta = g(\omega)$, obeying the continuity equation

$$
\frac{\partial f}{\partial t}
+ \frac{\partial}{\partial\theta}\left(f\,\dot\theta\right) = 0.
$$

Write $\mathrm{Im}[r e^{-i\theta}] = \frac{1}{2i}(r e^{-i\theta} - r^* e^{i\theta})$
and Fourier-expand the density,

$$
f(\theta,\omega;t)
= \frac{g(\omega)}{2\pi}\sum_{n=-\infty}^{\infty} a_n(\omega;t)\,e^{in\theta},
\qquad
a_0 = 1,\qquad a_{-n} = a_n^{*}.
$$

---

## 2. The Fourier hierarchy

Substituting the expansion into the continuity equation and collecting the
coefficient of $e^{in\theta}$ gives the infinite hierarchy

$$
\frac{\mathrm{d}a_n}{\mathrm{d}t}
= -i\,n\,\omega\,a_n
+ \frac{K\,n}{2}\left(r^{*}\,a_{n-1} - r\,a_{n+1}\right),
\qquad n\ge 1 .
$$

(The factor $n$ comes from the derivative $\partial_\theta$; the $K/2$ from the
identity $\sin(\theta_k-\theta_j)=\tfrac1{2i}(e^{i(\theta_k-\theta_j)}-e^{-i(\theta_k-\theta_j)})$.)
The order parameter is the first moment,

$$
r(t) = \iint e^{i\theta}\,f\,\mathrm{d}\theta\,\mathrm{d}\omega
      = \int g(\omega)\,a_{-1}(\omega;t)\,\mathrm{d}\omega
      = \int g(\omega)\,a_1^{*}(\omega;t)\,\mathrm{d}\omega .
$$

---

## 3. The OA ansatz

Ott & Antonsen (2008) observed that the *low-dimensional manifold*

$$
a_n(\omega;t) = \left[a(\omega;t)\right]^{n},
\qquad n\ge 1,
$$

is **invariant** under the hierarchy — if it holds at one instant it holds
forever. (The restriction to this manifold is exact, not an approximation.) With
$a_n = a^{n}$ the whole hierarchy collapses to the single equation

$$
\frac{\mathrm{d}a}{\mathrm{d}t}
= -i\,\omega\,a + \frac{K}{2}\left(r^{*} - r\,a^{2}\right),
\qquad
r(t) = \int g(\omega)\,a^{*}(\omega;t)\,\mathrm{d}\omega .
$$

---

## 4. Lorentzian frequencies: analytic continuation

For a Lorentzian (Cauchy) distribution

$$
g(\omega) = \frac{\gamma}{\pi}\,\frac{1}{(\omega-\mu)^{2}+\gamma^{2}},
$$

$g$ has simple poles at $\omega = \mu \pm i\gamma$. The integral defining $r$ can
be evaluated by contour integration: since $a(\omega;t)$ is analytic in the upper
half-plane and decays, closing the contour in the lower half-plane picks up the
single pole at $\omega = \mu - i\gamma$, giving

$$
r(t) = a^{*}(\mu - i\gamma;\,t).
$$

Evaluating the $a$-equation at this pole (with $\dot r = \dot a^{*}(\mu-i\gamma)$)
yields the **single-community OA equation**

$$
\frac{\mathrm{d}r}{\mathrm{d}t}
= (-\gamma + i\mu)\,r + \frac{K}{2}\left(1 - |r|^{2}\right) r,
$$

which is the complex ODE implemented by `OA`. In Cartesian form $r = x + iy$
($|r|^2 = x^2+y^2$):

$$
\frac{\mathrm{d}x}{\mathrm{d}t} = -\gamma x - \mu y + \frac{K}{2}(1-x^2-y^2)\,x,
\qquad
\frac{\mathrm{d}y}{\mathrm{d}t} = \mu x - \gamma y + \frac{K}{2}(1-x^2-y^2)\,y.
$$

The Cartesian form is preferred because the polar phase equation carries a
$1/\rho$ singularity as $\rho\to 0$ (the phase is undefined at the origin).

---

## 5. Multiple communities

Split the population into $C$ communities with fractions $\eta_c=N_c/N$,
Lorentzian parameters $(\gamma_c,\mu_c)$, order parameters $r_c$, and coupling
matrix $K_{cc'}$. Community $c$ sees the mean field
$\sum_{c'} K_{cc'}\,\eta_{c'}\,r_{c'}$, and the same ansatz applied per
community gives the coupled system

$$
\frac{\mathrm{d}r_c}{\mathrm{d}t}
= (-\gamma_c + i\mu_c)\,r_c
- \frac{1}{2}\sum_{c'} K_{cc'}\,\eta_{c'}\left(r_c^2\,\overline{r_{c'}} - r_{c'}\right),
$$

implemented by `OAGeneral`. Note the factor $K_{cc'}\,\eta_{c'}$ (coupling times
the *source* community fraction) and the $1/2$. The global order parameter is the
exact $\eta$-weighted sum

$$
R(t) = \sum_{c=1}^{C}\eta_c\,r_c(t).
$$

---

## API

All in `namespace MathEngine` (`src/MM/models/OA-Ansatz.hpp`):

| Function | Purpose |
|---|---|
| `OA(time, state, dstate, gamma, mu, K)` | single community, state `[x, y]` |
| `OAGeneral(time, state, dstate, gammas, mus, eta, K)` | multi-community, interleaved `[x_0,y_0,x_1,y_1,...]` |
| `OA_wrapper(OAParams)` / `OAGeneral_wrapper(OAGeneralParams)` | `MyFunc` for the solvers |

## Use cases

- Cheap bifurcation studies of the synchronisation transition in the
  thermodynamic limit (the Lorentzian OA result reproduces
  $\rho = \sqrt{1-2\gamma/K}$ exactly).
- Large communities: the reduction is $O(C)$ in communities rather than $O(N)$
  in oscillators.
- Compare finite-$N$ Kuramoto simulations against the OA prediction.

## References

- Ott & Antonsen (2008), *Chaos* **18**, 037113.
- Ott & Antonsen (2009), *Chaos* **19**, 023117.

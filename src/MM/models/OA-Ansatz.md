# Ott–Antonsen Ansatz

The Ott–Antonsen (OA) ansatz collapses the infinite-dimensional continuum
Kuramoto model onto a finite set of ODEs for the order parameter(s). This
document records the equations implemented in `OA-Ansatz.hpp`, together with the
derivation and the notational conventions.

---

## 1. All-to-all Kuramoto model

The (general) Kuramoto model with all-to-all coupling is

$$
\frac{\mathrm{d}\theta_j}{\mathrm{d}t}
= \omega_j + \frac{K}{N}\sum_{k=1}^{N}\sin(\theta_k - \theta_j),
\qquad j = 1,\dots,N.
$$

The order parameter is the first moment of the phase distribution,

$$
r(t) = \rho(t)\,e^{i\varphi(t)} = \frac{1}{N}\sum_{j=1}^{N}e^{i\theta_j(t)}.
$$

---

## 2. Continuum limit

In the limit $N\to\infty$ the phase distribution is described by a density
$f(\theta,\omega;t)$ obeying the continuity equation

$$
\frac{\partial f}{\partial t} + \frac{\partial}{\partial \theta}(f\,\dot\theta) = 0,
\qquad
\dot\theta = \omega + K\,\mathrm{Im}\!\big[e^{-i\theta} r(t)\big].
$$

Assume the intrinsic frequencies are drawn from a Lorentzian

$$
g(\omega) = \frac{\gamma}{\pi}\,\frac{1}{(\omega-\mu)^{2}+\gamma^{2}},
$$

with centre $\mu$ (mean frequency) and half-width $\gamma$.

---

## 3. The OA ansatz (single community)

Take the Fourier expansion

$$
f(\theta,\omega;t) = \frac{g(\omega)}{2\pi}\sum_{n=-\infty}^{\infty} a_n(\omega;t)\,e^{in\theta}
$$

and impose the OA closure $a_n = a^{n}$ for $n>0$ (with $a_{-n}=a_n^{*}$).
The order parameter is then recovered by analytic continuation of $a(\omega;t)$
to the pole of $g(\omega)$ in the appropriate half-plane,
$r(t) = a^{*}(\mu - i\gamma; t)$. The reduced system is the single complex ODE

$$
\frac{\mathrm{d}r}{\mathrm{d}t}
= (-\gamma + i\mu)\,r \;-\; \frac{K}{2}\left(r^{2}r^{*} - r\right),
$$

or, equivalently,

$$
\frac{\mathrm{d}r}{\mathrm{d}t}
= (-\gamma + i\mu)\,r \;+\; \frac{K}{2}\left(1 - |r|^{2}\right) r.
$$

> **Note the factor $K/2$.** It follows from writing
> $\sin(\theta_k-\theta_j) = \frac{1}{2i}(e^{i(\theta_k-\theta_j)} - e^{-i(\theta_k-\theta_j)})$.

### Cartesian form

Writing $r = x + iy$ (so $|r|^{2} = x^{2}+y^{2}$) gives the two real equations
implemented by `OA`:

$$
\begin{aligned}
\frac{\mathrm{d}x}{\mathrm{d}t}
&= -\gamma x - \mu y + \frac{K}{2}\left(1 - x^{2} - y^{2}\right)x,\\[0.5em]
\frac{\mathrm{d}y}{\mathrm{d}t}
&= \mu x - \gamma y + \frac{K}{2}\left(1 - x^{2} - y^{2}\right)y.
\end{aligned}
$$

The Cartesian form is preferred because the polar phase equation carries a
$1/\rho$ singularity as $\rho\to 0$ (the phase is undefined at the origin).

---

## 4. Multiple communities

Split the population into $C$ communities. Community $c$ holds a fraction

$$
\eta_c = \frac{N_c}{N}, \qquad \sum_{c=1}^{C}\eta_c = 1,
$$

of the oscillators and has its own Lorentzian parameters $\gamma_c$, $\mu_c$ and
order parameter

$$
r_c(t) = \frac{1}{N_c}\sum_{j\in c} e^{i\theta_j(t)}.
$$

With coupling strengths $K_{c,c'}$ between communities, the mean field seen by
community $c$ is $\sum_{c'} K_{c,c'}\,\eta_{c'}\,r_{c'}$, and the OA reduction
yields the coupled system

$$
\frac{\mathrm{d}r_c}{\mathrm{d}t}
= (-\gamma_c + i\mu_c)\,r_c
\;-\; \frac{1}{2}\sum_{c'=1}^{C}
   K_{c,c'}\,\eta_{c'}\left(r_c^{2}\,\overline{r_{c'}} - r_{c'}\right).
$$

This is `OAGeneral`. Note the combination $K_{c,c'}\,\eta_{c'}$ (coupling times
the *source* community fraction) together with the $1/2$. The global order
parameter is the $\eta$-weighted sum

$$
R(t) = \sum_{c=1}^{C}\eta_c\,r_c(t),
$$

which is an exact linear identity (no further ansatz is required to aggregate).

---

## 5. API reference

All symbols live in `namespace MathEngine`.

| Function | Purpose |
|---|---|
| `OA(time, state, dstate, gamma, mu, K)` | Single community, Cartesian state `[x, y]`. |
| `OAGeneral(time, state, dstate, gammas, mus, eta, K)` | Multi-community, state `[x_0, y_0, x_1, y_1, ...]`. |
| `OAOrder(state, eta)` | Global order parameter, returns `[Re, Im, rho]`. |
| `OAOrderPerCommunity(state)` | Per-community magnitudes `rho_c`, size `C`. |
| `OAOrder(results, eta)` | Global order-parameter time series (cols: `time, Re, Im, rho`). |
| `OA_wrapper(OAParams)` | `MyFunc` for the single-community model. |
| `OAGeneral_wrapper(OAGeneralParams)` | `MyFunc` for the multi-community model. |

Parameter structs:

```cpp
struct OAParams       { double gamma; double mu; double K; };
struct OAGeneralParams{ dVec gammas; dVec mus; dVec eta; dMatrix K; int C; };
```

For the single-community model pass `eta = {1.0}` to the order-parameter helpers.

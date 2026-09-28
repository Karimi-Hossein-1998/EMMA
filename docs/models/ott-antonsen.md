# Ott–Antonsen Dimensionality Reduction

The Ott–Antonsen (OA) ansatz collapses the infinite-dimensional continuum
Kuramoto model onto a finite set of ODEs for the order parameter(s). This page
summarises the model and its use; a full derivation is kept alongside the code in
[`src/MM/models/OA-Ansatz.md`](../../src/MM/models/OA-Ansatz.md).

## The reduction

For all-to-all coupling with intrinsic frequencies drawn from a Lorentzian
$g(\omega)=\dfrac{\gamma}{\pi}\dfrac{1}{(\omega-\mu)^2+\gamma^2}$, the OA ansatz
closes the Fourier hierarchy of the phase distribution, giving a single complex
ODE for the order parameter $r = x + iy$:

$$
\frac{\mathrm{d}r}{\mathrm{d}t}
= (-\gamma + i\mu)\,r + \frac{K}{2}\left(1 - |r|^2\right) r,
$$

i.e. the two real equations implemented by `OA`:

$$
\frac{\mathrm{d}x}{\mathrm{d}t} = -\gamma x - \mu y + \frac{K}{2}(1-x^2-y^2)\,x,
\qquad
\frac{\mathrm{d}y}{\mathrm{d}t} = \mu x - \gamma y + \frac{K}{2}(1-x^2-y^2)\,y.
$$

The Cartesian form is preferred because the polar phase equation diverges as
$1/\rho$ when $\rho\to 0$ (the phase is undefined at the origin).

## Multiple communities

Splitting the population into $C$ communities with fractions $\eta_c$,
Lorentzian parameters $(\gamma_c,\mu_c)$, and coupling matrix $K_{cc'}$ gives

$$
\frac{\mathrm{d}r_c}{\mathrm{d}t}
= (-\gamma_c + i\mu_c)\,r_c
- \frac{1}{2}\sum_{c'} K_{cc'}\,\eta_{c'}\left(r_c^2\,\overline{r_{c'}} - r_{c'}\right).
$$

The global order parameter is the exact $\eta$-weighted sum
$R = \sum_c \eta_c r_c$.

## API

All in `namespace MathEngine` (`src/MM/models/OA-Ansatz.hpp`):

| Function | Purpose |
|---|---|
| `OA(time, state, dstate, gamma, mu, K)` | single community, state `[x, y]` |
| `OAGeneral(time, state, dstate, gammas, mus, eta, K)` | multi-community, interleaved `[x_0,y_0,x_1,y_1,...]` |
| `OA_wrapper(OAParams)` / `OAGeneral_wrapper(OAGeneralParams)` | `MyFunc` for the solvers |

## Use cases

- Cheap bifurcation studies of the synchronisation transition in the
  thermodynamic limit.
- Large communities (the reduction is $O(C)$ in the number of communities rather
  than $O(N)$ in oscillators).
- Compare finite-$N$ Kuramoto simulations against the OA prediction.

## References

- Ott & Antonsen (2008), *Chaos* **18**, 037113.
- Ott & Antonsen (2009), *Chaos* **19**, 023117.

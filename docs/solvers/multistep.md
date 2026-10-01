# Multistep Methods: Adams–Bashforth and Adams–Bashforth–Moulton

Linear multistep methods reuse several past values of $f$ (rather than
re-evaluating inside one step), which makes them cheap per step. The toolkit
provides the Adams–Bashforth predictor and the Adams–Bashforth–Moulton
predictor–corrector, of orders 1–10.

## Adams–Bashforth (explicit predictor)

Using the past $k$ derivatives $f_n, f_{n-1}, \dots, f_{n-k+1}$:

$$
y_{n+1} = y_n + \mathrm{d}t\sum_{j=0}^{k-1} \beta_j\, f_{n-j},
$$

where the $\beta_j$ are the Adams–Bashforth coefficients (order $p=k$). For
example, order 4:

$$
y_{n+1} = y_n + \frac{\mathrm{d}t}{24}\left(55 f_n - 59 f_{n-1} + 37 f_{n-2} - 9 f_{n-3}\right).
$$

The coefficients are precomputed up to order 10 in
`src/MM/miscellany/ABM-Coefs/` and exposed via `abm_coefs::ab_coefs`.

### Derivation (Newton backward-difference interpolation)

The Adams–Bashforth coefficients follow by integrating the polynomial that
interpolates $f$ through the past $k$ points. Let
$s=(t-t_n)/\mathrm{d}t$ and expand $f$ in Newton backward-difference form,

$$
f(t_n + s\,\mathrm{d}t)
= f_n + s\,\nabla f_n + \frac{s(s+1)}{2!}\,\nabla^{2} f_n
+ \frac{s(s+1)(s+2)}{3!}\,\nabla^{3} f_n + \cdots,
$$

where $\nabla$ is the backward-difference operator
($\nabla f_n = f_n - f_{n-1}$, $\nabla^{2} f_n = f_n - 2 f_{n-1} + f_{n-2}$, …).
Integrating term by term over $s\in[0,1]$,

$$
\int_0^1 1\,\mathrm{d}s = 1,\quad
\int_0^1 s\,\mathrm{d}s = \frac12,\quad
\int_0^1 \frac{s(s+1)}{2}\,\mathrm{d}s = \frac5{12},\quad
\int_0^1 \frac{s(s+1)(s+2)}{6}\,\mathrm{d}s = \frac38,
$$

gives

$$
y_{n+1} = y_n + \mathrm{d}t\left[
f_n + \tfrac12\nabla f_n + \tfrac5{12}\nabla^{2} f_n
+ \tfrac38\nabla^{3} f_n + \cdots\right].
$$

Expanding the differences and collecting coefficients of
$f_n, f_{n-1}, f_{n-2}, f_{n-3}$ reproduces the order-4 table:

$$
y_{n+1} = y_n + \frac{\mathrm{d}t}{24}
\left(55 f_n - 59 f_{n-1} + 37 f_{n-2} - 9 f_{n-3}\right).
$$

The Adams–Moulton coefficients are derived identically, except the interpolating
polynomial is taken through $f_{n+1}$ *and* the past points, so the quadrature
starts at $s=1$ (extending back to $s=0$).

## Adams–Moulton (implicit corrector)

The corrector uses $f_{n+1}$ as well:

$$
y_{n+1} = y_n + \mathrm{d}t\sum_{j=-1}^{k-1} \beta_j^{*}\, f_{n+1-j},
$$

and has order $p=k+1$ (one higher than the explicit predictor of the same number
of stages). It is implicit, so in practice it is applied as a **corrector** to a
Bashforth prediction:

1. **Predict** $y_{n+1}^{(0)}$ with Adams–Bashforth.
2. **Evaluate** $f_{n+1}^{(0)} = f(t_{n+1}, y_{n+1}^{(0)})$.
3. **Correct** with the Moulton formula (optionally iterated `iterations` times,
   a P(EC)$^k$ scheme).

## Bootstrapping

Multistep methods need $k$ starting values. The implementation bootstraps the
first `order-1` steps with RK4 (which is itself documented in
[runge-kutta.md](runge-kutta.md)).

## Stability and error

Adams methods are only weakly stable for $k>2$: the linear test equation gives
stability regions much smaller than RK4 of comparable order. They are best suited
to non-stiff problems where the derivative is expensive to evaluate (one $f$
evaluation per step for AB, one or two for ABM), such as very large coupled ODE
systems.

- Adams–Bashforth order $k$: local error $\mathcal{O}(\mathrm{d}t^{k+1})$.
- Adams–Moulton order $k+1$: local error $\mathcal{O}(\mathrm{d}t^{k+2})$.

## Usage

```cpp
solverParams.solverFunc = adams_bashforth_wrapper();          // AB
solverParams.solverFunc = adams_bashforth_moulton_wrapper();  // ABM (P(EC)^k)
```

Selected in the GUI as `AB` and `ABM`; the Solver tab exposes the order (1–10)
and, for ABM, the number of corrector iterations.

## References

- Hairer, Nørsett & Wanner, *Solving Ordinary Differential Equations I* (1993).
- Butcher, *Numerical Methods for Ordinary Differential Equations* (2016).

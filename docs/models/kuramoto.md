# Kuramoto Model

The Kuramoto model is the canonical description of synchronisation in coupled
phase oscillators. It models, for instance, networks of fireflies, Josephson
junctions, neuronal populations, and power-grid generators. The MM toolkit ships
three variants (general dense, sparse, and modular) plus order-parameter helpers.

## Model

The general Kuramoto model with phase-lag is

$$
\frac{\mathrm{d}\theta_i}{\mathrm{d}t}
= \omega_i + \frac{K}{N}\sum_{j=1}^{N} A_{ij}\,\sin(\theta_j - \theta_i - \alpha),
\qquad i = 1,\dots,N,
$$

where $\theta_i$ is the phase of oscillator $i$, $\omega_i$ its natural frequency,
$K$ the coupling strength, $A_{ij}$ the adjacency (coupling) weight, and $\alpha$
a phase-lag (frustration) parameter. With $A_{ij}=1$ and $\alpha=0$ this reduces
to the classical all-to-all Kuramoto model.

The coupling is **normalised by** $N$ so the mean field is well-defined in the
thermodynamic limit $N\to\infty$.

---

## 1. Mean-field reduction

For all-to-all coupling the sum can be rewritten in terms of the complex order
parameter

$$
r(t) = \rho(t)\,e^{i\varphi(t)} = \frac{1}{N}\sum_{j=1}^{N} e^{i\theta_j(t)}.
$$

Multiplying by $e^{-i\theta_i}$ and taking the imaginary part,

$$
\mathrm{Im}\!\left[r\,e^{-i\theta_i}\right]
= \frac{1}{N}\sum_{j}\sin(\theta_j - \theta_i),
$$

so the coupling term is exactly

$$
\frac{K}{N}\sum_{j}\sin(\theta_j - \theta_i)
= K\,\rho\,\sin(\varphi - \theta_i),
\qquad \rho = |r|.
$$

The all-to-all model therefore reduces to the **mean-field equation**

$$
\frac{\mathrm{d}\theta_i}{\mathrm{d}t}
= \omega_i + K\rho\,\sin(\varphi - \theta_i).
$$

Each oscillator now only "sees" the collective variables $(\rho, \varphi)$ — the
coupling is mediated by the order parameter, not by pairwise interactions. Here
$\rho\in[0,1]$ quantifies phase coherence: $\rho\approx 0$ is incoherence and
$\rho\approx 1$ full synchrony. This is the quantity plotted live by the GUI.

---

## 2. Continuum limit

As $N\to\infty$ the ensemble is described by a phase density
$f(\theta,\omega;t)$, normalised so that $\int f\,\mathrm{d}\theta = g(\omega)$,
where $g(\omega)$ is the intrinsic-frequency distribution. Because phases are
conserved (oscillators are neither created nor destroyed), $f$ obeys the
continuity equation

$$
\frac{\partial f}{\partial t}
+ \frac{\partial}{\partial \theta}\left(f\,\dot\theta\right) = 0,
\qquad
\dot\theta = \omega + K\rho\,\sin(\varphi - \theta)
         = \omega + K\,\mathrm{Im}\!\left[r\,e^{-i\theta}\right].
$$

The order parameter is the first Fourier moment of the density,

$$
r(t) = \iint e^{i\theta}\,f(\theta,\omega;t)\,\mathrm{d}\theta\,\mathrm{d}\omega.
$$

---

## 3. Self-consistency equation and the synchronisation threshold

In the stationary (synchronised) state the population rotates rigidly with a
mean frequency $\Omega$ (for a symmetric $g$, $\Omega$ is the mean of $g$). Pass
to the rotating frame $\psi_i = \theta_i - \Omega t$ and, by rotational
invariance, set $\varphi=0$. The **locked** oscillators satisfy
$\dot\psi=0$, i.e.

$$
\omega_i - \Omega + K\rho\sin(\psi_i) = 0
\quad\Longrightarrow\quad
\sin\psi_i = \frac{\omega_i - \Omega}{K\rho},
$$

which is solvable only when $|\omega_i-\Omega|\le K\rho$. The **drifting**
oscillators ($|\omega_i-\Omega|>K\rho$) sweep the circle and contribute nothing
to $r$ on average. The order parameter is thus built entirely from locked
oscillators:

$$
\rho = \int_{|\omega-\Omega|\le K\rho} \cos\psi\,g(\omega)\,\mathrm{d}\omega,
\qquad \cos\psi = \sqrt{1-\left(\frac{\omega-\Omega}{K\rho}\right)^{2}}.
$$

Substituting $\omega = \Omega + K\rho\,x$ gives the **self-consistency equation**

$$
\rho = K\rho\int_{-1}^{1} \sqrt{1-x^{2}}\; g(\Omega + K\rho\,x)\,\mathrm{d}x,
$$

which, for $\rho>0$, becomes

$$
1 = K\int_{-1}^{1} \sqrt{1-x^{2}}\; g(\Omega + K\rho\,x)\,\mathrm{d}x.
$$

Taking the limit $\rho\to 0^{+}$ (the onset of synchronisation) and using
$\int_{-1}^{1}\sqrt{1-x^{2}}\,\mathrm{d}x = \pi/2$ yields the celebrated
**critical coupling**

$$
\boxed{\;K_c = \frac{2}{\pi\,g(\Omega)}\;}
\qquad\text{(and } K_c = \frac{2}{\pi\,g(0)} \text{ for a symmetric } g\text{).}
$$

For a Lorentzian $g(\omega)=\dfrac{\gamma}{\pi}\dfrac{1}{\omega^{2}+\gamma^{2}}$
(centred at $0$), $g(0)=1/(\pi\gamma)$, so

$$
K_c = 2\gamma.
$$

---

## 4. Bifurcation and the order-parameter scaling

Expanding the self-consistency equation in $\rho$ (for symmetric $g$, so
$g'(0)=0$) gives

$$
1 = K\left[\frac{\pi}{2}\,g(0) + \frac{\pi}{16}\,g''(0)\,K^{2}\rho^{2}
      + \mathcal{O}(\rho^{4})\right],
$$

where $\int_{-1}^{1}\sqrt{1-x^{2}}\,x^{2}\,\mathrm{d}x = \pi/8$. Using
$K_c=2/(\pi g(0))$,

$$
\rho^{2} = -\frac{16\,(K-K_c)}{\pi K_c^{4}\,g''(0)}
+ \mathcal{O}\!\left((K-K_c)^{2}\right).
$$

Since $g''(0)<0$ for a unimodal distribution, $\rho^{2}>0$ for $K>K_c$ and

$$
\rho \;\propto\; \sqrt{K - K_c},
$$

a **supercritical pitchfork bifurcation** at $K_c$: the incoherent state
($\rho=0$) is stable for $K<K_c$ and loses stability for $K>K_c$, where a
synchronised branch appears. For the Lorentzian distribution the Ott–Antonsen
reduction gives the *exact* result

$$
\rho = \sqrt{1 - \frac{K_c}{K}} = \sqrt{1 - \frac{2\gamma}{K}},
\qquad K > K_c,
$$

which reproduces $\rho\propto\sqrt{K-K_c}$ near threshold. (See
[ott-antonsen.md](ott-antonsen.md).)

---

## Variants

| Variant | Header | Coupling |
|---|---|---|
| General (dense) | `models/kuramoto/general.hpp` | dense $A_{ij}$ |
| Sparse | `models/kuramoto/sparse.hpp` | `SparsedMatrix` adjacency (only non-zero edges) |
| Modular | `models/kuramoto/special.hpp` | intra-module $K$, inter-module $Q$ (no explicit matrix) |

The modular variant groups oscillators into blocks of size `module_size` and uses

$$
\frac{\mathrm{d}\theta_i}{\mathrm{d}t}
= \omega_i + \frac{1}{N}\sum_{j\neq i} K_{m_i m_j}\,\sin(\theta_j-\theta_i-\alpha),
\qquad
K_{mm'}=
\begin{cases}
K_{\mathrm{intra}} & m=m',\\
K_{\mathrm{inter}} & m\neq m'.
\end{cases}
$$

---

## Order-parameter helpers

`models/kuramoto/order.hpp` provides `calculate_order` (global, per-node, and
per-module), returning triples $[\langle\sin\theta\rangle,\langle\cos\theta\rangle,\rho]$:

$$
\rho = \sqrt{\langle\sin\theta\rangle^2 + \langle\cos\theta\rangle^2},
\qquad
\langle f(\theta)\rangle = \frac{1}{N}\sum_{j} f(\theta_j).
$$

Local (weighted) variants use the adjacency as the weight.

## Wrapper pattern

Each variant exposes a `*_wrapper` that captures parameters by value and returns a
`MyFunc` for the ODE solvers:

```cpp
KuramotoParams p;                 // { K, N, omega, adj, alpha }
p.omega = ...; p.adj = ...;
auto rhs = MathEngine::kuramoto_general_wrapper(p);
solverParams.solverParams.derivative = rhs;
```

## Use cases

- Study the onset of synchronisation as $K$ crosses $K_c=2/(\pi g(0))$
  ($K_c=2\gamma$ for Lorentzian frequencies).
- Explore chimera states and partial synchrony with modular or small-world
  topologies (see [network-topology.md](../network-topology.md)).
- Phase-lag $\alpha$ introduces frustration and can stabilise splay or
  asynchronous states.

## References

- Kuramoto (1975), *Lecture Notes in Physics* **39**; also in Araki (ed.),
  *International Symposium on Mathematical Problems in Theoretical Physics*.
- Strogatz (2000), *Physica D* **143**, 1–20.
- Acebrón, Bonilla, Pérez Vicente, Ritort & Spigler (2005), *Rev. Mod. Phys.* **77**, 137.

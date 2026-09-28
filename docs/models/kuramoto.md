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
a phase-lag (frustration) parameter. With $A_{ij}=1$ and $\alpha=0$ this reduces to
the classical all-to-all Kuramoto model.

The coupling is **normalised by** $N$ so the mean field is well-defined in the
thermodynamic limit $N\to\infty$.

### Mean field and the order parameter

Writing the complex order parameter

$$
r(t) = \rho(t)\,e^{i\varphi(t)} = \frac{1}{N}\sum_{j=1}^{N} e^{i\theta_j(t)},
$$

the all-to-all model can be recast exactly as a mean-field equation:

$$
\frac{\mathrm{d}\theta_i}{\mathrm{d}t} = \omega_i + K\rho\,\sin(\varphi - \theta_i).
$$

$\rho\in[0,1]$ quantifies phase coherence: $\rho\approx 0$ is incoherence and
$\rho\approx 1$ full synchrony. This is the quantity plotted live by the GUI.

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

- Study the onset of synchronisation as $K$ crosses the critical value
  $K_c = 2/(\pi g(0))$ for a Lorentzian frequency distribution $g(\omega)$.
- Explore chimera states and partial synchrony with modular or small-world
  topologies (see [network-topology.md](../network-topology.md)).
- Phase-lag $\alpha$ introduces frustration and can stabilise splay or
  asynchronous states.

## References

- Kuramoto (1975), *Lecture Notes in Physics* **39**.
- Strogatz (2000), *Physica D* **143**, 1–20.
- Acebrón et al. (2005), *Rev. Mod. Phys.* **77**, 137.

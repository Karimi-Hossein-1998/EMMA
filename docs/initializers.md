# Initial Conditions & Distributions

`src/MM/initializers/initials.hpp` provides seeded random samplers for the
initial phases / frequencies of oscillator models. They are exposed in the GUI's
Initial Conditions tab.

## Distributions

All generators take a seed for reproducibility.

| Function | Distribution | Density / support |
|---|---|---|
| `random_uniform(N, a, b, seed)` | uniform | $[a,b)$ |
| `random_normal(N, μ, σ, seed)` | Gaussian | $\mathcal{N}(\mu,\sigma^2)$ |
| `random_cauchy(N, x₀, γ, seed)` | Cauchy/Lorentzian | $\dfrac{\gamma}{\pi}\dfrac{1}{(x-x_0)^2+\gamma^2}$ |
| `random_exponential(N, λ, seed)` | exponential | $\lambda e^{-\lambda x}$ |
| `random_circle(N, seed)` | uniform on the circle | $[-\pi,\pi)$ |

The Cauchy (Lorentzian) sampler is the natural frequency distribution that makes
the Ott–Antonsen ansatz exact (see [ott-antonsen.md](models/ott-antonsen.md)).

## Structured states

- **Splay** — $N$ phases equally spaced on the circle:
  $\theta_i = \dfrac{2\pi i}{N}$. This is a fixed point of the unperturbed
  all-to-all Kuramoto model and minimises the order parameter.
- **Splay perturbed** — splay plus uniform noise in $[-\text{amp},\text{amp}]$.
- **Modules** — generates per-module values (optionally identical across
  modules), used for modular Kuramoto topologies.

## API

```cpp
enum class InitState { Uniform, Normal, Cauchy, Exponential, Circle, Splay, SplayPerturbed, Modules };
enum class InitType  { Uniform, Normal, Cauchy, Exponential, Circle, Splay, SplayPerturbed };

template <FPNumber Num> struct InitializerParams { InitState initState; InitType moduleType;
    Num param1, param2; size_t N, numModules, moduleSize, seed; bool identical; };

template <FPNumber Num> Vec<Num> initialize_vector(const InitializerParams<Num>&);
```

## Use cases

- Uniform phases model a fully incoherent initial state ($\rho\approx 0$).
- Splay states test the stability of the incoherent fixed point.
- Lorentzian frequencies reproduce the textbook Kuramoto synchronisation
  threshold $K_c = 2/(\pi g(0))$.

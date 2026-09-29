# Models

The `src/MM/models` module hosts the mathematical models. Each model exposes a
`*_wrapper` (or, for molecular dynamics, an engine class) that plugs into the
app's run path.

| Model | Type | Document |
|---|---|---|
| Kuramoto | Phase oscillators (ODE) | [kuramoto.md](kuramoto.md) |
| Ott–Antonsen | Dimensionality-reduced phase oscillators (ODE) | [ott-antonsen.md](ott-antonsen.md) |
| Molecular dynamics | 2D classical N-body (symplectic integrator) | [molecular-dynamics.md](molecular-dynamics.md) |
| Random walk | Stochastic lattice/continuous walk (SoA engine) | [random-walk.md](random-walk.md) |

Two integration styles are used:

- **ODE models** (Kuramoto, Ott–Antonsen) produce a `MyFunc` derivative fed to the
  Runge–Kutta / Adams solvers (see [../solvers/](../solvers/)).
- **Molecular dynamics** is integrated with symplectic velocity-Verlet (not an
  RK/AB method) through a dedicated engine — see its document for the rationale.
- **Random walk** is integrated with its own stochastic `Step()` (also not an
  RK/AB method) through a dedicated engine.

See also the umbrella header `src/MM/models/models.phase-oscillators.hpp` and
`src/MM/models/molecular-dynamics.hpp`.

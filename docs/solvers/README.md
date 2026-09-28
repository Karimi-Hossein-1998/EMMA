# ODE Solvers

The `src/MM/solvers` module provides explicit numerical integrators for the
initial-value problem $\dot y = f(t,y)$, all exposing the `SolverFunc` interface.

| Family | Methods | Document |
|---|---|---|
| Runge–Kutta | RK1 (Euler), RK2 (midpoint), RK3, RK4, RK4 3/8, RK4 Gill, RK4 Ralston | [runge-kutta.md](runge-kutta.md) |
| Multistep | Adams–Bashforth (predictor), Adams–Bashforth–Moulton (predictor–corrector) | [multistep.md](multistep.md) |

Implementation layout:

```
solvers/ODE/rk/explicit/    rk1-solver.hpp, rk2-solver.hpp, rk3-solver.hpp,
                            rk4-solver.hpp, rk4-variants.hpp
solvers/ODE/multistep/      ab-solver.hpp, abm-solver.hpp
```

All solvers support an `onStep` callback that the GUI uses for live plotting
(the callback receives a `OneStepSolverResult` each step).

> **Note** — the molecular-dynamics model is **not** integrated with these
> solvers: it uses a symplectic velocity-Verlet integrator (see
> [../models/molecular-dynamics.md](../models/molecular-dynamics.md)).

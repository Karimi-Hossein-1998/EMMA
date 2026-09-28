# EMMA — Mathematical Modelling Toolkit Documentation

EMMA is a graphical front-end for the header-only **Mathematical Modelling
(MM)** toolkit living in `src/MM`. The toolkit provides numerical models
(phase oscillators and a 2D molecular-dynamics engine), ODE solvers, network
generators, initial-condition samplers, interpolators, and serialisation
utilities, all under a single namespace.

This directory documents every module: what it is for, the mathematical
conventions and derivations behind it, and the scientific context. It mirrors
the organisation of `src/MM`.

## Conventions

- Everything lives in `namespace MathEngine`.
- Fundamental containers (from `src/MM/typedefs/header.hpp`):
  - `Vec<T> = std::vector<T>`, with aliases `dVec` (doubles), `wVec` (size_t),
    `dcVec` (complex), `iVec`, `uVec`, `bVec`, `sVec`.
  - `Matrix<T>` (row-major, contiguous) with `dMatrix = Matrix<double>`.
- The ODE solvers exchange data through two structs:

  ```cpp
  using MyFunc  = std::function<void(double t, const dVec& y, dVec& dy/dt)>;
  struct ODESolverParameters { MyFunc derivative; CallBackFunc onStep; dVec initialConditions;
                               double t0, t1, dt; /* tolerances, order, ... */ };
  struct SolverResults { dMatrix solution; dVec timePoints; /* steps/errors histories */ };
  ```

- All code is header-only, C++23, and — for the molecular-dynamics engine —
  independent of any rendering backend (the GUI renders via raylib/ImPlot).

## Modules

| Module             | Description                                                            | Document                                   |
| ------------------ | ---------------------------------------------------------------------- | ------------------------------------------ |
| Models             | Kuramoto family, Ott–Antonsen reduction, molecular dynamics            | [models/](models/README.md)                |
| ODE solvers        | Runge–Kutta (RK1–RK4 + variants) and multistep (Adams) methods         | [solvers/](solvers/README.md)              |
| Initial conditions | Distributions, splay states, modules                                   | [initializers.md](initializers.md)         |
| Networks           | Random / Erdős–Rényi / small-world / modular / hierarchical topologies | [network-topology.md](network-topology.md) |
| Interpolators      | Lagrange and Newton divided-difference interpolation                   | [interpolators.md](interpolators.md)       |
| Types              | Vec / Matrix / solver callbacks                                        | [typedefs.md](typedefs.md)                 |
| Utility            | CSV / binary matrix & vector serialisation                             | [utility.md](utility.md)                   |

## The GUI

`src/AppState.hpp` implements the ImGui/raylib front-end. Its sidebar panels are
**model-aware**: the Model / Topology / Initial Conditions / Solver / Plot / Save
tabs adapt their controls to the selected `ModelType` (`Kuramoto`,
`OttAntonsen`, `MolecularDynamics`). See the model pages for the physics behind
each set of controls.

## Scientific references

- Kuramoto model: Kuramoto (1975); Strogatz (2000), *From Kuramoto to Crawford*.
- Ott–Antonsen ansatz: Ott & Antonsen (2008), *Chaos* **18**, 037113.
- Molecular dynamics: Allen & Tildesley, *Computer Simulation of Liquids* (2017);
  Frenkel & Smit, *Understanding Molecular Simulation* (2002).

# Random Walk

The random-walk model (`src/MM/models/RandomWalk/`, facade header
`src/MM/models/random-walk.hpp`) is a stochastic, Structure-of-Arrays engine that
is the direct analogue of the molecular-dynamics engine: it owns all walker state
in parallel arrays and is independent of any rendering backend.

## Model

`N` walkers start at a common origin — `startX`/`startY`, default `(0, 0)` — and
take one step per iteration. Each step is the product of a unit direction
(chosen by the move style) and a `stepSize` (distance per step, default `1.0`).

Every walker's position is accumulated **unwrapped** (the raw sum of its step
increments). The on-canvas coordinate is produced on demand via `WrapX`/`WrapY`
(for the periodic boundary) so that displacement and MSD remain physically
meaningful even as walkers wind around the torus.

```cpp
struct RandomWalkConfig
{
    size_t   numWalkers;              // number of walkers
    double   width, height;           // canvas size (for periodic/reflective)
    double   size;                    // draw size
    double   startX, startY;          // common starting position (default 0,0)
    double   stepSize;                // distance moved per step (default 1)
    WalkerMoveStyle moveStyle;        // step distribution
    BoundaryMode    boundary;         // edge behaviour
    uint64_t seed;                    // reproducible RNG seed
};
```

## Move styles

`WalkerMoveStyle` selects the step distribution (each unit step is then scaled by
`stepSize`):

| Style | Unit step |
|---|---|
| `Straight` | one of 4 cardinal directions |
| `Diagonal` | one of 4 diagonal directions |
| `StraightDiagonal` | one of 8 directions (4 cardinal + 4 diagonal) |
| `StraightWCenter` / `DiagonalWCenter` / `StraightDiagonalWCenter` | as above + a "stay" option |
| `StraightContinuous` | cardinal axis, continuous step in `[-1,1)` |
| `DiagonalContinuous` | diagonal, continuous step in `[-1,1)` |
| `StraightDiagonalContinuous` | independent continuous steps in `x` and `y` |

## Boundary modes

`BoundaryMode` controls the edges:

- `Periodic` — wrap around (torus); displacement is tracked unwrapped.
- `Reflective` — bounce off the walls.
- `Free` — no boundary; the view **auto-scales to fit the walkers**.

## Observables

`CollectObservables(rw, t)` (`analysis.hpp`) returns ensemble averages over the
`N` walkers at a given time:

- first/second moments of position: `⟨x⟩`, `⟨y⟩`, `⟨x²⟩`, `⟨y²⟩`, `⟨xy⟩`;
- variances/covariance: `σ_x² = ⟨x²⟩−⟨x⟩²`, `σ_y²`, `Cov(x,y) = ⟨xy⟩−⟨x⟩⟨y⟩`;
- distance/dispersion (relative to the starting point):
  `⟨|r|⟩`, `r_rms = √MSD`, `MSD = ⟨|r−r₀|²⟩`, `D = MSD/(4t)`
  (2D diffusion coefficient), and the radius of gyration
  `R_g² = ⟨x²+y²⟩ − ⟨x⟩² − ⟨y⟩²`.

### Diffusion scaling

For an isotropic unit-step walk the displacement converges (by the central-limit
theorem) to a 2D Gaussian of variance `MSD = t·⟨step²⟩`, hence:

- `MSD = stepSize² · t` (linear in time, quadratic in step size);
- `r_rms = √MSD ~ stepSize·√t`;
- the mean distance `⟨|r|⟩` follows a Rayleigh distribution,
  `⟨|r|⟩ → √(π/4)·√MSD ≈ 0.886·stepSize·√t`;
- `D = MSD/(4t) = stepSize²·⟨step²⟩/(4)` (constant).

These serve as validation targets: for `Straight` with `stepSize = 1`,
`⟨step²⟩ = 1` so `MSD = t` and `D = 1/4`. For `Diagonal`, `⟨step²⟩ = 2`.

## App integration

The GUI (`src/AppState.hpp`) adds RandomWalk as a first-class `ModelType`:

| Tab | Controls |
|---|---|
| Model | Move style (9 options), Move Size (step size), Walkers (N), Walker Size |
| Topology | Canvas Width/Height, Boundary Mode |
| Initial Condition | Start X, Start Y, Seed |
| Solver | Step size (dt), Steps, Stride (observable down-sampling) |
| Run | **Begin Simulation** (fresh run) and **Advance +N steps** (continue from current state) |
| Plot | MSD-vs-step (main) + trailing subplot (full resolution) + walker scatter (auto-scales for `Free`) |
| Save | `Observables.csv` (11 columns) and `FinalState.csv` |

The engine is persisted between runs, so **Advance** continues the walk from its
current positions and RNG state, applying the current move style and step size —
useful for evolving a lightweight system further without restarting.

## Reference

- Hughes, *Random Walks and Random Environments* (1995).
- Berg, *Random Walks in Biology* (1993).

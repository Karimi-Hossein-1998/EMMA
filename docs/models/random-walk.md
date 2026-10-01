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

Let the $t$-th step of a given walker be the random vector
$\mathbf{s}_t = (s_{t,x},\,s_{t,y})$, with independent, identically distributed
steps of mean zero, $\langle\mathbf{s}\rangle=0$, and isotropic second moment
$\langle s_x^2\rangle=\langle s_y^2\rangle=\tfrac12\langle\lVert\mathbf{s}\rVert^2\rangle$.
The position after $t$ steps is the sum $\mathbf{r}(t)=\sum_{\tau=1}^{t}\mathbf{s}_\tau$.
By the **central-limit theorem**, each Cartesian component converges to a Gaussian
of variance $t\langle s_x^2\rangle$, so the displacement is a 2D Gaussian

$$
P(\mathbf{r}) = \frac{1}{2\pi\sigma^2}\exp\!\left(-\frac{\lVert\mathbf{r}\rVert^{2}}{2\sigma^2}\right),
\qquad
\sigma^2 = t\,\langle s_x^2\rangle = \tfrac{t}{2}\langle\lVert\mathbf{s}\rVert^2\rangle.
$$

The mean-squared displacement is therefore

$$
\mathrm{MSD}(t) = \langle\lVert\mathbf{r}\rVert^2\rangle = 2\sigma^2 = t\,\langle\lVert\mathbf{s}\rVert^2\rangle,
$$

i.e. **linear in time** and quadratic in step size. The radial distance
$\lVert\mathbf{r}\rVert$ then follows a **Rayleigh distribution**

$$
P(\lVert\mathbf{r}\rVert) = \frac{\lVert\mathbf{r}\rVert}{\sigma^2}\exp\!\left(-\frac{\lVert\mathbf{r}\rVert^{2}}{2\sigma^2}\right),
$$

whose mean is

$$
\langle\lVert\mathbf{r}\rVert\rangle
= \sigma\sqrt{\frac{\pi}{2}}
= \sqrt{\frac{\pi}{4}}\,\sqrt{\mathrm{MSD}}
\approx 0.886\,\mathrm{stepSize}\,\sqrt{t}.
$$

Finally, Einstein's relation in $d$ dimensions, $\mathrm{MSD}=2d\,D\,t$, gives the
diffusion coefficient

$$
D = \frac{\mathrm{MSD}}{4t} = \frac{\langle\lVert\mathbf{s}\rVert^2\rangle}{4},
$$

a constant. These serve as validation targets: for `Straight` with
`stepSize = 1`, $\langle\lVert\mathbf{s}\rVert^2\rangle = 1$ so `MSD = t` and
`D = 1/4`; for `Diagonal`, $\langle\lVert\mathbf{s}\rVert^2\rangle = 2$.

## App integration

The GUI (`src/AppState.hpp`) adds RandomWalk as a first-class `ModelType`:

| Tab | Controls |
|---|---|
| Model | **Dimensions (2D/3D)**, Move style (9 / 21 options), Move Size, Walkers (N), Walker Size |
| Topology | Canvas Width/Height (+ Depth in 3D), Boundary Mode |
| Initial Condition | Start X, Start Y (+ Start Z in 3D), Seed |
| Solver | Step size (dt), Steps, Stride (observable down-sampling) |
| Run | **Begin Simulation** (fresh run); **Advanced (Continue Run)** → Move Style, Move Size, Steps, and **Advance +N steps** |
| Plot | MSD-vs-step (main) + trailing subplot (full resolution) + walker scatter (2D auto-scale) / 3D viewport |
| Save | `Observables.csv` (11 / 15 columns) and `FinalState.csv` (2 / 3 columns) |

The engine is persisted between runs, so **Advance** continues the walk from its
current positions and RNG state, applying the move style / move size / step count
set in the Run panel's *Advanced* section — useful for evolving a lightweight
system further without restarting.

## 3D variant

`RandomWalk3D` (`src/MM/models/RandomWalk3D/`, facade
`src/MM/models/random-walk3d.hpp`) is a separate mirror of the 2D engine with a
third axis. Positions are `posX/posY/posZ`, the box is
`width × height × depth`, and the boundary modes apply per axis.

### 3D move styles

`WalkerMoveStyle3D` adds the third dimension, decomposing the 26-neighbourhood by
the number of non-zero step components:

| Style | Directions | $\langle\lVert\mathbf{s}\rVert^2\rangle$ |
|---|---|---|
| `Straight` | 6 cardinal | 1 |
| `PlaneDiagonal` | 12 face diagonals | 2 |
| `Diagonal` | 8 body diagonals | 3 |
| `FullDiagonal` | 20 (plane + body) | 2.4 |
| `StraightPlaneDiagonal` | 18 | 5/3 |
| `StraightDiagonal` | 14 | 15/7 |
| `StraightFullDiagonal` | 26 (full neighbourhood) | 27/13 |

Each also has a `…WCenter` (adds a "stay") and a `…Continuous`
(magnitude in `[0,1)`) variant — 21 styles in total.

### 3D observables

`CollectObservables3D` returns the same moments extended to three dimensions
(`⟨z⟩`, `⟨z²⟩`, `Cov(xz)`, `Cov(yz)`, `σ_z²`, …). The diffusion coefficient uses
Einstein's relation in 3D:

$$
D = \frac{\mathrm{MSD}}{6t},
$$

and the radial distance follows the 3D Maxwell limit
$\langle\lVert\mathbf{r}\rVert\rangle \to \sqrt{8/(3\pi)}\,\sqrt{\mathrm{MSD}}
\approx 0.921\,\sqrt{\mathrm{MSD}}$.

### Rendering

In the app the walker positions are shown in a **3D viewport** (raylib
`Camera3D` rendered into a `RenderTexture2D`, embedded in an ImGui window with
mouse orbit/zoom), while the observables remain ordinary 2D ImPlot plots. The
trail is a **fading voxel cloud** (age-damped, matching the 2D fading trail).

## Reference

- Hughes, *Random Walks and Random Environments* (1995).
- Berg, *Random Walks in Biology* (1993).

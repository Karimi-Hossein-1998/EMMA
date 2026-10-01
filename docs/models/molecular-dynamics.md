# Molecular Dynamics (2D)

The molecular-dynamics (MD) model simulates a collection of $N$ classical
particles in two dimensions interacting through a pairwise potential. It is a
*condensed-matter* model — the analogue, for solids and liquids, of the Kuramoto
order parameter is the **bond-orientational order** $\psi_6$, which the engine
computes alongside the thermodynamics.

Files (all in `namespace MathEngine`, SDL/raylib-free):

| Header | Contents |
|---|---|
| `models/MolecularDynamics/potential.hpp` | `PotentialType` {LennardJones, WCA, Morse, SoftSphere, Yukawa, Coulomb2D} + shifted-force kernel |
| `models/MolecularDynamics/initial-conditions.hpp` | `InitialConditionType` {SquareLattice, HexagonalLattice, Random, TwoPhaseSlab, BinaryMixture} |
| `models/MolecularDynamics/engine.hpp` | `MolecularDynamics` (SoA) engine + `MDConfig` + `IntegratorType` |
| `models/MolecularDynamics/thermostats.hpp` | `ThermostatType` + barostat + pressure |
| `models/MolecularDynamics/analysis.hpp` | $\psi_4/\psi_6$, RDF, MSD, energies, observables |

---

## 1. Conventions and units

Particles have mass $m_i$, position $\mathbf{r}_i=(x_i,y_i)$, velocity
$\mathbf{v}_i$ and acceleration $\mathbf{a}_i$. State is stored in
**Structure-of-Arrays** layout (`posX`, `posY`, `velX`, `velY`, `accX`, `accY`,
`mass`, …) for cache efficiency and OpenMP vectorisation.

The energy scale is set by the potential parameters $\varepsilon$ (and $\sigma$
for the length scale). The temperature is defined through the 2D kinetic-energy
equipartition, i.e. with Boltzmann's constant set to $k_B=1$,

$$
T = \frac{1}{N}\sum_{i=1}^{N} \frac{1}{2}\,m_i\,\lVert \mathbf{v}_i\rVert^2
  = \frac{E_{\mathrm{kin}}}{N}.
$$

This is the convention used by `ComputeKineticEnergy()` and the thermostats.

**Why $T=E_{\mathrm{kin}}/N$ (equipartition).** In the canonical ensemble the
velocity distribution is the Maxwell–Boltzmann distribution

$$
P(\mathbf{v}_i) \propto \exp\!\left(-\frac{m_i \lVert\mathbf{v}_i\rVert^{2}}{2T}\right),
$$

a product of two independent Gaussians of variance $T/m_i$ per component. The
mean kinetic energy of a *single* component is therefore

$$
\left\langle \tfrac{1}{2}m_i v_{i\alpha}^{2}\right\rangle
= \frac{T}{2},
\qquad \alpha \in \{x,y\},
$$

so a particle with $d=2$ translational degrees of freedom carries
$\langle \tfrac12 m_i \lVert\mathbf{v}_i\rVert^{2}\rangle = 2\cdot T/2 = T$. This
is the equipartition theorem: each quadratic degree of freedom contributes
$T/2$ (with $k_B=1$). Summing over $N$ particles gives
$E_{\mathrm{kin}} = N T$, i.e. $T = E_{\mathrm{kin}}/N$.

---

## 2. Pair potentials

All potentials are expressed through `PairForceEnergy(r^2, type, params, fFactor, pe)`,
which returns a scalar $f$ such that the force **on** particle $i$ (located at the
origin of the separation vector $\mathbf{r}_{ij}=\mathbf{r}_j-\mathbf{r}_i$) is

$$
\mathbf{F}_{ij} = f\,(x_{ij},\,y_{ij}), \qquad f = \frac{1}{r}\,\frac{\mathrm{d}U}{\mathrm{d}r}.
$$

### 2.1 Lennard-Jones 12-6

$$
U(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right],
\qquad
\frac{\mathrm{d}U}{\mathrm{d}r} = \frac{24\varepsilon}{r}\left[\left(\frac{\sigma}{r}\right)^{6} - 2\left(\frac{\sigma}{r}\right)^{12}\right].
$$

The $r^{-6}$ term is the attractive van der Waals (London dispersion) tail; the
$r^{-12}$ term is a convenient model of the Pauli repulsion. The well minimum is
at $r_{\min}=2^{1/6}\sigma$ with depth $-\varepsilon$.

### 2.2 WCA (Weeks–Chandler–Andersen)

The purely repulsive core, obtained by truncating LJ at its minimum and shifting
the energy up by $\varepsilon$ so $U(r_{\min})=0$:

$$
U_{\mathrm{WCA}}(r) =
\begin{cases}
4\varepsilon\left[\left(\dfrac{\sigma}{r}\right)^{12} - \left(\dfrac{\sigma}{r}\right)^{6}\right] + \varepsilon, & r < 2^{1/6}\sigma,\\[0.5em]
0, & r \ge 2^{1/6}\sigma.
\end{cases}
$$

Because $U$ is cut at its minimum, the force is already exactly zero at the
cutoff — WCA needs no force shifting.

### 2.3 Morse

$$
U(r) = D\left[\left(1 - e^{-\alpha(r-r_0)}\right)^2 - 1\right]
     = D\left(e^{-2\alpha(r-r_0)} - 2 e^{-\alpha(r-r_0)}\right),
\qquad
\frac{\mathrm{d}U}{\mathrm{d}r} = 2D\alpha\,(1-e^{-\alpha(r-r_0)})\,e^{-\alpha(r-r_0)}.
$$

Here $D\equiv\varepsilon$ is the well depth, $r_0\equiv\sigma$ the equilibrium
bond length, and $\alpha$ the stiffness (width) of the well. The Morse potential
models chemical bonds and admits the discrete vibrational spectrum
$E_n = -D + \hbar\omega(n+\tfrac12)$ in quantum mechanics.

### 2.4 Soft sphere (power-law repulsion)

$$
U(r) = \varepsilon\left(\frac{\sigma}{r}\right)^{n},
\qquad
\frac{\mathrm{d}U}{\mathrm{d}r} = -\frac{n\,\varepsilon}{r}\left(\frac{\sigma}{r}\right)^{n}.
$$

A purely repulsive inverse-power interaction; the exponent $n$ (default 9) sets the
softness. It models steric / soft-particle repulsion (inverse-power fluids,
Hertzian contacts, granular matter).

### 2.5 Yukawa (screened Coulomb)

$$
U(r) = \varepsilon\,\frac{\sigma}{r}\,e^{-\kappa r},
\qquad
\frac{\mathrm{d}U}{\mathrm{d}r}
= -\varepsilon\,\sigma\,e^{-\kappa r}\left(\frac{\kappa}{r} + \frac{1}{r^2}\right).
$$

The Debye–Hückel screened Coulomb interaction; $\kappa$ is the inverse screening
length (set by the ionic strength of the background). It describes charged
colloids in an electrolyte, dusty plasmas, and Yukawa fluids.

### 2.6 Coulomb 2D (soft-core logarithmic)

$$
U(r) = \varepsilon\,\ln\!\frac{r+\sigma}{r},
\qquad
\frac{\mathrm{d}U}{\mathrm{d}r} = -\frac{\varepsilon\,\sigma}{r\,(r+\sigma)}.
$$

The two-dimensional Coulomb (logarithmic) interaction, regularised by a soft core
of width $\sigma$ that removes the $\ln$ singularity at $r=0$. In 2D the bare
Coulomb force decays as $1/r$; this soft-core form describes point vortices, 2D
plasmas, and charged discs. Here $\varepsilon$ is the product of the two "charge"
strengths.

### 2.7 Shifted-force cutoff

A plain truncation makes the force discontinuous at $r_c$, which leaks energy in
NVE. The engine uses a **linear shifted-force** cutoff: for $r<r_c$,

$$
f_{\mathrm{sf}}(r) = f(r) - \frac{f_c}{r_c}\,r,
\qquad
U_{\mathrm{sf}}(r) = U(r) - U(r_c) - \frac{f_c}{2 r_c}\left(r^2 - r_c^2\right),
$$

with $f_c = \mathrm{d}U/\mathrm{d}r\,\big|_{r_c}$ and $U_c=U(r_c)$. Then
$f_{\mathrm{sf}}(r_c)=0$ and $U_{\mathrm{sf}}(r_c)=0$ continuously. The shift is
linear in $r$ (equivalently in $r^2$), so it needs **no square root** in the inner
loop — unlike the constant-force shift — keeping Lennard-Jones/WCA branch-free.

---

## 3. Initial conditions

`MakeInitialState` places particles and assigns Maxwell–Boltzmann velocities.

- **Square / Hexagonal lattice** — grid sized to the box via `GridDims` so the
  lattice always fills the box without PBC wrap-overlap.
- **Random** — uniform positions with a minimum separation (rejection sampling).
- **Two-phase slab** — a dense band in the middle (a phase-separation seed).
- **Binary mixture** — two interleaved species with mass/radius ratios.

Velocities are drawn as $v_\alpha \sim \mathcal{N}(0,\, T/m_i)$ per component,
the centre-of-mass velocity is subtracted (total momentum $=0$), and the ensemble
is rescaled so the kinetic temperature equals the target exactly.

---

## 4. Integration

### 4.1 Velocity-Verlet (default)

The velocity-Verlet (a.k.a. leap-frog in staggered form) is symplectic and
time-reversible, which keeps the total energy bounded in NVE:

$$
\begin{aligned}
\mathbf{v}(t+\tfrac{\mathrm{d}t}{2}) &= \mathbf{v}(t) + \tfrac{\mathrm{d}t}{2}\,\mathbf{a}(t),\\
\mathbf{r}(t+\mathrm{d}t) &= \mathbf{r}(t) + \mathrm{d}t\,\mathbf{v}(t+\tfrac{\mathrm{d}t}{2}),\\
\mathbf{a}(t+\mathrm{d}t) &= \mathbf{F}(\mathbf{r}(t+\mathrm{d}t))/m,\\
\mathbf{v}(t+\mathrm{d}t) &= \mathbf{v}(t+\tfrac{\mathrm{d}t}{2}) + \tfrac{\mathrm{d}t}{2}\,\mathbf{a}(t+\mathrm{d}t).
\end{aligned}
$$

**Why not RK4?** RK4 is not symplectic: over long runs it dissipates energy. The
MD engine therefore exposes its own integrator rather than routing through the
generic RK/AB solvers. (`StepLeapfrog` is the equivalent position-Verlet form.)

**Symplectic nature.** Velocity-Verlet is a second-order symplectic integrator:
it exactly preserves a *shadow Hamiltonian* $\tilde{H}$ close to $H$, so the
energy only oscillates with amplitude $\mathcal{O}(\mathrm{d}t^2)$ instead of
drifting linearly. With the shifted-force cutoff, the measured NVE drift over
$5\times10^3$ steps is $\sim 10^{-4}$ (see the validation harness).

---

## 5. Thermostats (processes)

Applied on top of (or instead of) an NVE step; target temperature $T_0$.

| Thermostat | Update | Notes |
|---|---|---|
| Rescale | $\mathbf{v}\leftarrow \sqrt{T_0/T}\,\mathbf{v}$ | instantaneous; exact but non-physical |
| Berendsen | $\lambda=\sqrt{1+\frac{\mathrm{d}t}{\tau}\left(\frac{T_0}{T}-1\right)}$, $\mathbf{v}\leftarrow\lambda\mathbf{v}$ | weak coupling, relaxes to $T_0$ over $\tau$ |
| Andersen | with prob. $\nu\,\mathrm{d}t$, redraw $\mathbf{v}$ from MB($T_0$) | stochastic collision bath |
| Langevin | $\mathbf{v}\leftarrow e^{-\gamma\,\mathrm{d}t}\mathbf{v}+\sqrt{(1-e^{-2\gamma\,\mathrm{d}t})\,T_0/m}\,\boldsymbol{\xi}$ | friction + FDT noise (exact OU update) |
| Nosé–Hoover | extended variable $\xi$: $\dot\xi=\dfrac{2E_{\mathrm{kin}}-N_{\mathrm{dof}}T_0}{Q}$, $\dot{\mathbf{v}}=\mathbf{F}/m-\xi\mathbf{v}$, $Q=N_{\mathrm{dof}}T_0\tau^2$ | canonical sampling |

### Langevin thermostat — fluctuation–dissipation theorem

The Langevin thermostat couples each particle to a heat bath through a friction
force and a fluctuating force,

$$
m\,\frac{\mathrm{d}\mathbf{v}}{\mathrm{d}t}
= -\gamma\,\mathbf{v} + \boldsymbol{\eta}(t),
\qquad
\langle \boldsymbol{\eta}(t)\rangle = 0,
\qquad
\langle \eta_\alpha(t)\,\eta_\beta(t')\rangle
= 2\gamma\,T\,\delta_{\alpha\beta}\,\delta(t-t'),
$$

where $\gamma$ is the friction and $\boldsymbol{\eta}$ is white noise. The
fluctuation–dissipation theorem (FDT) fixes the noise amplitude to $2\gamma T$ —
the same $\gamma$ that appears in the damping — so that the stationary
distribution of $\mathbf{v}$ is the Maxwell–Boltzmann distribution at temperature
$T$. Equivalently, each velocity component is an Ornstein–Uhlenbeck process whose
stationary variance is $T/m$.

The engine uses the **exact** solution of the free-particle Ornstein–Uhlenbeck
process over one step (unconditionally stable for any $\gamma$),

$$
\mathbf{v}(t+\mathrm{d}t)
= e^{-\gamma\,\mathrm{d}t}\,\mathbf{v}(t)
+ \sqrt{\frac{T}{m}\left(1-e^{-2\gamma\,\mathrm{d}t}\right)}\;\boldsymbol{\xi},
\qquad \xi_\alpha \sim \mathcal{N}(0,1),
$$

whose variance is exactly $T/m$ for any step size.

### Nosé–Hoover thermostat — the extended system

The Nosé–Hoover thermostat is deterministic: it appends a single "friction"
variable $\xi$ whose dynamics drive the kinetic temperature toward $T_0$. It
derives from the extended Hamiltonian

$$
\mathcal{H}_{\mathrm{NH}}
= \sum_i \frac{\lVert\mathbf{p}_i\rVert^{2}}{2m_i} + V(\mathbf{q})
+ \frac{p_s^{2}}{2Q} + (N_{\mathrm{dof}}+1)\,T_0\,s,
$$

where $s$ is an extra "position" (with momentum $p_s$ and inertia $Q$) that
rescales time, and $N_{\mathrm{dof}}$ is the number of degrees of freedom. After
the standard change of variables to the friction coefficient
$\xi = p_s/Q$ and the thermostat inertia $Q = N_{\mathrm{dof}}\,T_0\,\tau^{2}$
(with relaxation time $\tau$), Hamilton's equations become

$$
\dot{\mathbf{v}}_i = \frac{\mathbf{F}_i}{m_i} - \xi\,\mathbf{v}_i,
\qquad
\dot{\xi} = \frac{\sum_i m_i\lVert\mathbf{v}_i\rVert^{2}
             - N_{\mathrm{dof}}\,T_0}{Q}
         = \frac{2 E_{\mathrm{kin}} - N_{\mathrm{dof}}\,T_0}{Q}.
$$

These are exactly the equations integrated by `StepNoseHoover`. The extended
system preserves the canonical distribution in the $(q,p,s,p_s)$ phase space, so
the physical marginals sample the canonical ensemble at temperature $T_0$.

---

## 6. Barostat (pressure control)

The Berendsen (weak-coupling) barostat rescales positions and the box so the
instantaneous pressure relaxes toward $P_0$ over a time $\tau_P$:

$$
\mu = \left[1 - \frac{\mathrm{d}t}{\tau_P}(P_0 - P)\right]^{1/2},
\qquad
(\mathbf{r}_i,\,L_x,\,L_y) \leftarrow \mu\,(\mathbf{r}_i,\,L_x,\,L_y).
$$

### Virial theorem derivation

Define the *virial* $G = \sum_i \mathbf{p}_i\cdot\mathbf{r}_i$. Its time
derivative is

$$
\frac{\mathrm{d}G}{\mathrm{d}t}
= \sum_i \left(\dot{\mathbf{p}}_i\cdot\mathbf{r}_i
             + \mathbf{p}_i\cdot\dot{\mathbf{r}}_i\right)
= \sum_i \mathbf{F}_i^{\mathrm{tot}}\cdot\mathbf{r}_i
  + \sum_i \frac{\lVert\mathbf{p}_i\rVert^{2}}{m_i}
= \sum_i \mathbf{F}_i^{\mathrm{tot}}\cdot\mathbf{r}_i
  + 2 E_{\mathrm{kin}},
$$

using Newton's second law $\dot{\mathbf{p}}_i = \mathbf{F}_i^{\mathrm{tot}}$. In
a stationary state $\langle \mathrm{d}G/\mathrm{d}t\rangle = 0$ (the virial is
bounded), and equipartition gives $\langle E_{\mathrm{kin}}\rangle = \tfrac{d}{2}N T$, so

$$
\left\langle \sum_i \mathbf{F}_i^{\mathrm{tot}}\cdot\mathbf{r}_i \right\rangle
= -2\langle E_{\mathrm{kin}}\rangle = -d\,N T .
$$

Split the total force into internal pairwise forces and the force exerted by the
walls, $\mathbf{F}_i^{\mathrm{tot}} = \mathbf{F}_i^{\mathrm{int}} + \mathbf{F}_i^{\mathrm{wall}}$.
The wall contribution is the (negative) pressure work,
$\langle \sum_i \mathbf{F}_i^{\mathrm{wall}}\cdot\mathbf{r}_i\rangle = -d\,P V$,
so

$$
\langle W\rangle - d P V = -d N T,
\qquad
W := \sum_i \mathbf{r}_i\cdot\mathbf{F}_i^{\mathrm{int}},
$$

i.e.

$$
P = \rho T + \frac{1}{dV}\langle W\rangle,
\qquad \rho = N/V.
$$

For pairwise forces, $W = \sum_{i<j}(\mathbf{r}_i-\mathbf{r}_j)\cdot\mathbf{f}_{ij}$
(Newton's third law), and with $\mathbf{r}_{ij}=\mathbf{r}_j-\mathbf{r}_i$ the
2D pressure is evaluated as

$$
P = \frac{N T - \tfrac{1}{2}\sum_{i<j}\mathbf{r}_{ij}\cdot\mathbf{f}_{ij}}{V},
\qquad V = L_x L_y.
$$

Repulsive forces have $\mathbf{r}_{ij}\cdot\mathbf{f}_{ij}<0$, so the minus sign
makes the correction positive: repulsion raises the pressure above the ideal-gas
value.

---

## 7. Analysis

### 7.1 Bond-orientational order $\psi_n$

The local bond-orientational order of particle $i$ is

$$
\psi_n(i) = \frac{1}{N_i}\sum_{j\in\mathrm{nbr}(i)} e^{i n\theta_{ij}},
$$

where $\theta_{ij}=\mathrm{atan2}(y_j-y_i,\,x_j-x_i)$ and neighbours are those
within a cutoff (default $1.4\sigma$). The global order parameter is

$$
\Psi_n = \left|\frac{1}{N}\sum_{i=1}^{N}\psi_n(i)\right|.
$$

$\Psi_6\approx 1$ for a perfect triangular/hexatic crystal and $\approx 0$ for a
liquid; $\Psi_4$ plays the same role for a square lattice. This is the
condensed-matter analogue of the Kuramoto $\rho$.

### 7.2 Radial distribution function

$$
g(r) = \frac{\langle \text{number of pairs at distance }r \rangle}
            {\tfrac{1}{2}N\,\rho\,(2\pi r\,\mathrm{d}r)},
\qquad \rho = N/V.
$$

$g(r)$ peaks at the lattice spacings of a solid and tends to 1 for an ideal gas.

### 7.3 Mean-squared displacement

$$
\mathrm{MSD}(t) = \left\langle \lVert \mathbf{r}_i(t) - \mathbf{r}_i(0)\rVert^2 \right\rangle,
$$

computed with the minimum-image convention (valid while the displacement is
$< L/2$). In 2D the diffusion coefficient follows from Einstein's relation
$\mathrm{MSD}(t)\sim 4 D t$.

---

## 8. API reference (namespace `MathEngine`)

**Enums** — `PotentialType`, `InitialConditionType`, `IntegratorType`, `ThermostatType`.

**Configuration**

```cpp
struct MDConfig { size_t numParticles; double width, height; double mass, radius,
                  massRatio, radiusRatio; double sigma, epsilon, cutoffCoeff, morseAlpha;
                  double temperature, restitution, minSeparation; size_t seed;
                  bool periodicBoundaryCondition, bounce, hardSphereCollisions;
                  PotentialType potential; InitialConditionType initialCondition; };
```

**Engine** — `MolecularDynamics`:

| Member | Purpose |
|---|---|
| `Step(dt)` / `StepLeapfrog(dt)` | NVE integrators |
| `CalculateAccelerations()` | O(N²) force loop (OpenMP when available) |
| `ComputeKineticEnergy()` / `ComputePotentialEnergy()` / `ComputeVirialSum()` | observables |
| `SetTemperature(T)` / `ScaleVelocities(s)` / `ScalePositionsAndBox(s)` | control |

**Processes** — `ApplyVelocityRescale`, `ApplyBerendsen`, `ApplyAndersen`,
`ApplyLangevin`, `StepNoseHoover`, `ApplyBerendsenBarostat`, `ComputePressure`.

**Analysis** — `ComputeBondOrientationalOrder`, `ComputeRDF`, `ComputeMSD`,
`CollectObservables` (returns `Observables{time, temperature, kineticEnergy,
potentialEnergy, totalEnergy, pressure, psi4, psi6, msd}`).

## Use cases

- **Crystallisation / melting** — cool a liquid below the melting point and watch
  $\Psi_6$ rise toward 1; the GUI streams $\Psi_6$ to the order-parameter plot.
- **Equation-of-state** — vary density/temperature and read off the virial
  pressure; use the barostat for NPT runs.
- **Diffusion** — extract $D$ from the MSD slope.
- **Structure** — $g(r)$ and $\Psi_4/\Psi_6$ discriminate solid / liquid /
  hexatic phases and reveal the lattice symmetry.
- **Soft / screened matter** — SoftSphere (inverse-power fluids), Yukawa (charged
  colloids in an electrolyte), and Coulomb2D (2D plasmas, point vortices) extend
  the model beyond van-der-Waals fluids.

## References

- Allen & Tildesley, *Computer Simulation of Liquids* (2017).
- Frenkel & Smit, *Understanding Molecular Simulation* (2002).
- Weeks, Chandler & Andersen (1971), *J. Chem. Phys.* **54**, 5237.
- Berendsen et al. (1984), *J. Chem. Phys.* **81**, 3684.
- Hoover (1985), *Phys. Rev. A* **31**, 1695; Nosé (1984), *J. Chem. Phys.* **81**, 511.

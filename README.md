# Mathematical Modelling GUI

This project aims to create a **_G_**_raphical_ **_U_**_ser_ **_I_**_nterface_ for the project [Mathematical Modelling Toolkit](https://github.com/Karimi-Hossein-1998/math-mod-tk.git).

## Features

+ Added Fonts: *Inter*, *BonaNova Regular*, *Brawler Regular*, *ClearSans Regular*.
+ Added *Style Editor* (from *ImGui*).
+ Added Model Picker. (*Kuramoto*, *Ott-Antonsen*, and *Molecular Dynamics* model families).
+ Added Initial Value Picker. (For *Kuramoto* models, two initial arrays are needed. Initial phases, and intrinsic frequencies.), (The two *Kuramoto* variants need an adjacency matrix.)
+ Added Network Structure Picker. (For *Kuramoto* models).
+ Added Solver Model Picker (Kinda Complete!). RK1 (Euler), RK2 (midpoint), RK3, RK4 (standard), RK4 (3/8), RK4 (Gill), RK4 (Ralston).
+ Added Live Plotting (Kinda Complete!).
+ Added Solver Picker.
+ Added Plotting.
- Added *Contiguous* Matrix container from [NumCpp](https://github.com/Karimi-Hossein-1998/NumCpp).
- Added *Multistep* solvers. *Adams-Bashforth* **predictor** and *Adams-Bashforth-Moulton* **predictor-corrector** methods are now available.
- Added **SIMD** directives.
- Added **Chronometer** for simulation runs.
- Added Saving to Files.
- Added the **Ott–Antonsen** dimensionality-reduction model (single- and multi-community).
- Added a **2D Molecular Dynamics** model: Lennard-Jones / WCA / Morse potentials (shifted-force cutoff), lattice / random / slab / binary initial conditions, symplectic velocity-Verlet & leapfrog integrators, thermostats (Rescale, Berendsen, Andersen, Langevin, Nosé–Hoover) plus a Berendsen barostat, and bond-orientational order ($\psi_4$/$\psi_6$) / RDF / MSD analysis.
- Added a **Random Walk** model (Structure-of-Arrays engine): 9 move styles (straight/diagonal/continuous ± center), periodic / reflective / free boundary modes, and full moment & diffusion observables ($\langle x\rangle,\langle y\rangle,\langle x^2\rangle,\langle y^2\rangle$, variances, MSD, $D=\mathrm{MSD}/4t$). The `Free` boundary auto-scales the walker view to fit.
- Added **model-aware** sidebars and plotting (each panel and plot adapts to the selected model).

## Documentation

Full documentation — mathematical conventions, derivations, and use cases — lives
in [`docs/`](docs/README.md), including the
[molecular-dynamics reference](docs/models/molecular-dynamics.md).

## TODO

+ [ ] Add Image and GIF Save Options.
+ [ ] Add running Simulations in Loops, Forward-Backward Loops, and Loop inside Loops.




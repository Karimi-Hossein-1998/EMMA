# Mathematical Modelling GUI

This project aims to create a **_G_**_raphical_ **_U_**_ser_ **_I_**_nterface_ for the project [Mathematical Modelling Toolkit](https://github.com/Karimi-Hossein-1998/math-mod-tk.git).

## Features

+ Added Fonts: *Inter*, *BonaNova Regular*, *Brawler Regular*, *ClearSans Regular*.
+ Added *Style Editor* (from *ImGui*).
+ Added Model Picker. (*Kuramoto*, *Ott-Antonsen*, *Molecular Dynamics*, and *Random Walk* model families).
+ Added Initial Value Picker. (For *Kuramoto* models, two initial arrays are needed. Initial phases, and intrinsic frequencies.), (The two *Kuramoto* variants need an adjacency matrix.)
+ Added Network Structure Picker. (For *Kuramoto* models).
+ Added Solver Model Picker (Kinda Complete!). RK1 (Euler), RK2 (midpoint), RK3, RK4 (standard), RK4 (3/8), RK4 (Gill), RK4 (Ralston).
+ Added Live Plotting (Kinda Complete!).
+ Added Solver Picker.
+ Added Plotting.
- Added *Contiguous* Matrix container from [NumCpp](https://github.com/Karimi-Hossein-1998/NumCpp).
- Added *Multistep* solvers. *Adams-Bashforth* **predictor** and *Adams-Bashforth-Moulton* **predictor-corrector** methods are now available.
- Added **OpenMP** parallelisation of the molecular-dynamics force loops (native builds; the Emscripten web build stays single-threaded).
- Added **Chronometer** for simulation runs.
- Added Saving to Files.
- Added the **Ott–Antonsen** dimensionality-reduction model (single- and multi-community).
- Added **2D and 3D Molecular Dynamics** models: Lennard-Jones / WCA / Morse / Soft-sphere / Yukawa / Coulomb potentials (shifted-force cutoff), lattice (square/hexagonal in 2D, SC/BCC/FCC in 3D) / random / slab / binary initial conditions, symplectic velocity-Verlet & leapfrog integrators, thermostats (Rescale, Berendsen, Andersen, Langevin, Nosé–Hoover) plus a Berendsen barostat, and bond-orientational order ($\psi_4$/$\psi_6$ in 2D, Steinhardt $Q_4$/$Q_6$ in 3D) / RDF / MSD analysis.
- Added **2D and 3D Random Walk** models (Structure-of-Arrays engines): 9 move styles in 2D / 21 in 3D (straight/diagonal/continuous ± center), periodic / reflective / free boundary modes, and full moment & diffusion observables ($\langle x\rangle,\langle y\rangle,\langle x^2\rangle,\langle y^2\rangle$, variances, MSD, $D=\mathrm{MSD}/4t$ in 2D, $D=\mathrm{MSD}/6t$ in 3D). The `Free` boundary auto-scales the walker view to fit.
- Added **model-aware** sidebars and plotting (each panel and plot adapts to the selected model).

## Building

Native (multi-threaded):

```sh
cmake -S . -B build
cmake --build build
```

Web (Emscripten → `EMMA.html`, `EMMA.js`, `EMMA.wasm`):

```sh
cmake -S . -B build-web \
  -DCMAKE_TOOLCHAIN_FILE="$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake"
cmake --build build-web
```

The web build is single-threaded (OpenMP and `-march=native` are automatically
disabled for Emscripten); deploy `build-web/EMMA.{html,js,wasm}` to GitHub Pages.

## Documentation

Full documentation — mathematical conventions, derivations, and use cases — lives
in [`docs/`](docs/README.md), including the
[molecular-dynamics reference](docs/models/molecular-dynamics.md).

## TODO

+ [ ] Add Image and GIF Save Options.
+ [ ] Add running Simulations in Loops, Forward-Backward Loops, and Loop inside Loops.




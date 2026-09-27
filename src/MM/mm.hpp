#pragma once
// -----------------------------------------------------------------------------
// Umbrella header for the Mathematical Modelling (MM) toolkit.
//
// Aggregates the actively-maintained components. The legacy adaptive-step and
// delay (DDE) solver headers (and the old `SolverParameters`-based utilities)
// are intentionally excluded: they still reference the removed `SolverParameters`
// type (now `ODESolverParameters`) and are not part of the current build.
// -----------------------------------------------------------------------------

#include "initializers/initials.hpp"
#include "models/models.phase-oscillators.hpp"
#include "network/topology.hpp"

#include "solvers/ODE/rk/explicit/rk1-solver.hpp"
#include "solvers/ODE/rk/explicit/rk2-solver.hpp"
#include "solvers/ODE/rk/explicit/rk3-solver.hpp"
#include "solvers/ODE/rk/explicit/rk4-solver.hpp"
#include "solvers/ODE/rk/explicit/rk4-variants.hpp"
#include "solvers/ODE/multistep/ab-solver.hpp"
#include "solvers/ODE/multistep/abm-solver.hpp"

#include "utility/write.hpp"

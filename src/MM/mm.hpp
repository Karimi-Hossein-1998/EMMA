#pragma once
// -----------------------------------------------------------------------------
// Umbrella header for the Mathematical Modelling (MM) toolkit.
//
// Aggregates the actively-maintained components. The legacy adaptive-step and
// delay (DDE) solver headers (which referenced the removed `SolverParameters`
// type) have been removed; re-add them if they are ever revived.
// -----------------------------------------------------------------------------

#include "initializers/initials.hpp"
#include "models/models.phase-oscillators.hpp"
#include "models/molecular-dynamics.hpp"
#include "models/random-walk.hpp"
#include "models/random-walk3d.hpp"
#include "network/topology.hpp"

#include "solvers/ODE/rk/explicit/rk1-solver.hpp"
#include "solvers/ODE/rk/explicit/rk2-solver.hpp"
#include "solvers/ODE/rk/explicit/rk3-solver.hpp"
#include "solvers/ODE/rk/explicit/rk4-solver.hpp"
#include "solvers/ODE/rk/explicit/rk4-variants.hpp"
#include "solvers/ODE/multistep/ab-solver.hpp"
#include "solvers/ODE/multistep/abm-solver.hpp"

#include "utility/write.hpp"

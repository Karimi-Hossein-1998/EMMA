#pragma once
#include "ODE/rk/explicit/rk1-solver.hpp"
#include "ODE/rk/explicit/rk2-solver.hpp"
#include "ODE/rk/explicit/rk3-solver.hpp"
#include "ODE/rk/explicit/rk4-solver.hpp"
#include "ODE/rk/explicit/rk4-variants.hpp"
#include "ODE/multistep/ab-solver.hpp"
#include "ODE/multistep/abm-solver.hpp"
#include "ODE/adaptive/adaptive-common.hpp"
#include "ODE/adaptive/adaptive-custom.hpp"

namespace MathEngine
{ // MathEngine namespace

// Main solver function that selects the appropriate solver based on the parameters
SolverResults solve_ode(const SolverParameters& Params, std::string method = "RK4")
{
    if (method == "RK1")
        return rk1_solver(Params);
    else if (method == "RK2")
        return rk2_solver(Params);
    else if (method == "RK3")
        return rk3_solver(Params);
    else if (method == "RK4")
        return rk4_solver(Params);
    else if (method == "RK4 -3/8")
        return rk4_38_solver(Params);
    else if (method == "RK4 -r")
        return rk4_ralston_solver(Params);
    else if (method == "RK4 -g")
        return rk4_gill_solver(Params);
    else if (method == "AB")
        return adams_bashforth_solver(Params);
    else if (method == "ABM")
        return adams_bashforth_moulton_solver(Params);
    else if (method == "Heun-Euler")
        return heun_euler_solver(Params);
    else if (method == "Bogachi-Shampine")
        return bogacki_shampine_solver(Params);
    else if (method == "Dormand-Prince")
        return dormand_prince_solver(Params);
    else if (method == "RK12")
        return heun_euler_solver(Params);
    else if (method == "RK23")
        return bogacki_shampine_solver(Params);
    else if (method == "RK34")
        return fehlberg34_solver(Params);
    else if (method == "RK45")
        return dormand_prince_solver(Params);
    else if (method == "RK12 -c")
        return rk12_custom_solver(Params);
    else if (method == "RK23 -c")
        return rk23_custom_solver(Params);
    else if (method == "RK34 -c")
        return rk34_custom_solver(Params);
    else if (method == "RK45 -c")
        return rk45_custom_solver(Params);
    else
        throw std::invalid_argument("Unknown solver method: " + method);
}

} // End MathEngine namespace

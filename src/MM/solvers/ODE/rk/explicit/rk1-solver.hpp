#pragma once
#include "../../../../typedefs/header.hpp"

// Euler's method (RK1)
// The simplest explicit method for numerical integration, using a single
// derivative evaluation per step. It's first-order and has limited accuracy.
inline SolverResults rk1_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f  = Params.derivative;
    const auto&  y0 = Params.initial_conditions;
    const double t0 = Params.t0;
    const double t1 = Params.t1;
    const double dt = Params.dt;
    const size_t N  = y0.size();

    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = dMatrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);

    solution[0]    = y0;
    time_points[0] = t0;

    auto  y = y0;
    dVec  k1(N, 0.0);

    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];

        k1 = f(t, y);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt * k1[j];
        
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for rk1_solver
inline dMatrix rk1_solver(
    MyFunc      deriv,
    const dVec& y0,
    double      t0,
    double      t1,
    double      dt
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    return rk1_solver(params).solution;
}

// Function aliases using templates
template<typename... Args>
auto euler_solver(Args&&... args) -> decltype(rk1_solver(std::forward<Args>(args)...))
{
    return rk1_solver(std::forward<Args>(args)...);
}

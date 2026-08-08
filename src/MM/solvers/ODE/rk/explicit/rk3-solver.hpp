#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
// Kutta's third-order method (RK3)
// A standard third-order Runge-Kutta method, offering a balance between accuracy
// and computational cost, better than RK2 but less complex than RK4.
inline SolverResults rk3_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.derivative;
    const auto&  y0   = Params.initialConditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();

    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = dMatrix(num_steps + 1, dVec(N));
    auto         timePoints = dVec(num_steps + 1);

    solution[0]    = y0;
    timePoints[0] = t0;

    auto y        = y0;
    auto y_temp   = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0);
    auto dt_half  = dt * 0.5;
    auto dt_2     = dt * 2.0;
    auto dt_sixth = dt / 6.0;
    auto dt_46    = dt_sixth * 4;

    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        k1 = f(t, y);
        for (size_t j = 0; j < N; ++j) 
            y_temp[j] = y[j] + dt_half * k1[j];
        
        k2 = f(t + dt_half, y_temp);
        for (size_t j = 0; j < N; ++j) 
            y_temp[j] = y[j] - dt * k1[j] + dt_2 * k2[j];
        
        k3 = f(t + dt, y_temp);
        
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k3[j]) + dt_46 * k2[j];
        
        solution[i + 1]    = y;
        timePoints[i + 1] = t + dt;
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.timePoints = timePoints;
    return results;
}

inline SolverResults rk3_solver_callback(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.derivative;
    const auto&  y0   = Params.initialConditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();

    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = dMatrix(num_steps + 1, dVec(N));
    auto         timePoints = dVec(num_steps + 1);

    solution[0]    = y0;
    timePoints[0] = t0;

    auto y        = y0;
    auto y_temp   = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0);
    auto dt_half  = dt * 0.5;
    auto dt_2     = dt * 2.0;
    auto dt_sixth = dt / 6.0;
    auto dt_46    = dt_sixth * 4;

    // Main integration loop
    OneStepSolverResult stepRes;
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        k1 = f(t, y);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k1[j];

        k2 = f(t + dt_half, y_temp);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] - dt * k1[j] + dt_2 * k2[j];

        k3 = f(t + dt, y_temp);

        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k3[j]) + dt_46 * k2[j];

        if (Params.onStep)
        {
            stepRes.sol = y;
            stepRes.timePoint = t;
            stepRes.stepSize = dt;
            Params.onStep(stepRes);
        }

        solution[i + 1]    = y;
        timePoints[i + 1] = t + dt;
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.timePoints = timePoints;
    return results;
}

// Basic interface wrapper for rk3_solver
inline dMatrix rk3_solver(
    MyFunc       deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initialConditions = y0;
    params.t0                = t0;
    params.t1                = t1;
    params.dt                = dt;
    return rk3_solver(params).solution;
}
} // End namespace MathEngine

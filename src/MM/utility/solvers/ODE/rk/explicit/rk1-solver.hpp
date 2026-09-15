#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
// Euler's method (RK1)
// The simplest explicit method for numerical integration, using a single
// derivative evaluation per step. It's first-order and has limited accuracy.
namespace funcWrapper
{
template <bool EnableCallBack>
inline SolverResults rk1(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f  = Params.derivative;
    const auto&  y0 = Params.initialConditions;
    const double t0 = Params.t0;
    const double t1 = Params.t1;
    const double dt = Params.dt;
    const size_t N  = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults{};

    // Initialize solution storage
    const size_t num_steps  = static_cast<size_t>(std::round((t1 - t0) / dt));
    auto         solution   = dMatrix(num_steps + 1, dVec(N));
    auto         timePoints = dVec(num_steps + 1);

    solution[0]   = y0;
    timePoints[0] = t0;

    auto y = y0;
    dVec k1(N, 0.0);

    OneStepSolverResult stepRes;
    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        f(t, y, k1);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt * k1[j];

        solution[i + 1] = y;
        const double t_next = (i+1)==num_steps?t1:t0+(i+1)*dt;
        timePoints[i + 1] = t_next;
        if constexpr (EnableCallBack)
        {
            if (Params.onStep)
            {
                stepRes.sol = y;
                stepRes.timePoint = t0+(i+1)*dt;
                stepRes.stepSize = dt;
                Params.onStep(stepRes);
            }
        }
    }

    auto results       = SolverResults{};
    results.solution   = solution;
    results.timePoints = timePoints;
    return results;
}
} // End funcWrapper namespace

inline SolverResults rk1_solver(const SolverParameters& Params)
{
    return funcWrapper::rk1<false>(Params);
}
inline SolverResults rk1_solver_callback(const SolverParameters& Params)
{
    return funcWrapper::rk1<true>(Params);
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
    params.initialConditions = y0;
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

// inline SolverFunc rk1_wrapper()
// {
//     return [](const SolverParameters& Params) {return rk1_solver(Params);};
// }
} // End namespace MathEngine

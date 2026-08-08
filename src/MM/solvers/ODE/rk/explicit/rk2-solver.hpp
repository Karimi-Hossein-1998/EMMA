#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
namespace funcWrapper
{
template <bool EnableCallBack>
inline SolverResults rk2(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f  = Params.derivative;
    const auto&  y0 = Params.initialConditions;
    const double t0 = Params.t0;
    const double t1 = Params.t1;
    const double dt = Params.dt;
    const size_t N  = y0.size();

    // Initialize solution storage
    const size_t num_steps  = static_cast<size_t>(std::round((t1 - t0) / dt));
    auto         solution   = dMatrix(num_steps + 1, dVec(N));
    auto         timePoints = dVec(num_steps + 1);

    solution[0]   = y0;
    timePoints[0] = t0;

    auto y       = y0;
    auto y_temp  = dVec(N, 0.0);
    dVec k1(N, 0.0), k2(N, 0.0);
    auto dt_half = dt * 0.5;

    OneStepSolverResult stepRes;
    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        f(t, y, k1);

        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k1[j];
        f(t + dt_half, y_temp, k2);

        for (size_t j = 0; j < N; ++j)
            y[j] += dt * k2[j];

        solution[i + 1] = y;
        const double t_next = (i+1)==num_steps?t1:t0+(i+1)*dt;
        timePoints[i + 1] = t_next;
        if constexpr (EnableCallBack)
        {
            if (Params.onStep)
            {
                stepRes.sol = y;
                stepRes.timePoint = timePoints[i+1];
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
} // End namespace funcWrapper

inline SolverResults rk2_solver(const SolverParameters& Params)
{
    return funcWrapper::rk2<false>(Params);
}
inline SolverResults rk2_solver_callback(const SolverParameters& Params)
{
    return funcWrapper::rk2<true>(Params);
}
// Basic interface wrapper for rk2_solver
inline dMatrix rk2_solver(
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
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    return rk2_solver(params).solution;
}

template<typename... Args>
auto midpoint_solver(Args&&... args) -> decltype(rk2_solver(std::forward<Args>(args)...))
{
    return rk2_solver(std::forward<Args>(args)...);
}
} // End namespace MathEngine

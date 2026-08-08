#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace}
// Classical fourth-order Runge-Kutta method (RK4)
// A widely used, well-balanced method for numerical integration, known for its
// accuracy and stability. It requires four derivative evaluations per step.
namespace funcWrapper
{
template <bool EnableCallBack>
inline SolverResults rk4(const SolverParameters& Params)
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

    auto y        = y0;
    auto y_temp   = dVec(N, 0.0);
    dVec k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto dt_half  = dt * 0.5;
    auto dt_third = dt / 3.0;
    auto dt_sixth = dt_third * 0.5;

    // Main integration loop
    OneStepSolverResult stepRes;
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        f(t, y, k1);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k1[j];
        f(t + dt_half, y_temp, k2);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k2[j];
        f(t + dt_half, y_temp, k3);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt * k3[j];
        f(t + dt, y_temp, k4);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k4[j]) + dt_third * (k2[j] + k3[j]);

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

inline SolverResults rk4_solver(const SolverParameters& Params)
{
    return funcWrapper::rk4<false>(Params);
}
inline SolverResults rk4_solver_callback(const SolverParameters& Params)
{
    return funcWrapper::rk4<true>(Params);
}
// Basic interface wrapper for rk4_solver
inline dMatrix rk4_solver(
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
    return rk4_solver(params).solution;
}

template<typename... Args>
auto classical_rk4_solver(Args&&... args) -> decltype(rk4_solver(std::forward<Args>(args)...))
{
    return rk4_solver(std::forward<Args>(args)...);
}
} // End namespace MathEngine

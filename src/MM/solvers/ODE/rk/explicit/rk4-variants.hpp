#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEnggine namespace
// 3/8 Rule RK4
// A variant of the classical RK4 method using a different set of coefficients
// for the final weighted average, trading some accuracy for different stability.
//////////// CALLBACKS ///////////
namespace funcWrapper
{
template <bool EnableCallBack>
inline SolverResults rk4_38(const SolverParameters& Params)
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

    auto y      = y0;
    auto y_temp = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto dt_third = dt / 3.0;
    auto dt_2_3   = dt_third * 2.0;
    auto dt_8     = dt / 8.0;
    auto dt_3_8   = dt_8 * 3.0;

    // Main integration loop
    OneStepSolverResult stepRes;
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        f(t, y, k1);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_third * k1[j];
        f(t + dt_third, y_temp, k2);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] - dt_third * k1[j] + dt * k2[j];
        f(t + dt_2_3, y_temp, k3);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt * (k1[j] - k2[j] + k3[j]);
        f(t + dt, y_temp, k4);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_8 * (k1[j] + k4[j]) + dt_3_8 * (k2[j] + k3[j]);

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

template <bool EnableCallBack>
inline SolverResults rk4_ralston(const SolverParameters& Params)
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

    auto y      = y0;
    auto y_temp = dVec(N, 0.0);
    dVec k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto c2 = 0.4, c3 = 0.45573725, c4 = 1;
    auto a21 = 0.4, a31 = 0.15875964, a32 = 0.29697761, a41 = 0.2181004, a42 = -3.05096516, a43 = 3.83286476;
    auto b1 = 0.17476028, b2 = -0.55148066, b3 = 1.2055356, b4 = 0.17118478;
	auto dtc2 = dt*c2, dtc3 = dt*c3;
    auto dta21 = dt*a21, dta31 = dt*a31, dta32 = dt*a32, dta41 = dt*a41, dta42 = dt*a42, dta43 = dt*a43;
    auto dtb1 = dt*b1, dtb2 = dt*b2, dtb3 = dt*b3, dtb4 = dt*b4;
    // Main integration loop
    OneStepSolverResult stepRes;
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        f(t, y, k1);

        for(size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dta21 * k1[j];
        f(t + dtc2, y_temp, k2);

        for(size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dta31 * k1[j] + dta32 * k2[j];
        f(t + dtc3, y_temp, k3);

        for(size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dta41 * k1[j] + dta42 * k2[j] + dta43 * k3[j];
        f(t + dt, y_temp, k4);

        for (size_t j = 0; j < N; ++j)
            y[j] += dtb1 *  k1[j] + dtb2 * k2[j] + dtb3 * k3[j] + dtb4 * k4[j];

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

    auto results        = SolverResults{};
    results.solution    = solution;
    results.timePoints = timePoints;
    return results;
}

template <bool EnableCallBack>
inline SolverResults rk4_gill(const SolverParameters& Params)
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

    auto y            = y0;
    auto dt_half      = dt * 0.5;
    auto dt_sixth     = dt / 6.0;
    auto sqrt2        = std::sqrt(2.0);
    auto sqrt2inverse = sqrt2 / 2.0;
    auto dt_s2inverse = dt * sqrt2inverse;
    auto a            = sqrt2inverse - 0.5;
    auto dt_a         = dt * a;
    auto b            = 1.0 - sqrt2inverse;
    auto dt_b         = dt * b;
    auto c            = 1.0 + sqrt2inverse;
    auto dt_c         = dt * c;
    auto dy_2         = dt_sixth * b * 2.0;
    auto dy_3         = dt_sixth * c * 2.0;
    auto y_temp       = dVec(N, 0.0);
    dVec k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    // Main integration loop
    OneStepSolverResult stepRes;
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        f(t, y, k1);

        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_half * k1[j];
        f(t + dt_half, y_temp, k2);

        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_a * k1[j] + dt_b * k2[j];
        f(t + dt_half, y_temp, k3);

        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_c * k3[j] - dt_s2inverse * k2[j];
        f(t + dt, y_temp, k4);

        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k4[j]) + dy_2 * k2[j] + dy_3 * k3[j];

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

inline SolverResults rk4_38_solver(const SolverParameters& Params)
{
    return funcWrapper::rk4_38<false>(Params);
}
inline SolverResults rk4_38_solver_callback(const SolverParameters& Params)
{
    return funcWrapper::rk4_38<true>(Params);
}

inline SolverResults rk4_ralston_solver(const SolverParameters& Params)
{
    return funcWrapper::rk4_ralston<false>(Params);
}
inline SolverResults rk4_ralston_solver_callback(const SolverParameters& Params)
{
    return funcWrapper::rk4_ralston<true>(Params);
}

inline SolverResults rk4_gill_solver(const SolverParameters& Params)
{
    return funcWrapper::rk4_gill<false>(Params);
}
inline SolverResults rk4_gill_solver_callback(const SolverParameters& Params)
{
    return funcWrapper::rk4_gill<true>(Params);
}

} // End namespace MathEngine

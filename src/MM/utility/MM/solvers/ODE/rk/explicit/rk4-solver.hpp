#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace}
// Classical fourth-order Runge-Kutta method (RK4)
// A widely used, well-balanced method for numerical integration, known for its
// accuracy and stability. It requires four derivative evaluations per step.
namespace funcWrapper
{
template <bool EnableCallBack, FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> rk4(const ODESolverParameters<Num,MF,CBF>& Params)
{
    // Extract parameters for clarity
    const MF&   f  = Params.derivative;
    const CBF&  cb = Params.onStep;
    const auto& y0 = Params.initialConditions;
    const Num   t0 = Params.t0;
    const Num   t1 = Params.t1;
    const Num   dt = Params.dt;
    const size_t N = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults<Num>{};

    // Initialize solution storage
    const size_t num_steps  = static_cast<size_t>(std::round((t1 - t0) / dt));
    Matrix<Num> solution(num_steps + 1, N);
    Vec<Num>    timePoints(num_steps + 1);

    solution.SetRow(0,y0);
    timePoints[0] = t0;

    Vec<Num> y(y0), y_temp(N,Num{}), k1(N, Num{}), k2(N, Num{}), k3(N, Num{}), k4(N, Num{});
    Num dt_half  = dt * 0.5;
    Num dt_third = dt / 3.0;
    Num dt_sixth = dt_third * 0.5;

    // Main integration loop
    OneStepSolverResult<Num> stepRes{
        .sol = y0,
        .timePoint = t0,
        .stepSize = dt
    };
    if constexpr (EnableCallBack)
    {
        if constexpr (!std::is_same_v<CBF,DefaultCallBackFunc>)
        {
            cb(stepRes);
        }
    }
    for (size_t i = 0; i < num_steps; ++i)
    {
        const Num t = timePoints[i];

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

        solution.SetRow(i+1, y);
        const Num t_next = (i+1)==num_steps?t1:t0+(i+1)*dt;
        timePoints[i + 1] = t_next;
        if constexpr (EnableCallBack)
        {
            if constexpr (!std::is_same_v<CBF,DefaultCallBackFunc>)
            {
                stepRes.sol       = y;
                stepRes.timePoint = t0+(i+1)*dt;
                stepRes.stepSize  = dt;
                cb(stepRes);
            }
        }
    }

    SolverResults<Num> results{
        .solution   = solution,
        .timePoints = timePoints
    };
    return results;
}
} // End namespace funcWrapper

template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> rk4_solver(const ODESolverParameters<Num,MF,CBF>& Params)
{
    return funcWrapper::rk4<false>(Params);
}
template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> rk4_solver_callback(const ODESolverParameters<Num,MF,CBF>& Params)
{
    return funcWrapper::rk4<true>(Params);
}
// Basic interface wrapper for rk4_solver
template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline Matrix<Num> rk4_solver(
    MF    deriv,
    const Vec<Num>& y0,
    Num   t0,
    Num   t1,
    Num   dt
)
{
    ODESolverParameters<Num,MF> params{
        .derivative        = deriv,
        .initialConditions = y0,
        .t0                = t0,
        .t1                = t1,
        .dt                = dt
    };
    return rk4_solver(params).solution;
}

template<typename... Args>
auto classical_rk4_solver(Args&&... args) -> decltype(rk4_solver(std::forward<Args>(args)...))
{
    return rk4_solver(std::forward<Args>(args)...);
}
} // End namespace MathEngine

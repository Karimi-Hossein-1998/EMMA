#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
namespace funcWrapper
{

template <bool EnableCallBack, FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> rk2(const ODESolverParameters<Num,MF,CBF>& Params)
{
    // Extract parameters for clarity
    const MF&    f  = Params.derivative;
    const CBF&   cb = Params.onStep;
    const auto&  y0 = Params.initialConditions;
    const Num    t0 = Params.t0;
    const Num    t1 = Params.t1;
    const Num    dt = Params.dt;
    const size_t N  = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults<Num>{};

    // Initialize solution storage
    const size_t num_steps = static_cast<size_t>(std::round((t1 - t0) / dt));
    Matrix<Num> solution   = Matrix<Num>(num_steps + 1,N);
    Vec<Num>    timePoints = Vec<Num>(num_steps + 1);

    solution.SetRow(0, y0);
    timePoints[0] = t0;

    Vec<Num> y(y0), y_temp(N,Num{}),k1(N, Num{}), k2(N, Num{});
    Num dt_half = static_cast<Num>(dt * 0.5);

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
            y[j] += dt * k2[j];

        solution[i + 1] = y;
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
inline SolverResults<Num> rk2_solver(const ODESolverParameters<Num,MF,CBF>& Params)
{
    return funcWrapper::rk2<false>(Params);
}

template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> rk2_solver_callback(const ODESolverParameters<Num,MF,CBF>& Params)
{
    return funcWrapper::rk2<true>(Params);
}
// Basic interface wrapper for rk2_solver
template <FPNumber Num, MyFunc<Num> MF>
inline Matrix<Num> rk2_solver(
    MF    deriv,
    const Vec<Num>& y0,
    Num   t0,
    Num   t1,
    Num   dt
)
{
    ODESolverParameters<Num, MF> params{
        .derivative        = deriv,
        .initialConditions = y0,
        .t0                = t0,
        .t1                = t1,
        .dt                = dt
    };
    return rk2_solver(params).solution;
}

template<typename... Args>
auto midpoint_solver(Args&&... args) -> decltype(rk2_solver(std::forward<Args>(args)...))
{
    return rk2_solver(std::forward<Args>(args)...);
}
} // End namespace MathEngine

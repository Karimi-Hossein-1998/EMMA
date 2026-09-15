#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../namespaces/abm-coefs.hpp"
#include "../rk/explicit/rk4-solver.hpp"

namespace MathEngine
{ // namespace MathEngine
namespace funcWrapper
{ // namespace funcWrapper
// General Adams-Bashforth-Moulton predictor-corrector solver (order 1-10)
template <bool EnableCallBack, FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> adams_bashforth_moulton(const ODESolverParameters<Num,MF,CBF>& Params)
{
    // Extract parameters for clarity
    const MF&    f               = Params.derivative;
    const CBF&   cb              = Params.onStep;
    const auto&  y0              = Params.initialConditions;
    const Num    t0              = Params.t0;
    const Num    t1              = Params.t1;
    const Num    dt              = Params.dt;
    const int    order           = Params.order;
    const int    corrector_iters = Params.iterations;
    const size_t N               = y0.size();
    const int    max_order       = 10;

    if (order < 1 || order > max_order)
        throw std::invalid_argument("ABM order must be between 1 and 10");

    // Initialize solution and time storage
    const size_t num_steps = static_cast<size_t>((t1 - t0) / dt);
    Matrix<Num> solution(num_steps + 1, N);
    Vec<Num>    timePoints(num_steps + 1, Num{});
    Vec<Num>    y(y0),y_temp(N,Num{});
    solution.SetRow(0, y0);
    timePoints[0] = t0;

    // size_t boot_strap{std::min(static_cast<size_t>(order-1),num_steps)};
    Num t1rk = t0 + (order-1)*dt;
    ODESolverParameters<Num,MF,CBF> rk_params{
        .derivative        = f,
        .onStep            = cb,
        .initialConditions = y0,
        .t0                = t0,
        .t1                = t1rk,
        .dt                = dt
    };
    auto res_rk = rk4<EnableCallBack>(rk_params);    
    // Bootstrap with RK4 for the first (order-1) steps
    for (size_t i = 1; i < order; ++i)
    {
        solution.SetRow(i, res_rk.solution[i]);
        timePoints[i] = t0+i*dt;
    }

    // Store derivative history
    Matrix<Num> f_hist(order, N);
    for (int i = 0; i < order; ++i)
    {
        y = std::move(Vec<Num>(std::from_range,solution[i]));
        f(timePoints[i], y, y_temp);
        f_hist.SetRow(i,y_temp);
    }
    // Main ABM loop
    Matrix<Num> f_corr(order + 1, N);
    OneStepSolverResult<Num> oneStep{
        .sol       = y,
        .timePoint = timePoints[std::min(static_cast<size_t>(order),num_steps)],
        .stepSize  = dt
    };
    for (size_t i = order-1; i < num_steps; ++i)
    {
        // Predictor (Adams-Bashforth)
        solution.SetRow(i + 1,solution[i]);
        for (size_t j = 0; j < N; ++j)
        {
            Num increment = Num{};
            for (int k = 0; k < order; ++k)
                increment += abm_coefs::ab_coefs[order - 1][k] * f_hist[k,j];
            solution[i + 1,j] += dt * increment;
        }

        // Corrector (Adams-Moulton)
        for (int k = 0; k < order; ++k)
            f_corr.SetRow(k,f_hist[k]);
        for (int iter = 0; iter < corrector_iters; ++iter)
        {
            y = std::move(Vec<Num>(std::from_range,solution[i+1]));
            f(timePoints[i+1],y,y_temp);
            f_corr.SetRow(order, y_temp);
            for (size_t j = 0; j < N; ++j)
            {
                Num increment = Num{};
                for (int k = 0; k < order + 1; ++k)
                    increment += abm_coefs::am_coefs[order - 1][k] * f_corr[k,j];
                solution[i + 1,j] = solution[i,j] + dt * increment;
            }
        }
        y = std::move(Vec<Num>(std::from_range,solution[i+1]));
        timePoints[i + 1] = t0 + (i+1)*dt;
        if constexpr (EnableCallBack)
        {
            if constexpr (!std::is_same_v<CBF,DefaultCallBackFunc>)
            {
                oneStep.sol       = y;
                oneStep.timePoint = timePoints[i+1];
                oneStep.stepSize  = dt;
                cb(oneStep);
            }
        }
        // Update history by shifting and adding the new derivative
        for (int k = 0; k < order - 1; ++k)
            f_hist.SetRow(k,f_hist[k + 1]);
        f(timePoints[i + 1], y,y_temp);
        f_hist.SetRow(order - 1,y_temp);
    }

    SolverResults<Num> results{
        .solution   = solution,
        .timePoints = timePoints
    };
    return results;
}
} // End funcWrapper namespace
template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> adams_bashforth_moulton_solver(const ODESolverParameters<Num,MF,CBF>& params)
{
    return funcWrapper::adams_bashforth_moulton<false>(params);
}
template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline SolverResults<Num> adams_bashforth_moulton_solver_callback(const ODESolverParameters<Num,MF,CBF>& params)
{
    return funcWrapper::adams_bashforth_moulton<true>(params);
}
// Basic interface wrapper for Adams-Bashforth-Moulton solver
template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF>
inline Matrix<Num> adams_bashforth_moulton_solver(
    MF              deriv,
    const Vec<Num>& y0,
    Num             t0,
    Num             t1,
    Num             dt,
    int             order,
    int             corrector_iters
)
{
    ODESolverParameters<Num,MF,CBF> params{
        .derivative        = deriv,
        .initialConditions = y0,
        .t0                = t0,
        .t1                = t1,
        .dt                = dt,
        .order             = order,
        .iterations        = corrector_iters
    };
    return adams_bashforth_moulton_solver_wrapper<false>(params).solution;
}
} // End MathEngine namespace

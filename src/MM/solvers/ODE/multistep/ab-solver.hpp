#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../namespaces/abm-coefs.hpp"
#include "../rk/explicit/rk4-solver.hpp"

namespace MathEngine
{ // namespace MathEngine
namespace funcWrapper
{ // namespace funcWrapper
// General Adams-Bashforth-Moulton predictor-corrector solver (order 1-10)
template <bool EnableCallBack>
inline SolverResults adams_bashforth(const ODESolverParameters& Params)
{
    // Extract parameters for clarity
    const MyFunc&       f         = Params.derivative;
    const CallBackFunc& cb        = Params.onStep;
    const auto&         y0        = Params.initialConditions;
    const double        t0        = Params.t0;
    const double        t1        = Params.t1;
    const double        dt        = Params.dt;
    const int           order     = Params.order;
    const size_t        N         = y0.size();
    const int           max_order = 10;

    if (order < 1 || order > max_order)
        throw std::invalid_argument("ABM order must be between 1 and 10");

    // Initialize solution and time storage
    const size_t num_steps = static_cast<size_t>((t1 - t0) / dt);
    Matrix<double> solution(num_steps + 1, N);
    Vec<double>    timePoints(num_steps + 1, double{});
    Vec<double>    y(y0),y_temp(N,double{});

    double*       __restrict__ yptr   = y.data();
    double*       __restrict__ solptr = solution.ptr();
    const double* __restrict__ y0ptr  = y0.data();
    #pragma omp simd
    for (size_t i = 0; i < N; ++i)
        solptr[i] = y0ptr[i];
    // solution.SetRow(0, y0);
    timePoints[0] = t0;
    double t1rk = t0 + (order-1)*dt;
    ODESolverParameters rk_params{
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
        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            solptr[i*N+j] = res_rk.solution[i,j];
        // solution.SetRow(i, res_rk.solution[i]);
        timePoints[i + 1] = t0+(i+1)*dt;
    }

    // Store derivative history
    Matrix<double> f_hist(order, N);
    for (int i = 0; i < order; ++i)
    {
        y = std::move(Vec<double>(std::from_range,solution[i]));
        f(timePoints[i], y, y_temp);
        f_hist.SetRow(i,y_temp);
    }
    // Main ABM loop
    Matrix<double> f_corr(order + 1, N);
    OneStepSolverResult oneStep{
        .sol       = y,
        .timePoint = timePoints[std::min(static_cast<size_t>(order),num_steps)],
        .stepSize  = dt
    };
    for (size_t i = order - 1; i < num_steps; ++i)
    {
        // Predictor (Adams-Bashforth)
        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            solptr[(i+1)*N+j] = solptr[i*N+j];
        // solution.SetRow(i + 1,solution[i]);
        for (size_t j = 0; j < N; ++j)
        {
            double increment = double{};
            for (int k = 0; k < order; ++k)
                increment += abm_coefs::ab_coefs[order - 1][k] * f_hist[k,j];
            solution[i + 1,j] += dt * increment;
        }
        timePoints[i + 1] = t0 + (i+1)*dt;

        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            yptr[j] = solptr[(i+1)*N+j];
        // y = std::move(Vec<double>(std::from_range,solution[i+1]));
        if constexpr (EnableCallBack)
        {
            oneStep.sol       = y;
            oneStep.timePoint = timePoints[i+1];
            oneStep.stepSize  = dt;
            cb(oneStep);
        }
        // Update history by shifting and adding the new derivative
        for (int k = 0; k < order - 1; ++k)
            f_hist.SetRow(k,f_hist[k + 1]);
        f(timePoints[i + 1], y,y_temp);
        f_hist.SetRow(order - 1,y_temp);
    }

    SolverResults results{
        .solution   = solution,
        .timePoints = timePoints
    };
    return results;
}
} // End funcWrapper namespace
inline SolverResults adams_bashforth_solver(const ODESolverParameters& params)
{
    return funcWrapper::adams_bashforth<false>(params);
}
inline SolverResults adams_bashforth_solver_callback(const ODESolverParameters& params)
{
    return funcWrapper::adams_bashforth<true>(params);
}
// Basic interface wrapper for Adams-Bashforth-Moulton solver
inline Matrix<double> adams_bashforth_solver(
    MyFunc             deriv,
    const Vec<double>& y0,
    double             t0,
    double             t1,
    double             dt,
    uint8_t            order
)
{
    ODESolverParameters params{
        .derivative        = deriv,
        .initialConditions = y0,
        .t0                = t0,
        .t1                = t1,
        .dt                = dt,
        .order             = order,
    };
    return funcWrapper::adams_bashforth<false>(params).solution;
}
} // End MathEngine namespace


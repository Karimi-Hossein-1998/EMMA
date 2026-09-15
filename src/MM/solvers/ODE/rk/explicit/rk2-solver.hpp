#pragma once
#include "../../../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
namespace funcWrapper
{

template <bool EnableCallBack>
inline SolverResults rk2(const ODESolverParameters& Params)
{
    // Extract parameters for clarity
    const MyFunc&       f  = Params.derivative;
    const CallBackFunc& cb = Params.onStep;
    const Vec<double>&  y0 = Params.initialConditions;
    const double        t0 = Params.t0;
    const double        t1 = Params.t1;
    const double        dt = Params.dt;
    const size_t N  = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults{};

    // Initialize solution storage
    const size_t num_steps = static_cast<size_t>(std::round((t1 - t0) / dt));
    Matrix<double> solution   = Matrix<double>(num_steps + 1,N);
    Vec<double>    timePoints = Vec<double>(num_steps + 1);

    Vec<double> y(y0), y_temp(N,double{}),k1(N, double{}), k2(N, double{});
    double dt_half = static_cast<double>(dt * 0.5);

    // double*       __restrict__ y_tptr = y_temp.data();
    // double*       __restrict__ k1ptr  = k1.data();
    // double*       __restrict__ k2ptr  = k2.data();
    // double*       __restrict__ k3ptr  = k3.data();
    // double*       __restrict__ k4ptr  = k4.data();
    double*       __restrict__ yptr   = y.data();
    double*       __restrict__ solptr = solution.ptr();
    const double* __restrict__ y0ptr  = y0.data();
    #pragma omp simd
    for (size_t i = 0; i < N; ++i)
        solptr[i] = y0ptr[i];
    // solution.SetRow(0, y0);
    timePoints[0] = t0;
    // Main integration loop
    OneStepSolverResult stepRes{
        .sol = y0,
        .timePoint = t0,
        .stepSize = dt
    };
    if constexpr (EnableCallBack)
    {
        cb(stepRes);
    }
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = timePoints[i];

        f(t, y, k1);

        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k1[j];
        f(t + dt_half, y_temp, k2);

        for (size_t j = 0; j < N; ++j)
            y[j] += dt * k2[j];

        #pragma omp simd
        for (size_t j = 0; j<N; ++j)
            solptr[(i+1)*N+j] = y[j];
        // solution.SetRow(i+1,y);
        const double t_next = (i+1)==num_steps?t1:t0+(i+1)*dt;
        timePoints[i + 1] = t_next;
        if constexpr (EnableCallBack)
        {
            stepRes.sol       = y;
            stepRes.timePoint = t_next;
            stepRes.stepSize  = dt;
            cb(stepRes);
        }
    }

    SolverResults results{
        .solution   = solution,
        .timePoints = timePoints
    };
    return results;
}
} // End namespace funcWrapper

inline SolverResults rk2_solver(const ODESolverParameters& Params)
{
    return funcWrapper::rk2<false>(Params);
}

inline SolverResults rk2_solver_callback(const ODESolverParameters& Params)
{
    return funcWrapper::rk2<true>(Params);
}
// Basic interface wrapper for rk2_solver
inline Matrix<double> rk2_solver(
    MyFunc             deriv,
    const Vec<double>& y0,
    double             t0,
    double             t1,
    double             dt
)
{
    ODESolverParameters params{
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

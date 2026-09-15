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
inline SolverResults rk4(const ODESolverParameters& Params)
{
    // Extract parameters for clarity
    const MyFunc&       f  = Params.derivative;
    const CallBackFunc& cb = Params.onStep;
    const auto&         y0 = Params.initialConditions;
    const double        t0 = Params.t0;
    const double        t1 = Params.t1;
    const double        dt = Params.dt;
    const size_t        N = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults{};

    // Initialize solution storage
    const size_t num_steps  = static_cast<size_t>(std::round((t1 - t0) / dt));
    Matrix<double> solution(num_steps + 1, N);
    Vec<double>    timePoints(num_steps + 1);

    Vec<double> y(y0), y_temp(N,double{}), k1(N, double{}), k2(N, double{}), k3(N, double{}), k4(N, double{});
    double dt_half  = dt * 0.5;
    double dt_third = dt / 3.0;
    double dt_sixth = dt_third * 0.5;

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
        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k1[j];
        f(t + dt_half, y_temp, k2);
        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k2[j];
        f(t + dt_half, y_temp, k3);
        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt * k3[j];
        f(t + dt, y_temp, k4);
        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k4[j]) + dt_third * (k2[j] + k3[j]);

        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            solptr[(i+1)*N+j] = y[j];
        // solution.SetRow(i+1, y);
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

inline SolverResults rk4_solver(const ODESolverParameters& Params)
{
    return funcWrapper::rk4<false>(Params);
}
inline SolverResults rk4_solver_callback(const ODESolverParameters& Params)
{
    return funcWrapper::rk4<true>(Params);
}
// Basic interface wrapper for rk4_solver
inline Matrix<double> rk4_solver(
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
    return rk4_solver(params).solution;
}

template<typename... Args>
auto classical_rk4_solver(Args&&... args) -> decltype(rk4_solver(std::forward<Args>(args)...))
{
    return rk4_solver(std::forward<Args>(args)...);
}
} // End namespace MathEngine

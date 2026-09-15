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
inline SolverResults rk4_38(const ODESolverParameters& Params)
{
    // Extract parameters for clarity
    const MyFunc&        f  = Params.derivative;
    const CallBackFunc& cb = Params.onStep;
    const auto&         y0 = Params.initialConditions;
    const double        t0 = Params.t0;
    const double        t1 = Params.t1;
    const double        dt = Params.dt;
    const size_t        N  = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults{};

    // Initialize solution storage
    const size_t num_steps = static_cast<size_t>(std::round((t1 - t0) / dt));
    Matrix<double> solution(num_steps + 1, N);
    Vec<double> timePoints(num_steps + 1);

    Vec<double>  y(y0), y_temp(N,0.0), k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    double dt_third = dt / 3.0;
    double dt_2_3   = dt_third * 2.0;
    double dt_8     = dt / 8.0;
    double dt_3_8   = dt_8 * 3.0;

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

        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            solptr[(i+1)*N+j] = y[j];
        // solution.SetRow(i + 1, y);
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

    SolverResults  results{
        .solution   = solution,
        .timePoints = timePoints
    };
    return results;
}

template <bool EnableCallBack>
inline SolverResults rk4_ralston(const ODESolverParameters& Params)
{
    // Extract parameters for clarity
    const MyFunc&       f  = Params.derivative;
    const CallBackFunc& cb = Params.onStep;
    const auto&         y0 = Params.initialConditions;
    const double        t0 = Params.t0;
    const double        t1 = Params.t1;
    const double        dt = Params.dt;
    const size_t        N  = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults{};

    // Initialize solution storage
    const size_t num_steps = static_cast<size_t>(std::round((t1 - t0) / dt));
    Matrix<double> solution(num_steps + 1, N);
    Vec<double> timePoints(num_steps + 1);

    Vec<double>  y(y0), y_temp(N,0.0), k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    double c2 = 0.4, c3 = 0.45573725, c4 = 1.0;
    double a21 = 0.4, a31 = 0.15875964, a32 = 0.29697761, a41 = 0.2181004, a42 = -3.05096516, a43 = 3.83286476;
    double b1 = 0.17476028, b2 = -0.55148066, b3 = 1.2055356, b4 = 0.17118478;
	double dtc2 = dt*c2, dtc3 = dt*c3;
    double dta21 = dt*a21, dta31 = dt*a31, dta32 = dt*a32, dta41 = dt*a41, dta42 = dt*a42, dta43 = dt*a43;
    double dtb1 = dt*b1, dtb2 = dt*b2, dtb3 = dt*b3, dtb4 = dt*b4;
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

        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            solptr[(i+1)*N+j] = y[j];
        // solution.SetRow(i + 1, y);
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

template <bool EnableCallBack>
inline SolverResults rk4_gill(const ODESolverParameters& Params)
{
    // Extract parameters for clarity
    const MyFunc&       f  = Params.derivative;
    const CallBackFunc& cb = Params.onStep;
    const auto&         y0 = Params.initialConditions;
    const double        t0 = Params.t0;
    const double        t1 = Params.t1;
    const double        dt = Params.dt;
    const size_t        N  = y0.size();
    if (N==0 || dt<=1e-13 || t0>t1) return SolverResults{};

    // Initialize solution storage
    const size_t num_steps = static_cast<size_t>(std::round((t1 - t0) / dt));
    Matrix<double> solution(num_steps + 1, N);
    Vec<double> timePoints(num_steps + 1);

    Vec<double> y(y0), y_temp(N,0.0), k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    double dt_half      = dt * 0.5;
    double dt_sixth     = dt / 6.0;
    double sqrt2        = std::sqrt(2.0);
    double sqrt2inverse = sqrt2 / 2.0;
    double dt_s2inverse = dt * sqrt2inverse;
    double a            = sqrt2inverse - 0.5;
    double dt_a         = dt * a;
    double b            = 1.0 - sqrt2inverse;
    double dt_b         = dt * b;
    double c            = 1.0 + sqrt2inverse;
    double dt_c         = dt * c;
    double dy_2         = dt_sixth * b * 2.0;
    double dy_3         = dt_sixth * c * 2.0;
    // Main integration loop
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

        #pragma omp simd
        for (size_t j = 0; j < N; ++j)
            solptr[(i+1)*N+j] = y[j];
        // solution.SetRow(i + 1, y);
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

inline SolverResults rk4_38_solver(const ODESolverParameters& Params)
{
    return funcWrapper::rk4_38<false>(Params);
}
inline SolverResults rk4_38_solver_callback(const ODESolverParameters& Params)
{
    return funcWrapper::rk4_38<true>(Params);
}

inline SolverResults rk4_ralston_solver(const ODESolverParameters& Params)
{
    return funcWrapper::rk4_ralston<false>(Params);
}
inline SolverResults rk4_ralston_solver_callback(const ODESolverParameters& Params)
{
    return funcWrapper::rk4_ralston<true>(Params);
}

inline SolverResults rk4_gill_solver(const ODESolverParameters& Params)
{
    return funcWrapper::rk4_gill<false>(Params);
}
inline SolverResults rk4_gill_solver_callback(const ODESolverParameters& Params)
{
    return funcWrapper::rk4_gill<true>(Params);
}

} // End namespace MathEngine

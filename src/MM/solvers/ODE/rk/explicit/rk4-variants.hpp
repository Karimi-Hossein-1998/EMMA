#pragma once
#include "../../../../typedefs/header.hpp"

// 3/8 Rule RK4
// A variant of the classical RK4 method using a different set of coefficients
// for the final weighted average, trading some accuracy for different stability.
inline SolverResults rk4_38_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.derivative;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();

    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = dMatrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);

    solution[0]    = y0;
    time_points[0] = t0;

    auto y            = y0;
    auto y_temp       = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto dt_third     = dt / 3.0;
    auto dt_2_3       = dt * 2.0 / 3.0;
    auto dt_8         = dt / 8.0;
    auto dt_3_8       = dt * 3.0 / 8.0;

    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];

        k1 = f(t, y);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_third * k1[j];
        k2 = f(t + dt_third, y_temp);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] - dt_third * k1[j] + dt * k2[j];
        k3 = f(t + dt_2_3, y_temp);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt * (k1[j] - k2[j] + k3[j]);
        k4 = f(t + dt, y_temp);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_8 * (k1[j] + k4[j]) + dt_3_8 * (k2[j] + k3[j]);
        
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for 3/8 Rule RK4
inline dMatrix rk4_38_solver(
    MyFunc       deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    return rk4_38_solver(params).solution;
}

// Ralston's RK4 Method
// An RK4 variant optimized for a smaller error bound with a specific set of
// coefficients for the stage evaluations and the final update.
inline SolverResults rk4_ralston_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.derivative;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();

    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = dMatrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);

    solution[0]    = y0;
    time_points[0] = t0;

    auto y            = y0;
    auto y_temp       = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto dt_2_5       = dt * 2.0 / 5.0;
    auto sqrt_5       = std::sqrt(5.0);
    auto dt_1_4       = dt * 0.25;
    auto dt_1_40      = dt_1_4 * 0.1;
    auto k3_t         = dt_1_4 * ( sqrt_5 - 1.0 );
    auto k3_1         = dt_1_40 * ( 11.0 - sqrt_5 );
    auto k3_2         = dt_1_40 * ( 25.0 - 13.0 * sqrt_5 );
    auto k4_1         = dt_1_40 * ( 10.0 - 7.0 * sqrt_5 );
    auto k4_3         = dt_1_4 * ( 5.0 + sqrt_5 );
    auto dt_1_20      = dt_1_4 * 0.2;
    auto dy_1         = dt_1_20 * ( 5.0 - sqrt_5 );
    auto dy_3         = dt_1_20 * ( 5.0 + 3.0 * sqrt_5 );

    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        
        k1 = f(t, y);
        
        for(size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_2_5 * k1[j];
        k2 = f(t + dt_2_5, y_temp);

        for(size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + k3_1 * k1[j] + k3_2 * k2[j];
        k3 = f(t + k3_t, y_temp);

        for(size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + k4_1 * ( k1[j] + k2[j]) + k4_3 * k3[j];
        k4 = f(t + dt, y_temp);

        for (size_t j = 0; j < N; ++j)
            y[j] += dy_1 * ( k1[j] + k2[j] + k4[j]) + dy_3 * k3[j];

        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for Ralston's RK4
inline dMatrix rk4_ralston_solver(
    MyFunc       deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    return rk4_ralston_solver(params).solution;
}

// Gill's RK4 Method
// A memory-efficient RK4 variant that reduces the number of temporary
// variables required by reusing storage during the stage calculations.
inline SolverResults rk4_gill_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.derivative;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();

    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = dMatrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);

    solution[0]    = y0;
    time_points[0] = t0;

    auto y            = y0;
    auto dt_half      = dt * 0.5;
    auto dt_sixth     = dt / 6.0;
    auto sqrt2        = std::sqrt(2.0);
    auto sqrt2inverse = sqrt2 / 2.0;
    auto dt_s2inverse = dt * sqrt2inverse;
    auto a            = 1.0 - sqrt2inverse;
    auto dt_a_half    = dt_half * a;
    auto b            = 1.0 + sqrt2inverse;
    auto dt_b_half    = dt_half * b;
    auto dt_b         = dt * b;
    auto dy_2         = dt_sixth * a * 2.0;
    auto dy_3         = dt_sixth * b * 2.0;
    auto y_temp       = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        
        k1 = f(t, y);
        
        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_half * k1[j];
        k2 = f(t + dt_half, y_temp);

        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_a_half * k1[j] + dt_b_half * k2[j];
        k3 = f(t + dt_half, y_temp);

        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_b * k3[j] - dt_s2inverse * k2[j];
        k4 = f(t + dt, y_temp);

        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k4[j]) + dy_2 * k2[j] + dy_3 * k3[j];

        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for Gill's RK4
inline dMatrix rk4_gill_solver(
    MyFunc       deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    return rk4_gill_solver(params).solution;
}
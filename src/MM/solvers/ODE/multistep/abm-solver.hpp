#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../namespaces/abm-coefs.hpp"
#include "../rk/explicit/rk4-solver.hpp"

// General Adams-Bashforth-Moulton predictor-corrector solver (order 1-10)
inline SolverResults adams_bashforth_moulton_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f       = Params.derivative;
    const auto&  y0      = Params.initial_conditions;
    const double t0      = Params.t0;
    const double t1      = Params.t1;
    const double dt      = Params.dt;
    const int    order   = Params.order;
    const int    corrector_iters = Params.iterations;
    const size_t N       = y0.size();
    const int    max_order = 10;

    if (order < 1 || order > max_order)
        throw std::invalid_argument("ABM order must be between 1 and 10");

    // Initialize solution and time storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = dMatrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1, 0.0);

    solution[0]    = y0;
    time_points[0] = t0;

    // Bootstrap with RK4 for the first (order-1) steps
    for (size_t i = 0; i < static_cast<size_t>(order - 1) && i < num_steps; ++i) {
        SolverParameters rk_params;
        rk_params.derivative         = f;
        rk_params.initial_conditions = solution[i];
        rk_params.t0                 = time_points[i];
        rk_params.t1                 = time_points[i] + dt;
        rk_params.dt                 = dt;
        solution[i + 1]    = rk4_solver(rk_params).solution.back();
        time_points[i + 1] = time_points[i] + dt;
    }

    // Store derivative history
    auto f_hist = dMatrix(order, dVec(N));
    for (int i = 0; i < order; ++i)
        f_hist[i] = f(time_points[i], solution[i]);

    // Main ABM loop
    auto f_corr = dMatrix(order + 1, dVec(N));
    for (size_t i = order - 1; i < num_steps; ++i)
    {
        // Predictor (Adams-Bashforth)
        solution[i + 1] = solution[i];
        for (size_t j = 0; j < N; ++j)
        {
            double increment = 0.0;
            for (int k = 0; k < order; ++k)
                increment += abm_coefs::ab_coefs[order - 1][k] * f_hist[k][j];
            solution[i + 1][j] += dt * increment;
        }
        time_points[i + 1] = time_points[i] + dt;

        // Corrector (Adams-Moulton)
        for (int k = 0; k < order; ++k)
            f_corr[k] = f_hist[k];
        for (int iter = 0; iter < corrector_iters; ++iter)
        {
            f_corr[order] = f(time_points[i + 1], solution[i + 1]);
            for (size_t j = 0; j < N; ++j)
            {
                double increment = 0.0;
                for (int k = 0; k < order + 1; ++k)
                    increment += abm_coefs::am_coefs[order - 1][k] * f_corr[k][j];
                solution[i + 1][j] = solution[i][j] + dt * increment;
            }
        }

        // Update history by shifting and adding the new derivative
        for (int k = 0; k < order - 1; ++k)
            f_hist[k] = f_hist[k + 1];
        f_hist[order - 1] = f(time_points[i + 1], solution[i + 1]);
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for Adams-Bashforth-Moulton solver
inline dMatrix adams_bashforth_moulton_solver(
    MyFunc      deriv,
    const dVec& y0,
    double      t0,
    double      t1,
    double      dt,
    int         order,
    int         corrector_iters
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    params.order              = order;
    params.iterations         = corrector_iters;
    return adams_bashforth_moulton_solver(params).solution;
} 
#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../namespaces/abm-coefs.hpp"
#include "../../ODE/multistep/ab-solver.hpp"
#include "../../ODE/multistep/abm-solver.hpp"
#include "../../ODE/rk/explicit/rk4-solver.hpp"

// Linear interpolation to retrieve solution value at a delayed time t_delay
// Solution history is stored at discrete time points t_0, t_1, ..., t_n
// If t_delay is before stored history, uses history_f(t_delay)
inline dVec interpolate_delay(
    const dMatrix&        solution,
    const dVec&          time_points,
    const double         t_delay,
    size_t               current_index,
    const Func&          history_f,
    const double&        t0
)
{
    size_t N = solution[0].size();
    if (t_delay <= t0)
    {
        return history_f(t_delay);
    }

    // Find the two nearest time points for interpolation
    int idx_lower = -1;
    for (int i = current_index; i >= 0; --i)
    {
        if (time_points[i] <= t_delay)
        {
            idx_lower = i;
            // std::cout << "Lowwer index is: " << idx_lower << " | time = " << time_points[i] << '\n';
            break;
        }
    }

    if (idx_lower < 0)
        return history_f(t_delay);

    if (idx_lower == current_index)
        return solution[current_index]; // Can't interpolate beyond current solution

    // Linear interpolation between idx_lower and idx_lower+1
    const double t_lower  = time_points[idx_lower];
    const double t_upper  = time_points[idx_lower + 1];
    const double dt_local = t_upper - t_lower;
    const double alpha    = (t_delay - t_lower) / dt_local;
    if (std::fabs(dt_local) < 1.0e-14 || std::fabs(alpha) < 1.0e-14)
        return solution[idx_lower];

    dVec result(N);
    for (size_t j = 0; j < N; ++j)
    {
        result[j] = (1.0 - alpha) * solution[idx_lower][j] + alpha * solution[idx_lower + 1][j];
    }

    return result;
}

// General explicit Adams-Bashforth solver (order 1-10)
// A multi-step method that uses a history of derivative evaluations to improve
// efficiency. Requires a single-step method like RK4 for bootstrapping.
inline SolverResults adams_bashforth_solver_dde(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  y0        = Params.initial_conditions;
    const double t0        = Params.t0;
    const double t1        = Params.t1;
    const double user_dt   = Params.dt;
    const int    order     = Params.order;
    const auto   delays    = Params.delay_times;
    const size_t N         = y0.size();
    const int    max_order = 10;
    const size_t num_delay = delays.size();
    const auto   history_f = Params.history_function;
    const double tau_max = (num_delay > 0) ? *std::max_element(delays.begin(), delays.end()) : 0.0;
    const double tau_min = (num_delay > 0) ? *std::min_element(delays.begin(), delays.end()) : 0.0;

    if (tau_max <= 1.0e-6 * user_dt)
    {
        std::cout << "Delays are negligible; Max{delay} = " << tau_max << std::endl << "Using ODE solver." << std::endl;
        return adams_bashforth_solver(Params); // If delays are negligible, use the ODE solver
    }
    if (num_delay == 0)
    {
        std::cout << "Delays not provided! Using ODE solver.\n";
        return adams_bashforth_solver(Params); // If no delays, use the ODE solver
    }

    if (order < 1 || order > max_order)
        throw std::invalid_argument("AB order must be between 1 and 10");

    // Setting the RHS
    DelayFunc delay_f        = nullptr;
    DelayFuncMulti delay_f_m = nullptr;
    if (num_delay == 1)
        delay_f = Params.delay_derivative;
    else
        delay_f_m = Params.delay_derivative_multi;

    if (delay_f == nullptr && delay_f_m == nullptr)
        throw std::invalid_argument("DDE derivative function must be provided for non-zero delays");
    
    // Initialize solution storage
    // Determine step size to capture delay history adequately
    const size_t num_per_steps = std::max(static_cast<size_t>(8), static_cast<size_t>(std::ceil(1.0 / user_dt)));
    const double min_dt        = tau_min / num_per_steps;
    const double dt            = ( tau_min < Params.dt_scale*user_dt ) ? tau_min : std::min(user_dt, min_dt);
    const size_t num_pre_steps = static_cast<size_t>(std::ceil(tau_max / dt));
    const size_t num_steps     = static_cast<size_t>(std::ceil((t1 - t0) / dt)) + num_pre_steps + order;
    auto         solution      = dMatrix(num_steps + 1, dVec(N));
    auto         time_points   = dVec(num_steps + 1, 0.0);

    // Initialize pre-history using history function
    double pre_time = t0 - dt*num_pre_steps;
    size_t pre_idx  = 0;
    
    for (pre_idx = 0; pre_idx < num_pre_steps && pre_time < t0; ++pre_idx)
    {
        solution[pre_idx] = history_f(pre_time);
        time_points[pre_idx] = pre_time;
        pre_time += dt;
    }
    
    pre_idx              = num_pre_steps;
    solution[pre_idx]    = y0;
    pre_time             = t0;
    time_points[pre_idx] = pre_time;

    dMatrix f_hist(order, dVec(N));
    // f_hist stores derivatives in chronological order: [oldest, ..., newest]
    size_t hist_idx = 0;

    // Bootstrap with RK4 for the first `order` steps using DDE derivative
    dVec y_delayed(N);
    dMatrix y_delayed_all(num_delay, dVec(N));
    // Ensure we don't overrun allocated arrays: compute bootstrap upper index
    size_t bootstrap_upper = pre_idx;
    if (pre_idx < num_steps)
        bootstrap_upper = std::min(num_steps, pre_idx + static_cast<size_t>(order));
    for (size_t i = pre_idx; i < bootstrap_upper; ++i)
    {
        // Use RK4 for single step with DDE derivative
        const double t_curr = time_points[i];
        dVec y_curr = solution[i], y_next(N);
        dVec k1(N), k2(N), k3(N), k4(N);
        dVec y_k2(N), y_k3(N), y_k4(N);

        if (num_delay == 1)
        {
            if (delay_f)
            {
                y_delayed = interpolate_delay(solution, time_points, t_curr - delays[0], i, history_f, t0);
                k1 = delay_f(t_curr, y_curr, y_delayed);
                f_hist[hist_idx] = k1;

                // Compute k2
                y_k2 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k2[j] += 0.5 * dt * k1[j];

                y_delayed = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[0], i, history_f, t0);
                k2 = delay_f(t_curr + 0.5 * dt, y_k2, y_delayed);

                // Compute k3
                y_k3 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k3[j] += 0.5 * dt * k2[j];

                y_delayed = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[0], i, history_f, t0);
                k3 = delay_f(t_curr + 0.5 * dt, y_k3, y_delayed);

                // Compute k4
                y_k4 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k4[j] += dt * k3[j];

                y_delayed = interpolate_delay(solution, time_points, t_curr + dt - delays[0], i, history_f, t0);
                k4 = delay_f(t_curr + dt, y_k4, y_delayed);
            }
        }
        else
        {
            if (delay_f_m)
            {
                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr - delays[d], i, history_f, t0);
                }
                k1 = delay_f_m(t_curr, y_curr, y_delayed_all);
                f_hist[hist_idx] = k1;

                // Compute k2
                y_k2 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k2[j] += 0.5 * dt * k1[j];

                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[d], i, history_f, t0);
                }
                k2 = delay_f_m(t_curr + 0.5 * dt, y_k2, y_delayed_all);

                // Compute k3
                y_k3 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k3[j] += 0.5 * dt * k2[j];

                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[d], i, history_f, t0);
                }
                k3 = delay_f_m(t_curr + 0.5 * dt, y_k3, y_delayed_all);

                // Compute k4
                y_k4 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k4[j] += dt * k3[j];

                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr + dt - delays[d], i, history_f, t0);
                }
                k4 = delay_f_m(t_curr + dt, y_k4, y_delayed_all);
            }
        }
        
        // Update solution
        y_next = y_curr;
        for (size_t j = 0; j < N; ++j)
            y_next[j] += (dt / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        
        solution[i + 1] = y_next;
        time_points[i + 1] = t_curr + dt;

        hist_idx += 1;
    }

    // Bootstrap always fills f_hist with exactly `order` derivatives by construction
    assert(hist_idx == static_cast<size_t>(order));

    // Main Adams-Bashforth loop
    for (size_t i = pre_idx; i < num_steps; ++i)
    {
        const double t_curr = time_points[i];
        const dVec& y_curr = solution[i];
        dVec y_next = y_curr;

        // Apply Adams-Bashforth formula
        for (size_t j = 0; j < N; ++j)
        {
            double increment = 0.0;
            for (int k = 0; k < order; ++k)
                increment += abm_coefs::ab_coefs[order - 1][k] * f_hist[k][j];
            y_next[j] += dt * increment;
        }
        
        solution[i + 1] = y_next;
        time_points[i + 1] = t_curr + dt;

        // Update history by shifting and adding the new derivative
        for (int k = 0; k < order - 1; ++k)
            f_hist[k] = f_hist[k + 1];
        
        // Compute new derivative with delayed arguments
        const double t_next = time_points[i + 1];
        if (num_delay == 1 && delay_f)
        {
            y_delayed = interpolate_delay(solution, time_points, t_next - delays[0], i + 1, history_f, t0);
            f_hist[order - 1] = delay_f(t_next, y_next, y_delayed);
        }
        else if (num_delay > 1 && delay_f_m)
        {
            for (size_t d = 0; d < num_delay; ++d)
            {
                y_delayed_all[d] = interpolate_delay(solution, time_points, t_next - delays[d], i + 1, history_f, t0);
            }
            f_hist[order - 1] = delay_f_m(t_next, y_next, y_delayed_all);
        }
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for Adams-Bashforth solver
inline dMatrix adams_bashforth_solver_dde(
    MyFunc       deriv,
    const dVec&  y0,
    double       t0,
    double       t1,
    double       dt,
    int          order
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    params.order              = order;
    return adams_bashforth_solver_dde(params).solution;
}


// General explicit Adams-Bashforth-Moulton predictor-corrector solver for DDEs (order 1-10)
inline SolverResults adams_bashforth_moulton_solver_dde(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  y0        = Params.initial_conditions;
    const double t0        = Params.t0;
    const double t1        = Params.t1;
    const double user_dt   = Params.dt;
    const int    order     = Params.order;
    const int    corrector_iters = Params.iterations;
    const auto   delays    = Params.delay_times;
    const size_t N         = y0.size();
    const int    max_order = 10;
    const size_t num_delay = delays.size();
    const auto   history_f = Params.history_function;
    const double tau_max = (num_delay > 0) ? *std::max_element(delays.begin(), delays.end()) : 0.0;
    const double tau_min = (num_delay > 0) ? *std::min_element(delays.begin(), delays.end()) : 0.0;

    if (tau_max <= 1.0e-6 * user_dt)
    {
        std::cout << "Delays are negligible; Max{delay} = " << tau_max << std::endl << "Using ODE solver." << std::endl;
        return adams_bashforth_moulton_solver(Params);
    }
    if (num_delay == 0)
    {
        std::cout << "Delays not provided! Using ODE solver.\n";
        return adams_bashforth_moulton_solver(Params);
    }

    if (order < 1 || order > max_order)
        throw std::invalid_argument("ABM order must be between 1 and 10");

    // Setting the RHS
    DelayFunc      delay_f   = nullptr;
    DelayFuncMulti delay_f_m = nullptr;
    if (num_delay == 1)
        delay_f = Params.delay_derivative;
    else
        delay_f_m = Params.delay_derivative_multi;

    if (delay_f == nullptr && delay_f_m == nullptr)
        throw std::invalid_argument("DDE derivative function must be provided for non-zero delays");
    
    // Initialize solution storage
    const size_t num_per_steps = std::max(static_cast<size_t>(8), static_cast<size_t>(std::ceil(1.0 / user_dt)));
    const double min_dt        = tau_min / num_per_steps;
    const double dt            = (tau_min < Params.dt_scale * user_dt) ? tau_min : std::min(user_dt, min_dt);
    const size_t num_pre_steps = static_cast<size_t>(std::ceil(tau_max / dt));
    const size_t num_steps     = static_cast<size_t>(std::ceil((t1 - t0) / dt)) + num_pre_steps + order;
    auto         solution      = dMatrix(num_steps + 1, dVec(N));
    auto         time_points   = dVec(num_steps + 1, 0.0);

    // Initialize pre-history using history function
    double pre_time = t0 - dt * num_pre_steps;
    size_t pre_idx  = 0;
    
    for (pre_idx = 0; pre_idx < num_pre_steps && pre_time < t0; ++pre_idx)
    {
        solution[pre_idx] = history_f(pre_time);
        time_points[pre_idx] = pre_time;
        pre_time += dt;
    }
    
    pre_idx              = num_pre_steps;
    solution[pre_idx]    = y0;
    pre_time             = t0;
    time_points[pre_idx] = pre_time;

    dMatrix f_hist(order, dVec(N));
    // f_hist stores derivatives in chronological order: [oldest, ..., newest]
    size_t hist_idx = 0;

    // Bootstrap with RK4 for the first `order` steps using DDE derivative
    dVec y_delayed(N);
    dMatrix y_delayed_all(num_delay, dVec(N));
    // Ensure we don't overrun allocated arrays: compute bootstrap upper index
    size_t bootstrap_upper = pre_idx;
    if (pre_idx < num_steps)
        bootstrap_upper = std::min(num_steps, pre_idx + static_cast<size_t>(order));
    for (size_t i = pre_idx; i < bootstrap_upper; ++i)
    {
        // Use RK4 for single step with DDE derivative
        const double t_curr = time_points[i];
        dVec y_curr = solution[i], y_next(N);
        dVec k1(N), k2(N), k3(N), k4(N);
        dVec y_k2(N), y_k3(N), y_k4(N);

        if (num_delay == 1)
        {
            if (delay_f)
            {
                y_delayed = interpolate_delay(solution, time_points, t_curr - delays[0], i, history_f, t0);
                k1 = delay_f(t_curr, y_curr, y_delayed);
                f_hist[hist_idx] = k1;

                // Compute k2
                y_k2 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k2[j] += 0.5 * dt * k1[j];

                y_delayed = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[0], i, history_f, t0);
                k2 = delay_f(t_curr + 0.5 * dt, y_k2, y_delayed);

                // Compute k3
                y_k3 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k3[j] += 0.5 * dt * k2[j];

                y_delayed = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[0], i, history_f, t0);
                k3 = delay_f(t_curr + 0.5 * dt, y_k3, y_delayed);

                // Compute k4
                y_k4 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k4[j] += dt * k3[j];

                y_delayed = interpolate_delay(solution, time_points, t_curr + dt - delays[0], i, history_f, t0);
                k4 = delay_f(t_curr + dt, y_k4, y_delayed);
            }
        }
        else
        {
            if (delay_f_m)
            {
                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr - delays[d], i, history_f, t0);
                }
                k1 = delay_f_m(t_curr, y_curr, y_delayed_all);
                f_hist[hist_idx] = k1;

                // Compute k2
                y_k2 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k2[j] += 0.5 * dt * k1[j];

                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[d], i, history_f, t0);
                }
                k2 = delay_f_m(t_curr + 0.5 * dt, y_k2, y_delayed_all);

                // Compute k3
                y_k3 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k3[j] += 0.5 * dt * k2[j];

                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr + 0.5 * dt - delays[d], i, history_f, t0);
                }
                k3 = delay_f_m(t_curr + 0.5 * dt, y_k3, y_delayed_all);

                // Compute k4
                y_k4 = y_curr;
                for (size_t j = 0; j < N; ++j) y_k4[j] += dt * k3[j];

                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_curr + dt - delays[d], i, history_f, t0);
                }
                k4 = delay_f_m(t_curr + dt, y_k4, y_delayed_all);
            }
        }
        
        // Update solution
        y_next = y_curr;
        for (size_t j = 0; j < N; ++j)
            y_next[j] += (dt / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        
        solution[i + 1] = y_next;
        time_points[i + 1] = t_curr + dt;

        hist_idx += 1;
    }

    // Bootstrap always fills f_hist with exactly `order` derivatives by construction
    assert(hist_idx == static_cast<size_t>(order));

    // Main Adams-Bashforth-Moulton predictor-corrector loop
    auto f_corr = dMatrix(order + 1, dVec(N));
    for (size_t i = pre_idx; i < num_steps; ++i)
    {
        const double t_curr = time_points[i];
        const double t_next = t_curr + dt;
        const dVec   y_curr = solution[i];

        // Predictor (Adams-Bashforth)
        // std::cout << "time = " << t_curr << '\n';
        for (size_t j = 0; j < N; ++j)
        {
            double increment = 0.0;
            // std::cout << "Y[" << j << "] = " << y_curr[j] << '\n';
            for (int k = 0; k < order; ++k)
                increment += abm_coefs::ab_coefs[order - 1][k] * f_hist[k][j];
            // std::cout << "correction iteration = Main | i = " << i << " | increment = " << increment*dt << '\n';
            solution[i + 1][j] = y_curr[j] + dt * increment;
        }

        time_points[i + 1] = t_next;

        // Corrector (Adams-Moulton)
        for (int k = 0; k < order; ++k)
            f_corr[k] = f_hist[k];
        
        for (int iter = 0; iter < corrector_iters; ++iter)
        {
            // Evaluate DDE derivative at (t_next, solution[i + 1])
            if (num_delay == 1 && delay_f)
            {
                y_delayed = interpolate_delay(solution, time_points, t_next - delays[0], i + 1, history_f, t0);
                f_corr[order] = delay_f(t_next, solution[i + 1], y_delayed);
            }
            else if (num_delay > 1 && delay_f_m)
            {
                for (size_t d = 0; d < num_delay; ++d)
                {
                    y_delayed_all[d] = interpolate_delay(solution, time_points, t_next - delays[d], i + 1, history_f, t0);
                }
                f_corr[order] = delay_f_m(t_next, solution[i + 1], y_delayed_all);
            }
            
            // Apply Adams-Moulton formula
            for (size_t j = 0; j < N; ++j)
            {
                double increment = 0.0;
                for (int k = 0; k < order + 1; ++k)
                    increment += abm_coefs::am_coefs[order - 1][k] * f_corr[k][j];
                // std::cout << "correction iteration = " << iter+1 << " | i = " << i << " | increment = " << increment*dt << '\n';
                // std::cout << "y_delay[" << j << "] = " << y_delayed[j] << " | f_corr[" << j << "] = " << f_corr[order][j] << '\n';
                solution[i + 1][j] = y_curr[j] + dt * increment;
            }
        }        

        // Update history by shifting and adding the new derivative
        for (int k = 0; k < order - 1; ++k)
            f_hist[k] = f_hist[k + 1];
        
        // Compute final derivative for next step's history
        if (num_delay == 1 && delay_f)
        {
            y_delayed = interpolate_delay(solution, time_points, t_next - delays[0], i + 1, history_f, t0);
            f_hist[order - 1] = delay_f(t_next, solution[i + 1], y_delayed);
        }
        else if (num_delay > 1 && delay_f_m)
        {
            for (size_t d = 0; d < num_delay; ++d)
            {
                y_delayed_all[d] = interpolate_delay(solution, time_points, t_next - delays[d], i + 1, history_f, t0);
            }
            f_hist[order - 1] = delay_f_m(t_next, solution[i + 1], y_delayed_all);
        }
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for Adams-Bashforth-Moulton DDE solver
inline dMatrix adams_bashforth_moulton_solver_dde(
    MyFunc         deriv,
    DelayFunc      dderiv,
    DelayFuncMulti dderivm,
    const dVec&    y0,
    double         t0,
    double         t1,
    double         dt,
    int            order,
    Func           history_f,
    int            corrector_iters
)
{
    auto params                   = SolverParameters{};
    params.derivative             = deriv;
    params.delay_derivative       = dderiv;
    params.delay_derivative_multi = dderivm;
    params.initial_conditions     = y0;
    params.t0                     = t0;
    params.t1                     = t1;
    params.dt                     = dt;
    params.order                  = order;
    params.history_function       = history_f;
    params.iterations             = corrector_iters;
    return adams_bashforth_moulton_solver_dde(params).solution;
}
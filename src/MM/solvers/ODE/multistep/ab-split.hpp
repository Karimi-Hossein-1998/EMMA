#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../namespaces/abm-coefs.hpp"

// Split API for Adams-Bashforth bootstrapping and predictor-only integration
// - ab_bootstrap_rk4: generates a history of `order` states for times
//   t in [t0 - (order-1)*dt, ..., t0] (increasing order). It steps "backward"
//   from t0 using a backward RK4 step (implemented by a negative step size),
//   but returns the history in time-increasing order for convenience.
//   The returned SolverResults contains:
//     - solution: matrix of `order` states (earliest -> t0)
//     - time_points: corresponding times (earliest -> t0)
//     - step_sizes: populated with dt repeated (size = order - 1)
//     - steps_history: populated with the derivative evaluations at each history time
// - ab_predictor_only: uses an existing history (as produced by `ab_bootstrap_rk4`)
//   to integrate forward from t0 to t1 using the explicit Adams-Bashforth formula
//   (fixed dt). The history must contain at least `order` entries ending at t0.
// - adams_bashforth_split: convenience wrapper that calls bootstrap then
//   predictor, returning the forward-in-time solution (from t0 to t1).

// Note: These functions are intentionally minimal and do not call the global
// validator. Callers may invoke `validate_solver_parameters` if they wish.

// This is not a main file and is not supposed to be maintained (unsless proven to be necessary and useful in the future. So, please stick to the ab_solver.hpp file for any changes regarding the Adams-Bashforth method or any improvements/adjustments).

inline SolverResults ab_bootstrap_rk4(const SolverParameters& Params)
{
    const auto&  f     = Params.derivative;
    const auto&  y0    = Params.initial_conditions; // y(t0)
    const double t0    = Params.t0;
    const double dt    = Params.dt;
    const int    order = std::max(1, Params.order);
    const size_t N     = y0.size();

    // Prepare containers: we'll build arrays backwards then reverse to increasing time
    std::vector<dVec> sol_rev(order);
    std::vector<double> times_rev(order);

    // last entry is the known y0 at t0
    sol_rev[order - 1] = y0;
    times_rev[order - 1] = t0;

    // local RK4 single-step (advances by step h; h can be negative for backward stepping)
    auto rk4_step = [&](double t, const dVec& y, double h) {
        dVec k1 = f(t, y);
        dVec y2(N);
        for (size_t j = 0; j < N; ++j) y2[j] = y[j] + 0.5 * h * k1[j];
        dVec k2 = f(t + 0.5 * h, y2);
        for (size_t j = 0; j < N; ++j) y2[j] = y[j] + 0.5 * h * k2[j];
        dVec k3 = f(t + 0.5 * h, y2);
        dVec y4(N);
        for (size_t j = 0; j < N; ++j) y4[j] = y[j] + h * k3[j];
        dVec k4 = f(t + h, y4);

        dVec y_new(N);
        for (size_t j = 0; j < N; ++j)
            y_new[j] = y[j] + (h / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        return y_new;
    };

    // Step backward (h = -dt) to generate earlier points
    const double h = -dt;
    for (int i = order - 2; i >= 0; --i)
    {
        const double t_next = times_rev[i + 1];
        const dVec& y_next = sol_rev[i + 1];
        dVec y_prev = rk4_step(t_next, y_next, h); // step backward
        sol_rev[i] = std::move(y_prev);
        times_rev[i] = t_next + h; // equal to t_next - dt
    }

    // Reverse to increasing order
    auto solution = dMatrix(order);
    auto time_points = dVec(order);
    for (int i = 0; i < order; ++i)
    {
        solution[i] = std::move(sol_rev[i]);
        time_points[i] = times_rev[i];
    }

    // Compute derivative history and step sizes
    dMatrix derivatives(order, dVec(N, 0.0));
    for (int i = 0; i < order; ++i)
        derivatives[i] = f(time_points[i], solution[i]);

    dVec step_sizes;
    if (order >= 2) step_sizes.assign(order - 1, dt);

    SolverResults results;
    results.solution = std::move(solution);
    results.time_points = std::move(time_points);
    results.step_sizes = std::move(step_sizes);
    results.steps_history = std::move(derivatives); // store derivative history here for convenience

    return results;
}

inline SolverResults ab_predictor_only(const SolverParameters& Params, const SolverResults& history)
{
    const auto& f   = Params.derivative;
    const double t0 = Params.t0;
    const double t1 = Params.t1;
    const double dt = Params.dt;
    const int order = std::max(1, Params.order);

    if (history.solution.size() < static_cast<size_t>(order))
        throw std::invalid_argument("History must contain at least 'order' states ending at t0");

    // We'll use the last `order` entries from history (assumed in increasing time)
    const size_t hist_sz = history.solution.size();
    const size_t start_idx = hist_sz - static_cast<size_t>(order);

    // Initialize f_hist: oldest -> newest
    dMatrix f_hist(order, dVec(history.solution[0].size(), 0.0));
    if (!history.steps_history.empty() && history.steps_history.size() >= hist_sz)
    {
        // Use provided derivatives if available
        for (int k = 0; k < order; ++k)
            f_hist[k] = history.steps_history[start_idx + k];
    }
    else
    {
        // Compute derivatives
        for (int k = 0; k < order; ++k)
            f_hist[k] = f(history.time_points[start_idx + k], history.solution[start_idx + k]);
    }

    // Determine number of forward steps (consistent with fixed-step non-adaptive solvers)
    const size_t num_steps = static_cast<size_t>((t1 - t0) / dt);

    // Result containers: include t0 as first point
    auto solution = dMatrix(num_steps + 1, dVec(history.solution.back().size()));
    auto time_points = dVec(num_steps + 1);

    // initial state at t0 is the latest history element
    solution[0] = history.solution.back();
    time_points[0] = history.time_points.back();

    dVec y_curr = solution[0];
    double t_curr = time_points[0];

    // AB predictor loop
    for (size_t step = 0; step < num_steps; ++step)
    {
        // compute increment using coefficients: f_hist is oldest->newest
        dVec increment(y_curr.size(), 0.0);
        for (int k = 0; k < order; ++k)
        {
            const auto& fk = f_hist[k];
            const double coef = abm_coefs::ab_coefs[order - 1][k];
            for (size_t j = 0; j < y_curr.size(); ++j)
                increment[j] += coef * fk[j];
        }

        dVec y_next = y_curr;
        for (size_t j = 0; j < y_next.size(); ++j)
            y_next[j] += dt * increment[j];

        double t_next = t_curr + dt;

        solution[step + 1] = y_next;
        time_points[step + 1] = t_next;

        // update history: shift left and append newest derivative
        for (int k = 0; k < order - 1; ++k)
            f_hist[k] = f_hist[k + 1];
        f_hist[order - 1] = f(t_next, y_next);

        // prepare for next iteration
        y_curr = std::move(y_next);
        t_curr = t_next;
    }

    SolverResults results;
    results.solution = std::move(solution);
    results.time_points = std::move(time_points);
    // populate step_sizes as constant dt steps
    results.step_sizes = dVec(num_steps, dt);
    return results;
}

inline SolverResults adams_bashforth_split(const SolverParameters& Params)
{
    // Bootstrap then predictor
    auto history = ab_bootstrap_rk4(Params);
    return ab_predictor_only(Params, history);
}

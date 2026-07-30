#pragma once
#include "../../typedefs/header.hpp"

/**
 * @file validator.hpp
 * @brief Lightweight solver parameter validator.
 *
 * Provides inexpensive sanity checks for user-supplied `SolverParameters`.
 * By default the checks are "cheap" (no derivative evaluation) and suitable
 * for use in production code paths where performance is important.
 *
 * If `deep_check == true` the validator evaluates the derivative once at
 * `t0` to ensure that the user-provided `derivative` has the correct
 * behaviour and returns a vector of the expected size.
 *
 * Behavior:
 *  - Returns `true` when parameters look acceptable.
 *  - Returns `false` on failure and, if `out_err` is provided, writes a
 *    clear, actionable error message naming the offending parameter and
 *    suggesting a fix.
 *
 * Complexity: O(1) (unless `deep_check == true`, then 1 derivative evaluation).
 *
 * Examples:
 *  std::string err;
 *  if (!validate_solver_parameters(params, false, &err)) throw std::invalid_argument(err);
 *
 */

// Validate solver parameters and provide detailed, actionable error messages.
inline bool validate_solver_parameters(
    const SolverParameters& params,
    bool deep_check = false,
    std::string* out_err = nullptr,
    size_t max_steps = 100000000UL
)
{
    // Check dt
    if (!std::isfinite(params.dt) || params.dt <= 0.0)
    {
        if (out_err)
        {
            std::ostringstream ss;
            ss << "Invalid parameter 'dt' = " << params.dt
               << ": step size must be finite and > 0. Suggestion: use dt = (t1 - t0)/N for some reasonable N.";
            *out_err = ss.str();
        }
        return false;
    }

    // Check t0/t1
    if (!std::isfinite(params.t0) || !std::isfinite(params.t1))
    {
        if (out_err)
        {
            std::ostringstream ss;
            ss << "Invalid parameter: 't0' or 't1' is not finite (t0=" << params.t0 << ", t1=" << params.t1 << ")";
            *out_err = ss.str();
        }
        return false;
    }

    // Check initial conditions
    if (params.initial_conditions.empty())
    {
        if (out_err) *out_err = "Invalid parameter 'initial_conditions': vector is empty. Provide initial state values (e.g., {1.0, 0.0}).";
        return false;
    }

    // Forward-integration assumption
    if (params.t1 < params.t0)
    {
        if (out_err)
        {
            std::ostringstream ss;
            ss << "Invalid time interval: t1 < t0 (t0=" << params.t0 << ", t1=" << params.t1 << ")";
            ss << ". The validator assumes forward integration (t1 >= t0).";
            ss << " Suggestion: swap t0/t1 for forward integration or adjust your caller for backward integration.";
            *out_err = ss.str();
        }
        return false;
    }

    // Compute number of steps (ceiling so remainder short final step is accounted for)
    const double interval = params.t1 - params.t0;
    size_t num_steps = 0;
    if (interval > 0.0)
    {
        num_steps = static_cast<size_t>(std::ceil(interval / params.dt));
        if (num_steps == 0) num_steps = 1; // very small interval
    }

    if (num_steps > max_steps)
    {
        if (out_err)
        {
            std::ostringstream ss;
            ss << "Invalid parameters: expected number of steps (" << num_steps << ") exceeds max_steps (" << max_steps << ").";
            ss << " Suggestion: increase 'dt' or reduce the interval (t1 - t0).";
            *out_err = ss.str();
        }
        return false;
    }

    // Optional deep check: evaluate derivative once at t0 to verify size and behaviour
    if (deep_check)
    {
        try
        {
            auto v = params.derivative(params.t0, params.initial_conditions);
            if (v.size() != params.initial_conditions.size())
            {
                if (out_err)
                {
                    std::ostringstream ss;
                    ss << "Derivative returned vector of size " << v.size()
                       << " at t0 = " << params.t0 << ", but expected size " << params.initial_conditions.size() << ".";
                    ss << " Suggestion: ensure your 'derivative' returns dVec of the same length as initial_conditions.";
                    *out_err = ss.str();
                }
                return false;
            }
        }
        catch (const std::exception& e)
        {
            if (out_err) {
                std::ostringstream ss;
                ss << "Derivative threw an exception during evaluation at t0 = " << params.t0 << ": " << e.what();
                ss << " Suggestion: run the derivative separately to reproduce the error.";
                *out_err = ss.str();
            }
            return false;
        }
        catch (...)
        {
            if (out_err) *out_err = "Derivative threw an unknown (non-std) exception during evaluation at t0";
            return false;
        }
    }

    return true;
}

// Convenience wrapper that throws std::invalid_argument on failure
inline void validate_solver_parameters_or_throw(
    const SolverParameters& params,
    bool deep_check = false,
    size_t max_steps = 100000000UL
)
{
    std::string errmsg;
    if (!validate_solver_parameters(params, deep_check, &errmsg, max_steps))
    {
        throw std::invalid_argument(errmsg);
    }
}

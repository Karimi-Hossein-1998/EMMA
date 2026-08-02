#pragma once
#include "../typedefs/header.hpp"

// Helper function for error estimation
double estimate_error(const dVec& y, const dVec& y_ref, const SolverParameters& Params)
{
    double       error   = 0.0;
    const double epsilon = std::max(1e-15, Params.absolute_tol);
    double       denom   = epsilon;
    if (Params.weightedError)
    {
        // Weighted component-wise error
        for (size_t i = 0; i < y.size(); ++i)
        {
            double diff = std::abs(y[i] - y_ref[i]);
            denom = Params.absolute_tol + std::max(std::abs(y_ref[i]), std::abs(y[i]));
            if (denom < epsilon) denom = epsilon; // avoid division by zero
            error = std::max(error, diff / denom);
        }
    }
    else if (Params.normError)
    {
        // L2 norm-based error
        double norm_diff = 0.0;
        double norm_sol  = 0.0;
        for (size_t i = 0; i < y.size(); ++i)
        {
            double diff = y[i] - y_ref[i];
            norm_diff += diff * diff;
            norm_sol += y_ref[i] * y_ref[i];
        }
        denom = Params.absolute_tol + std::sqrt(norm_sol);
        if (denom < epsilon) denom = epsilon; // avoid division by zero
        error = std::sqrt(norm_diff) / denom;
    }
    else
    {
        for (size_t i = 0; i < y.size(); ++i)
        {
            double diff = std::abs(y[i] - y_ref[i]);
            denom = Params.absolute_tol + std::abs(y_ref[i]);
            if (denom < epsilon) denom = epsilon; // avoid division by zero
            error = std::max(error, diff / denom);
        }
    }
    return error;
}

// Helper function for step size adjustment
double adjust_step_size(double current_dt, double error, const SolverParameters& Params, size_t p)
{
    if (error > Params.localTol && current_dt > Params.minDt)
    {
        return std::max(std::min(current_dt * Params.decreaseFactor * std::pow(Params.localTol / error, 1.0 / (p + 1.0)),
                                Params.maxDt),
                       Params.minDt);
    }
    else if (error < Params.localTol * Params.localTolErrorRatio || current_dt <= Params.minDt)
    {
        return std::max(std::min(current_dt * Params.increaseFactor * std::pow(Params.localTol / error, 1.0 / (p + 1.0)),
                                Params.maxDt),
                       Params.minDt);
    }
    return current_dt;
}

#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../utility/adaptive-solvers-utility.hpp"

// RK12 solver with restructured error estimation
SolverResults rk12_custom_solver(const SolverParameters& Params)
{
    auto deriv = Params.derivative;
    auto y     = Params.initial_conditions;
    auto t     = Params.t0;
    auto t1    = Params.t1;
    auto dt    = Params.dt;

    dMatrix result;
    dVec   Times;
    dVec   Errors;
    dVec   Steps;
    result.push_back(y);
    Times.push_back(t);

    dMatrix      Errors_History;
    dMatrix      Steps_History;
    wVec Failed_Trials;
    dVec        Error_Attempt(Params.max_trial, 0.0);
    dVec        Step_Attempt(Params.max_trial, 0.0);

    // Case 1: No error estimation - simple RK1
    if (!Params.error_estimate)
    {
        while (t < t1)
        {
            auto k1 = deriv(t, y);
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt * k1[i];
            
            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 2: Error estimation without variable steps
    else if (!Params.variable_steps)
    {
        auto dt_half = dt * 0.5;
        auto y_old   = y;
        auto y_rk2   = y_old;
        while (t < t1)
        {
            
            // Compute RK2 solution first (more efficient)
            auto k1 = deriv(t, y_old);

            // Second stage (midpoint)
            y_rk2 = y;
            y_old = y;
            for (size_t i = 0; i < y.size(); ++i)
                y_old[i] += dt_half * k1[i];
            auto k2 = deriv(t + dt_half, y_old);

            // Update using RK2
            for (size_t i = 0; i < y_rk2.size(); ++i)
                y_rk2[i] += dt * k2[i];

            // Compute RK1 solution (using k1 from RK2)
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt * k1[i];

            // Calculate error
            auto error = estimate_error(y, y_rk2, Params);
            Errors.push_back(error);
            Steps.push_back(dt);

            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 3: Full adaptive step size with error estimation
    else
    {
        while (t < t1)
        {
            bool   step_accepted = false;
            auto   y_old         = y;
            double current_dt    = dt;
            size_t trials        = 0;

            while (!step_accepted)
            {
                y = y_old;  // Reset to previous state

                // Compute RK2 solution first (more efficient)
                auto y_rk2   = y_old;
                auto dt_half = current_dt * 0.5;
                auto k1      = deriv(t, y_old);

                // Second stage (midpoint)
                auto y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += dt_half * k1[i];
                auto k2 = deriv(t + dt_half, y_mid);

                // Update using RK2
                for (size_t i = 0; i < y_rk2.size(); ++i)
                    y_rk2[i] += current_dt * k2[i];

                // Compute RK1 solution (using k1 from RK2)
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] += current_dt * k1[i];

                // Calculate error and adjust step size
                auto error = estimate_error(y, y_rk2, Params);
                Error_Attempt[trials] = error;
                Step_Attempt[trials] = current_dt;
                trials++;

                auto new_dt = adjust_step_size(current_dt, error, Params, 1);
                if (new_dt != current_dt)
                {
                    current_dt = new_dt;
                    continue;
                }

                step_accepted = true;
                dt = new_dt;
            }

            if (Params.attempts_history)
            {
                Errors_History.push_back(Error_Attempt);
                Steps_History.push_back(Step_Attempt);
                Failed_Trials.push_back(trials);
            }

            Errors.push_back(Error_Attempt[trials - 1]);
            Steps.push_back(current_dt);
            t += current_dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }

    SolverResults Results;
    Results.solution = result;
    Results.time_points = Times;
    if (Params.error_estimate)
    {
        Results.errors = Errors;
        if (Params.variable_steps)
        {
            Results.step_sizes = Steps;
            if (Params.attempts_history)
            {
                Results.errors_history = Errors_History;
                Results.steps_history = Steps_History;
                Results.failed_trials = Failed_Trials;
            }
        }
    }
    return Results;
}

// RK23 solver with restructured error estimation
SolverResults rk23_custom_solver(const SolverParameters& Params)
{
    auto deriv = Params.derivative;
    auto y     = Params.initial_conditions;
    auto t     = Params.t0;
    auto t1    = Params.t1;
    auto dt    = Params.dt;

    dMatrix result;
    dVec   Times;
    dVec   Errors;
    dVec   Steps;
    result.push_back(y);
    Times.push_back(t);

    dMatrix      Errors_History;
    dMatrix      Steps_History;
    wVec Failed_Trials;
    dVec        Error_Attempt(Params.max_trial, 0.0);
    dVec        Step_Attempt(Params.max_trial, 0.0);

    auto y_old = y;
    auto y_rk3 = y;
    auto y_mid = y;
    auto k1    = y;
    auto k2    = y;
    auto k3    = y;
    // Case 1: No error estimation - simple RK2
    if (!Params.error_estimate)
    {
        auto dt_half  = dt * 0.5;
        while (t < t1)
        {
            k1 = deriv(t, y);
            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_old);
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt * k2[i];
            
            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 2: Error estimation without variable steps
    else if (!Params.variable_steps)
    {
        auto dt_half  = dt * 0.5;
        auto dt_2     = dt * 2.0; 
        auto dt_23    = dt_2 / 3.0;
        auto dt_sixth = dt / 6.0;
        while (t < t1)
        {
            
            // Compute RK3 solution first (more efficient)
            y_old = y;
            y_rk3 = y;
            auto k1 = deriv(t, y_old);

            // Second stage (midpoint)
            y_old = y;
            for (size_t i = 0; i < y.size(); ++i)
                y_old[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_old);

            // Third stage
            y_old = y;
            for (size_t i = 0; i < y.size(); ++i)
                y_old[i] += - dt * k1[i] + dt_2 * k2[i];
            k3 = deriv(t + dt, y_old);

            // Update using RK3
            for (size_t i = 0; i < y_rk3.size(); ++i)
                y_rk3[i] += dt_sixth * (k1[i] + k3[i]) + dt_23 * k2[i];

            // Compute RK2 solution (using k1 and k2 from RK3)
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt * k2[i];

            // Calculate error
            auto error = estimate_error(y, y_rk3, Params);
            Errors.push_back(error);
            Steps.push_back(dt);

            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 3: Full adaptive step size with error estimation
    else
    {
        while (t < t1)
        {
            bool   step_accepted = false;
            y_old                = y;
            double current_dt    = dt;
            size_t trials        = 0;

            while (!step_accepted)
            {
                y = y_old;  // Reset to previous state

                // Compute RK3 solution first (more efficient)
                y_rk3        = y_old;
                auto dt_half = current_dt * 0.5;
                auto dt_2    = current_dt * 2.0;
                auto dt_23    = dt_2 / 3.0;
                auto dt_sixth = dt / 6.0;
                
                k1 = deriv(t, y_old);

                y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += dt_half * k1[i];
                k2 = deriv(t + dt_half, y_mid);

                y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += - current_dt * k1[i] + dt_2 * k2[i];
                k3 = deriv(t + current_dt, y_mid);

                for (size_t i = 0; i < y_rk3.size(); ++i)
                    y_rk3[i] += dt_sixth * (k1[i] + k3[i]) + dt_23 * k2[i];

                for (size_t i = 0; i < y.size(); ++i)
                    y[i] += current_dt * k2[i];

                // Calculate error and adjust step size
                auto error = estimate_error(y, y_rk3, Params);
                Error_Attempt[trials] = error;
                Step_Attempt[trials] = current_dt;
                trials++;

                auto new_dt = adjust_step_size(current_dt, error, Params, 2);
                if (new_dt != current_dt)
                {
                    current_dt = new_dt;
                    continue;
                }

                step_accepted = true;
                dt = new_dt;
            }

            if (Params.attempts_history)
            {
                Errors_History.push_back(Error_Attempt);
                Steps_History.push_back(Step_Attempt);
                Failed_Trials.push_back(trials);
            }

            Errors.push_back(Error_Attempt[trials - 1]);
            Steps.push_back(current_dt);
            t += current_dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }

    SolverResults Results;
    Results.solution = result;
    Results.time_points = Times;
    if (Params.error_estimate)
    {
        Results.errors = Errors;
        if (Params.variable_steps)
        {
            Results.step_sizes = Steps;
            if (Params.attempts_history)
            {
                Results.errors_history = Errors_History;
                Results.steps_history = Steps_History;
                Results.failed_trials = Failed_Trials;
            }
        }
    }
    return Results;
}

// RK34 solver with restructured error estimation
SolverResults rk34_custom_solver(const SolverParameters& Params)
{
    auto deriv = Params.derivative;
    auto y     = Params.initial_conditions;
    auto t     = Params.t0;
    auto t1    = Params.t1;
    auto dt    = Params.dt;

    dMatrix result;
    dVec   Times;
    dVec   Errors;
    dVec   Steps;
    result.push_back(y);
    Times.push_back(t);

    dMatrix      Errors_History;
    dMatrix      Steps_History;
    wVec Failed_Trials;
    dVec        Error_Attempt(Params.max_trial, 0.0);
    dVec        Step_Attempt(Params.max_trial, 0.0);

    auto y_temp = y;
    auto y_rk4  = y;
    auto y_old  = y;
    auto k1     = y;
    auto k2     = y;
    auto k3     = y;
    auto k4     = y;
    // Case 1: No error estimation - simple RK3
    if (!Params.error_estimate)
    {
        auto dt_half = dt * 0.5;
        auto dt_3_4  = dt * 0.75;
        auto dt_2_9  = dt * 2.0 / 9.0;
        auto dt_1_3  = dt / 3.0;
        auto dt_4_9  = dt * 4.0 / 9.0;
        while (t < t1)
        {
            k1 = deriv(t, y);
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_temp);

            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_3_4 * k2[i];
            k3 = deriv(t + dt_3_4, y_temp);

            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt_2_9 * k1[i] + dt_1_3 * k2[i] + dt_4_9 * k3[i];

            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 2: Error estimation without variable steps
    else if (!Params.variable_steps)
    {
        auto dt_half  = dt * 0.5;
        auto dt_3_4   = dt * 0.75;
        auto dt_2_9   = dt * 2.0 / 9.0;
        auto dt_1_3   = dt / 3.0;
        auto dt_4_9   = dt * 4.0 / 9.0;
        auto dt_7_24  = dt * 7.0 / 24.0;
        auto dt_1_4   = dt * 0.25;
        auto dt_1_8   = dt * 0.125;
        while (t < t1)
        {
            // Compute RK4 solution
            k1 = deriv(t, y);
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_temp);

            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_3_4 * k2[i];
            k3 = deriv(t + dt_3_4, y_temp);

            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt * (2.0/9.0 * k1[i] + 1.0/3.0 * k2[i] + 4.0/9.0 * k3[i]);
            k4 = deriv(t + dt, y_temp);

            // Update using RK4 (higher order)
            for (size_t i = 0; i < y_rk4.size(); ++i)
                y_rk4[i] = y[i] + dt_7_24 * k1[i] + dt_1_4 * k2[i] + dt_1_3 * k3[i] + dt_1_8 * k4[i];

            // Update using RK3 (lower order)
            for (size_t i = 0; i < y.size(); ++i)
                y[i] = y[i] + dt_2_9 * k1[i] + dt_1_3 * k2[i] + dt_4_9 * k3[i];

            // Calculate error
            auto error = estimate_error(y, y_rk4, Params);
            Errors.push_back(error);
            Steps.push_back(dt);

            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 3: Full adaptive step size with error estimation
    else
    {
        while (t < t1)
        {
            bool   step_accepted = false;
            y_old                = y;
            double current_dt    = dt;
            size_t trials        = 0;

            while (!step_accepted)
            {
                y      = y_old;  // Reset to previous state
                y_rk4  = y_old;
                auto dt_half  = current_dt * 0.5;
                auto dt_3_4   = current_dt * 0.75;
                auto dt_2_9   = current_dt * 2.0 / 9.0;
                auto dt_1_3   = current_dt / 3.0;
                auto dt_4_9   = current_dt * 4.0 / 9.0;
                auto dt_7_24  = current_dt * 7.0 / 24.0;
                auto dt_1_4   = current_dt * 0.25;
                auto dt_1_8   = current_dt * 0.125;

                // Compute RK4 solution
                k1 = deriv(t, y_old);
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_half * k1[i];
                k2 = deriv(t + dt_half, y_temp);

                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_3_4 * k2[i];
                k3 = deriv(t + dt_3_4, y_temp);

                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += current_dt * (2.0/9.0 * k1[i] + 1.0/3.0 * k2[i] + 4.0/9.0 * k3[i]);
                k4 = deriv(t + current_dt, y_temp);

                // Update using RK4 (higher order)
                for (size_t i = 0; i < y_rk4.size(); ++i)
                    y_rk4[i] = y_old[i] + dt_7_24 * k1[i] + dt_1_4 * k2[i] + dt_1_3 * k3[i] + dt_1_8 * k4[i];

                // Update using RK3 (lower order)
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] = y_old[i] + dt_2_9 * k1[i] + dt_1_3 * k2[i] + dt_4_9 * k3[i];

                // Calculate error and adjust step size
                auto error = estimate_error(y, y_rk4, Params);
                Error_Attempt[trials] = error;
                Step_Attempt[trials]  = current_dt;
                trials++;

                auto new_dt = adjust_step_size(current_dt, error, Params, 3);
                if (new_dt != current_dt)
                {
                    current_dt = new_dt;
                    continue;
                }

                step_accepted = true;
                dt = new_dt;
            }

            if (Params.attempts_history)
            {
                Errors_History.push_back(Error_Attempt);
                Steps_History.push_back(Step_Attempt);
                Failed_Trials.push_back(trials);
            }

            Errors.push_back(Error_Attempt[trials - 1]);
            Steps.push_back(current_dt);
            t += current_dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }

    SolverResults Results;
    Results.solution = result;
    Results.time_points = Times;
    if (Params.error_estimate)
    {
        Results.errors = Errors;
        if (Params.variable_steps)
        {
            Results.step_sizes = Steps;
            if (Params.attempts_history)
            {
                Results.errors_history = Errors_History;
                Results.steps_history = Steps_History;
                Results.failed_trials = Failed_Trials;
            }
        }
    }
    return Results;
}

// RK45 solver with restructured error estimation
SolverResults rk45_custom_solver(const SolverParameters& Params)
{
    auto deriv = Params.derivative;
    auto y     = Params.initial_conditions;
    auto t     = Params.t0;
    auto t1    = Params.t1;
    auto dt    = Params.dt;

    dMatrix result;
    dVec   Times;
    dVec   Errors;
    dVec   Steps;
    result.push_back(y);
    Times.push_back(t);

    dMatrix      Errors_History;
    dMatrix      Steps_History;
    wVec Failed_Trials;
    dVec        Error_Attempt(Params.max_trial, 0.0);
    dVec        Step_Attempt(Params.max_trial, 0.0);

    auto y_old = y;
    auto y_rk5 = y;
    auto y_mid = y;
    auto k1    = y;
    auto k2    = y;
    auto k3    = y;
    auto k4    = y;
    auto k5    = y;
    // Case 1: No error estimation - simple RK4
    if (!Params.error_estimate)
    {
        auto dt_half  = dt * 0.5;
        auto dt_third = dt / 3.0;
        auto dt_sixth = dt_third * 0.5;
        while (t < t1)
        {
            k1 = deriv(t, y);

            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_old);

            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt_half * k2[i];
            k3 = deriv(t + dt_half, y_old);

            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt * k3[i];
            k4 = deriv(t + dt, y_old);

            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt_third * (k2[i] + k3[i]) + dt_sixth * (k1[i] + k4[i]);
            
            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 2: Error estimation without variable steps
    else if (!Params.variable_steps)
    {
        auto dt_half  = dt * 0.5;
        auto dt_third = dt / 3.0;
        auto dt_sixth = dt_third * 0.5;
        auto dt_23    = dt_third * 2.0;
        auto coef_1   = dt * 7.0 / 90.0;
        auto coef_2   = dt * 32.0 / 90.0;
        auto coef_3   = dt * 12.0 / 90.0;
        auto dt_3     = dt * 3.0;
        auto dt_4     = dt * 4.0;
        while (t < t1)
        {
            // First stage
            k1 = deriv(t, y);

            // Second stage
            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_old);

            // Third stage
            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt_half * k2[i];
            k3 = deriv(t + dt_half, y_old);

            // Fourth stage
            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt * k3[i];
            k4 = deriv(t + dt, y_old);

            // Fifth stage
            y_old = y;
            for (size_t i = 0; i < y_old.size(); ++i)
                y_old[i] += dt * k1[i] - dt_3 * k2[i] + dt_4 * k3[i];
            k5 = deriv(t + dt_23, y_old);

            // Update using RK5
            for (size_t i = 0; i < y_rk5.size(); ++i)
                y_rk5[i] += coef_1 * k1[i] + coef_2 * (k3[i] + k5[i])  + coef_3 * k4[i];

            // Compute RK4 solution (using k1, k2, k3, k4 from RK5)
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt_third * (k2[i] + k3[i]) + dt_sixth * (k1[i] + k4[i]);

            // Calculate error
            auto error = estimate_error(y, y_rk5, Params);
            Errors.push_back(error);
            Steps.push_back(dt);

            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 3: Full adaptive step size with error estimation
    else
    {
        while (t < t1)
        {
            bool   step_accepted = false;
            y_old                = y;
            auto   current_dt    = dt;
            size_t trials        = 0;

            while (!step_accepted)
            {
                y = y_old;  // Reset to previous state

                // Compute RK5 solution first (more efficient)
                y_rk5         = y_old;
                auto dt_half  = current_dt * 0.5;
                auto dt_third = current_dt / 3.0;
                auto dt_sixth = dt_third * 0.5;
                auto dt_23    = dt_third * 2.0;
                auto coef_1   = current_dt * 7.0 / 90.0;
                auto coef_2   = current_dt * 32.0 / 90.0;
                auto coef_3   = current_dt * 12.0 / 90.0;
                auto dt_3     = current_dt * 3.0;
                auto dt_4     = current_dt * 4.0;
                // First stage
                k1 = deriv(t, y_old);

                // Second stage
                y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += dt_half * k1[i];
                k2 = deriv(t + dt_half, y_mid);

                // Third stage
                y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += dt_half * k2[i];
                k3 = deriv(t + dt_half, y_mid);

                // Fourth stage
                y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += current_dt * k3[i];
                k4 = deriv(t + current_dt, y_mid);

                // Fifth stage
                y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += current_dt * k1[i] - dt_3 * k2[i] + dt_4 * k3[i];
                k5 = deriv(t + dt_23, y_mid);

                // Update using RK5
                for (size_t i = 0; i < y_rk5.size(); ++i)
                    y_rk5[i] += coef_1 * k1[i] + coef_2 * (k3[i] + k5[i]) + coef_3 * k4[i];

                // Compute RK4 solution (using k1, k2, k3, k4 from RK5)
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] += dt_third * (k2[i] + k3[i]) + dt_sixth * (k1[i] + k4[i]);

                // Calculate error and adjust step size
                auto error = estimate_error(y, y_rk5, Params);
                Error_Attempt[trials] = error;
                Step_Attempt[trials] = current_dt;
                trials++;

                auto new_dt = adjust_step_size(current_dt, error, Params, 4);
                if (new_dt != current_dt)
                {
                    current_dt = new_dt;
                    continue;
                }

                step_accepted = true;
                dt = new_dt;
            }

            if (Params.attempts_history)
            {
                Errors_History.push_back(Error_Attempt);
                Steps_History.push_back(Step_Attempt);
                Failed_Trials.push_back(trials);
            }

            Errors.push_back(Error_Attempt[trials - 1]);
            Steps.push_back(current_dt);
            t += current_dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }

    SolverResults Results;
    Results.solution = result;
    Results.time_points = Times;
    if (Params.error_estimate)
    {
        Results.errors = Errors;
        if (Params.variable_steps)
        {
            Results.step_sizes = Steps;
            if (Params.attempts_history)
            {
                Results.errors_history = Errors_History;
                Results.steps_history = Steps_History;
                Results.failed_trials = Failed_Trials;
            }
        }
    }
    return Results;
}
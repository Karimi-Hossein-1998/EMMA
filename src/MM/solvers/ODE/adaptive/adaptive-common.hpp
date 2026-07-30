#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../utility/adaptive-solvers-utility.hpp"

// Heun-Euler RK1(2) solver
SolverResults heun_euler_solver(const SolverParameters& Params)
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
    auto y_heun = y;
    auto y_old  = y;
    auto k1     = y;
    auto k2     = y;
    // Case 1: No error estimation - simple Euler
    if (!Params.error_estimate)
    {
        while (t < t1)
        {
            k1 = deriv(t, y);
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
        y_old        = y;
        y_heun       = y_old;
        auto dt_half = dt * 0.5;
        while (t < t1)
        {
            // Compute Heun solution (2nd order)
            k1 = deriv(t, y);
            
            // Second stage
            y_heun = y;
            for (size_t i = 0; i < y_heun.size(); ++i)
                y_heun[i] += dt * k1[i];
            k2 = deriv(t + dt, y_heun);

            // Update using Heun
            for (size_t i = 0; i < y_heun.size(); ++i)
                y_heun[i] = y[i] + dt_half * (k1[i] + k2[i]);

            // Compute Euler solution (1st order)
            for (size_t i = 0; i < y.size(); ++i)
                y[i] = y[i] + dt * k1[i];

            // Calculate error
            auto error = estimate_error(y, y_heun, Params);
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

            while (!step_accepted && trials < Params.max_trial)
            {
                y            = y_old;  // Reset to previous state
                y_heun       = y_old;
                auto dt_half = current_dt * 0.5;
                // Compute Heun solution (2nd order)
                k1 = deriv(t, y_old);
                
                // Second stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += current_dt * k1[i];
                k2 = deriv(t + current_dt, y_temp);

                // Update using Heun
                for (size_t i = 0; i < y_heun.size(); ++i)
                    y_heun[i] = y_old[i] + dt_half * (k1[i] + k2[i]);

                // Compute Euler solution (1st order)
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] = y_old[i] + current_dt * k1[i];

                // Calculate error and adjust step size
                auto error = estimate_error(y, y_heun, Params);
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

template<typename... Args>
auto rk12_solver(Args&&... args) -> decltype(rk4_solver(std::forward<Args>(args)...))
{
    return heun_euler_solver(std::forward<Args>(args)...);
}

// Bogacki-Shampine RK2(3) solver
SolverResults bogacki_shampine_solver(const SolverParameters& Params)
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
    auto y_rk3  = y;
    auto y_old  = y;
    auto k1     = y;
    auto k2     = y;
    auto k3     = y;
    // Case 1: No error estimation - simple RK2
    if (!Params.error_estimate)
    {
        auto dt_half = dt * 0.5;
        while (t < t1)
        {
            k1 = deriv(t, y);
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_temp);
            
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
        y_rk3         = y;
        auto dt_half  = dt * 0.5; 
        auto dt_34    = dt * 0.75;
        auto dt_14    = dt * 0.25;
        auto dt_third = dt / 3.0;
        auto dt_29    = dt * 2.0 / 9.0;
        auto dt_49    = dt * 4.0 / 9.0;
        while (t < t1)
        {
            // Compute RK3 solution
            k1 = deriv(t, y);
            
            // Second stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_temp);

            // Third stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_34 * k2[i];
            k3 = deriv(t + dt_34, y_temp);

            // Update using RK3
            for (size_t i = 0; i < y_rk3.size(); ++i)
                y_rk3[i] = y[i] + dt_29 * k1[i] + dt_third * k2[i] + dt_49 * k3[i];

            // Compute RK2 solution
            for (size_t i = 0; i < y.size(); ++i)
                y[i] = y[i] + dt_14 * k1[i] + dt_34 * k3[i];

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
                auto dt_half  = current_dt * 0.5; 
                auto dt_34    = current_dt * 0.75;
                auto dt_14    = current_dt * 0.25;
                auto dt_third = current_dt / 3.0;
                auto dt_29    = current_dt * 2.0 / 9.0;
                auto dt_49    = current_dt * 4.0 / 9.0;

                // Compute RK3 solution
                y_rk3 = y_old;
                k1 = deriv(t, y_old);
                
                // Second stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_half * k1[i];
                k2 = deriv(t + dt_half, y_temp);

                // Third stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_34 * k2[i];
                k3 = deriv(t + dt_34, y_temp);

                // Update using RK3
                for (size_t i = 0; i < y_rk3.size(); ++i)
                    y_rk3[i] = y_old[i] + dt_29 * k1[i] + dt_third * k2[i] + dt_49 * k3[i];

                // Compute RK2 solution
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] = y_old[i] + dt_14 * k1[i] + dt_34 * k3[i];

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

template<typename... Args>
auto rk23_solver(Args&&... args) -> decltype(rk4_solver(std::forward<Args>(args)...))
{
    return bogacki_shampine_solver(std::forward<Args>(args)...);
}

// Fehlberg RK3(4) solver
SolverResults fehlberg34_solver(const SolverParameters& Params)
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
        auto dt_3_9  = dt / 3.0;
        auto dt_4_9  = dt_2_9 * 2;
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
                y[i] += dt_2_9 * k1[i] + dt_3_9 * k2[i] + dt_4_9 * k3[i];

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
        auto dt_1_4   = dt * 0.25;
        auto dt_2_9   = dt * 2.0 / 9.0;
        auto dt_1_3   = dt / 3.0;
        auto dt_4_9   = dt * 4.0 / 9.0;
        auto dt_7_24  = dt * 7.0 / 24.0;
        auto dt_1_4b  = dt * 0.25;
        auto dt_1_3b  = dt / 3.0;
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
                y_rk4[i] = y[i] + dt_7_24 * k1[i] + dt_1_4b * k2[i] + dt_1_3b * k3[i] + dt_1_8 * k4[i];

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
                auto dt_1_4   = current_dt * 0.25;
                auto dt_2_9   = current_dt * 2.0 / 9.0;
                auto dt_1_3   = current_dt / 3.0;
                auto dt_4_9   = current_dt * 4.0 / 9.0;
                auto dt_7_24  = current_dt * 7.0 / 24.0;
                auto dt_1_4b  = current_dt * 0.25;
                auto dt_1_3b  = current_dt / 3.0;
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
                    y_rk4[i] = y_old[i] + dt_7_24 * k1[i] + dt_1_4b * k2[i] + dt_1_3b * k3[i] + dt_1_8 * k4[i];

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

template<typename... Args>
auto rk34_solver(Args&&... args) -> decltype(rk4_solver(std::forward<Args>(args)...))
{
    return fehlberg34_solver(std::forward<Args>(args)...);
}

// Dormand-Prince RK4(5) solver
SolverResults dormand_prince_solver(const SolverParameters& Params)
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
    auto y_rk5  = y;
    auto y_old  = y;
    auto k1     = y;
    auto k2     = y;
    auto k3     = y;
    auto k4     = y;
    auto k5     = y;
    auto k6     = y;
    auto k7     = y;
    // Case 1: No error estimation - simple RK4
    if (!Params.error_estimate)
    {
        auto dt_half  = dt * 0.5;
        auto dt_third = dt / 3.0;
        auto dt_sixth = dt_third * 0.5;
        while (t < t1)
        {
            k1 = deriv(t, y);
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k1[i];
            k2 = deriv(t + dt_half, y_temp);
            
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k2[i];
            k3 = deriv(t + dt_half, y_temp);
            
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt * k3[i];
            k4 = deriv(t + dt, y_temp);
            
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt_third * (k1[i] + k4[i]) + dt_sixth * (k2[i] + k3[i]);
            
            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 2: Error estimation without variable steps
    else if (!Params.variable_steps)
    {
        auto dt_15     = dt / 5.0;
        auto dt_310    = dt * 0.3;
        auto dt_340    = dt * 0.075;
        auto dt_940    = dt_340 * 3.0;
        auto dt_4445   = dt * 44.0 / 45.0;
        auto dt_5615   = dt * 56.0 / 15.0;
        auto dt_329    = dt * 32.0 / 9.0;
        auto dt_45     = dt * 0.8; 
        auto coef1_k5  = dt * 19372.0 / 6561.0;
        auto coef2_k5  = dt * 25360.0 / 2187.0;
        auto coef3_k5  = dt * 64448.0 / 6561.0;
        auto coef4_k5  = dt * 212.0 / 729.0;
        auto dt_89     = dt * 8.0 / 9.0;
        auto coef1_k6  = dt * 9017.0 / 3168.0;
        auto coef2_k6  = dt * 355.0 / 33.0;
        auto coef3_k6  = dt * 46732.0 / 5247.0;
        auto coef4_k6  = dt * 49.0 / 176.0;
        auto coef5_k6  = dt * 5103.0 / 18656.0;
        auto coef1_k7  = dt * 35.0 / 384.0;
        auto coef2_k7  = dt * 500.0 / 1113.0;
        auto coef3_k7  = dt * 125.0 / 192.0;
        auto coef4_k7  = dt * 2187.0 / 6784.0;
        auto coef5_k7  = dt * 11.0 / 84.0;
        auto coef1_rk4 = dt * 5179.0 / 57600.0;
        auto coef2_rk4 = dt * 7571.0 / 16695.0;
        auto coef3_rk4 = dt * 393.0 / 640.0;
        auto coef4_rk4 = dt * 92097.0 / 339200.0;
        auto coef5_rk4 = dt * 187.0 / 2100.0;
        auto coef6_rk4 = dt / 40.0;
        while (t < t1)
        {
            // Compute RK5 solution
            k1 = deriv(t, y);
            
            // Second stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_15 * k1[i];
            k2 = deriv(t + dt_15, y_temp);

            // Third stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_340 * k1[i] + dt_940 * k2[i];
            k3 = deriv(t + dt_310, y_temp);

            // Fourth stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_4445 * k1[i] - dt_5615 * k2[i] + dt_329 * k3[i];
            k4 = deriv(t + dt_45, y_temp);

            // Fifth stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += coef1_k5 * k1[i] - coef2_k5 * k2[i] + coef3_k5 * k3[i] - coef4_k5 * k4[i];
            k5 = deriv(t + dt_89, y_temp);

            // Sixth stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += coef1_k6 * k1[i] - coef2_k6 * k2[i] + coef3_k6 * k3[i] + coef4_k6 * k4[i] - coef5_k6 * k5[i];
            k6 = deriv(t + dt, y_temp);

            // Update using RK5
            for (size_t i = 0; i < y_rk5.size(); ++i)
            {
                y_temp[i] = coef1_k7 * k1[i] + coef2_k7 * k3[i] + coef3_k7 * k4[i] - coef4_k7 * k5[i] + coef5_k7 * k6[i];
                y_rk5[i] = y[i] + y_temp[i];
            }

            // Seventh stage
            k7 = deriv(t + dt, y_rk5);

            // Compute RK4 solution
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += coef1_rk4 * k1[i] + coef2_rk4 * k3[i] + coef3_rk4 * k4[i] - coef4_rk4 * k5[i] + coef5_rk4 * k6[i] + coef6_rk4 * k7[i];

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
            auto   y_old         = y;
            double current_dt    = dt;
            size_t trials        = 0;

            while (!step_accepted)
            {
                auto dt_15     = current_dt / 5.0;
                auto dt_310    = current_dt * 0.3;
                auto dt_340    = current_dt * 0.075;
                auto dt_940    = dt_340 * 3.0;
                auto dt_4445   = current_dt * 44.0 / 45.0;
                auto dt_5615   = current_dt * 56.0 / 15.0;
                auto dt_329    = current_dt * 32.0 / 9.0;
                auto dt_45     = current_dt * 0.8; 
                auto coef1_k5  = current_dt * 19372.0 / 6561.0;
                auto coef2_k5  = current_dt * 25360.0 / 2187.0;
                auto coef3_k5  = current_dt * 64448.0 / 6561.0;
                auto coef4_k5  = current_dt * 212.0 / 729.0;
                auto dt_89     = current_dt * 8.0 / 9.0;
                auto coef1_k6  = current_dt * 9017.0 / 3168.0;
                auto coef2_k6  = current_dt * 355.0 / 33.0;
                auto coef3_k6  = current_dt * 46732.0 / 5247.0;
                auto coef4_k6  = current_dt * 49.0 / 176.0;
                auto coef5_k6  = current_dt * 5103.0 / 18656.0;
                auto coef1_k7  = current_dt * 35.0 / 384.0;
                auto coef2_k7  = current_dt * 500.0 / 1113.0;
                auto coef3_k7  = current_dt * 125.0 / 192.0;
                auto coef4_k7  = current_dt * 2187.0 / 6784.0;
                auto coef5_k7  = current_dt * 11.0 / 84.0;
                auto coef1_rk4 = current_dt * 5179.0 / 57600.0;
                auto coef2_rk4 = current_dt * 7571.0 / 16695.0;
                auto coef3_rk4 = current_dt * 393.0 / 640.0;
                auto coef4_rk4 = current_dt * 92097.0 / 339200.0;
                auto coef5_rk4 = current_dt * 187.0 / 2100.0;
                auto coef6_rk4 = current_dt / 40.0;

                // Compute RK5 solution
                y     = y_old;  // Reset to previous state
                y_rk5 = y_old;
                k1    = deriv(t, y_old);

                // Second stage
                auto y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_15 * k1[i];
                k2 = deriv(t + dt_15, y_temp);

                // Third stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_340 * k1[i] + dt_940 * k2[i];
                k3 = deriv(t + dt_310, y_temp);

                // Fourth stage
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_4445 * k1[i] - dt_5615 * k2[i] + dt_329 * k3[i];
                k4 = deriv(t + dt_45, y_temp);

                // Fifth stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += coef1_k5 * k1[i] - coef2_k5 * k2[i] + coef3_k5 * k3[i] - coef4_k5 * k4[i];
                k5 = deriv(t + dt_89, y_temp);

                // Sixth stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += coef1_k6 * k1[i] - coef2_k6 * k2[i] + coef3_k6 * k3[i] + coef4_k6 * k4[i] - coef5_k6 * k5[i];
                k6 = deriv(t + dt, y_temp);

                // Update using RK5
                y_temp = y_old;
                for (size_t i = 0; i < y_rk5.size(); ++i)
                {
                    y_temp[i] = coef1_k7 * k1[i] + coef2_k7 * k3[i] + coef3_k7 * k4[i] - coef4_k7 * k5[i] + coef5_k7 * k6[i];
                    y_rk5[i] = y[i] + y_temp[i];
                }

                // Seventh stage
                k7 = deriv(t + dt, y_rk5);

                // Compute RK4 solution
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] = y_old[i] + coef1_rk4 * k1[i] + coef2_rk4 * k3[i] + coef3_rk4 * k4[i] - coef4_rk4 * k5[i] + coef5_rk4 * k6[i] + coef6_rk4 * k7[i];

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

template<typename... Args>
auto rk45_solver(Args&&... args) -> decltype(rk4_solver(std::forward<Args>(args)...))
{
    return dormand_prince_solver(std::forward<Args>(args)...);
}

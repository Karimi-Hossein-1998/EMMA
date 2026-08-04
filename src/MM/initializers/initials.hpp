#pragma once
#include "../typedefs/header.hpp"

namespace MathEngine
{ // namespace MathEngine
// Uniform distribution in [min, max)
inline dVec  random_uniform(
    size_t   N,
    double   min_val,
    double   max_val,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) 
    {
        throw std::invalid_argument("[random_uniform_phases] Number of phases N cannot be zero.");
    }
    if (min_val > max_val) 
    {
        throw std::invalid_argument("[random_uniform_phases] Min value (" + std::to_string(min_val) + ") cannot be greater than max value (" + std::to_string(max_val) + ").");
    }

    // Original logic
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(min_val, max_val);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Normal (Gaussian) distribution
inline dVec  random_normal(
    size_t      N,
    double      mean,
    double      stddev,
    unsigned    seed
)
{
    // Error Handling
    if (N == 0) 
    {
        throw std::invalid_argument("[random_normal_phases] Number of phases N cannot be zero.");
    }
    if (stddev < 0.0) 
    {
        throw std::invalid_argument("[random_normal_phases] Standard deviation stddev (" + std::to_string(stddev) + ") cannot be negative.");
    }

    // Original logic
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist(mean, stddev);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Cauchy (Lorentzian) distribution
inline dVec  random_cauchy(
    size_t      N,
    double      location,
    double      scale,
    unsigned    seed
)
{
    // Error Handling
    if (N == 0) 
    {
        throw std::invalid_argument("[random_cauchy_phases] Number of phases N cannot be zero.");
    }
    if (scale <= 0.0) 
    {
        throw std::invalid_argument("[random_cauchy_phases] Scale parameter (" + std::to_string(scale) + ") must be positive.");
    }

    // Original logic
    std::mt19937 rng(seed);
    std::cauchy_distribution<double> dist(location, scale);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Exponential distribution
inline dVec  random_exponential(
    size_t      N,
    double      lambda, // rate parameter
    unsigned    seed
)
{
    // Error Handling
    if (N == 0) 
    {
        throw std::invalid_argument("[random_exponential_phases] Number of phases N cannot be zero.");
    }
    if (lambda <= 0.0) 
    {
        throw std::invalid_argument("[random_exponential_phases] Lambda (rate parameter) (" + std::to_string(lambda) + ") must be positive.");
    }

    // Original logic
    std::mt19937 rng(seed);
    std::exponential_distribution<double> dist(lambda);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Uniform distribution on the unit circle [0, 2*pi)
inline dVec  random_circle(
    size_t      N,
    unsigned    seed
)
{
    // Error Handling
    if (N == 0) 
    {
        throw std::invalid_argument("[random_circle_phases] Number of phases N cannot be zero.");
    }
    return random_uniform(N, -PI, PI, seed);
}

// Splay phases: equidistant around the circle [0, 2*pi)
inline dVec  splay(
    size_t N
)
{
    // Error Handling
    if (N == 0) 
    {
        throw std::invalid_argument("[splay_phases] Number of phases N cannot be zero.");
    }

    dVec  phases(N);
    if (N == 1) 
    {
        phases[0] = 0.0;
        return phases;
    }

    double delta = 2.0 * PI / static_cast<double>(N);
    for (size_t i = 0; i < N; ++i)
        phases[i] = i * delta;
    return phases;
}

// Splay phases with random perturbation in [-amplitude, amplitude]
inline dVec  splay_perturbed(
    size_t   N,
    double   amplitude,
    unsigned seed
)
{
    // Error Handling
    if (amplitude < 0.0) 
    {
        throw std::invalid_argument("[splay_phases_perturbed] Amplitude (" + std::to_string(amplitude) + ") cannot be negative.");
    }

    auto phases = splay(N);
    if (N == 0) return phases;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-amplitude, amplitude);
    for (auto& x : phases)
        x += dist(rng);
    return phases;
}

// Generate phases/frequencies for one module based on a given condition (uniform, normal, etc)
inline dVec  module_by_condition(
    size_t                module_size,
    const std::string&    dist_type,
    double                a,
    double                b,
    unsigned              seed
)
{
    // Error Handling
    if (module_size == 0) 
    {
        throw std::invalid_argument("[module_by_condition] module_size cannot be zero for distribution type '" + dist_type + "'.");
    }

    if (dist_type == "uniform")
    {
        return random_uniform(module_size, a, b, seed);
    }
    else if (dist_type == "normal")
    {
        return random_normal(module_size, a, b, seed);
    }
    else if (dist_type == "cauchy")
    {
        return random_cauchy(module_size, a, b, seed);
    }
    else if (dist_type == "exponential")
    {
        return random_exponential(module_size, a, seed);
    }
    else if (dist_type == "circle")
    {
        return random_circle(module_size, seed);
    }
    else if (dist_type == "splay")
    {
        return splay(module_size);
    }
    else if (dist_type == "splay_perturbed")
    {
        return splay_perturbed(module_size,a,seed);
    }
    else
    {
        throw std::invalid_argument("[module_by_condition] Unknown distribution type: '" + dist_type + "'.");
    }
}

// Generate phases/frequencies for one module and copy to all modules (identical modules)
inline dVec  identical_modules(
    size_t               N_per_module,
    size_t               num_modules,
    const std::string&   dist_type,
    double               a,
    double               b,
    unsigned             seed
)
{
    // Error Handling
    if (N_per_module == 0) 
    {
        throw std::invalid_argument("[identical_modules] N_per_module (nodes per module) cannot be zero.");
    }
    if (num_modules == 0) 
    {
        throw std::invalid_argument("[identical_modules] num_modules cannot be zero.");
    }

    dVec  base = module_by_condition(N_per_module, dist_type, a, b, seed);
    dVec  result(N_per_module * num_modules);
    for (size_t m = 0; m < num_modules; ++m)
    {
        for (size_t i = 0; i < N_per_module; ++i)
        {
            result[m * N_per_module + i] = base[i];
        }
    }
    return result;
}

inline dVec modules(
    size_t               N_per_module,
    size_t               num_modules,
    const std::string&   dist_type,
    double               a,
    double               b,
    unsigned             seed,
    bool                 identical
)
{
    if (N_per_module == 0)
    {
        throw std::invalid_argument("[identical_modules] N_per_module (nodes per module) cannot be zero.");
    }
    if (num_modules == 0)
    {
        throw std::invalid_argument("[identical_modules] num_modules cannot be zero.");
    }
    dVec results(N_per_module*num_modules);
    if (identical)
    {
    	results = identical_modules(N_per_module, num_modules, dist_type, a, b, seed);
    }
    else
    {
        dVec subResults(N_per_module);
		for (size_t mOdule=0; mOdule<num_modules; ++mOdule)
        {
            // seed += mOdule;
			subResults = module_by_condition(N_per_module,dist_type,a,b,seed+mOdule);
            for (size_t i=0; i<N_per_module; ++i)
            {
                results[i+mOdule*N_per_module]=subResults[i];
            }
        }
    }
	return results;
}
} // End namespace MathEngine

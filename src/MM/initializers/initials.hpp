#pragma once
#include "../typedefs/header.hpp"

namespace MathEngine
{ // namespace MathEngine
enum class InitState
{
    Uniform=0,
    Normal,
    Cauchy,
    Exponential,
    Circle,
    Splay,
    SplayPerturbed,
    Modules
};
enum class InitType
{
    Uniform=0,
    Normal,
    Cauchy,
    Exponential,
    Circle,
    Splay,
    SplayPerturbed
};
struct InitializerParams
{
    InitState initState  = InitState::Uniform;
    InitType  moduleType = InitType::Uniform;
    double    param1     = 0.0;
    double    param2     = 1.0;
    size_t    N          = 50;
    size_t    numModules = 1;
    size_t    moduleSize = 50;
    size_t    seed       = 41;
    bool      identical  = false;
};
// Uniform distribution in [min, max)
inline dVec random_uniform(
    size_t   N,
    double   min_val,
    double   max_val,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (min_val > max_val) std::swap(min_val,max_val);
    // Original logic
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(min_val, max_val);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Normal (Gaussian) distribution
inline dVec random_normal(
    size_t   N,
    double   mean,
    double   stddev,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (stddev <= 0.0) return dVec(N,mean);;
    // Original logic
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist(mean, stddev);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Cauchy (Lorentzian) distribution
inline dVec random_cauchy(
    size_t   N,
    double   location,
    double   scale,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (scale <= 0.0) return dVec(N,location);
    // Original logic
    std::mt19937 rng(seed);
    std::cauchy_distribution<double> dist(location, scale);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Exponential distribution
inline dVec random_exponential(
    size_t   N,
    double   lambda, // rate parameter
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (lambda <= 0.0) return dVec(N,0.0);
    // Original logic
    std::mt19937 rng(seed);
    std::exponential_distribution<double> dist(lambda);
    dVec  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Uniform distribution on the unit circle [-pi, pi)
inline dVec random_circle(
    size_t   N,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    return random_uniform(N, -PI, PI, seed);
}

// Splay phases: equidistant around the circle [0, 2*pi)
inline dVec splay(size_t N)
{
    // Error Handling
    if (N == 0) return {};
    if (N == 1) return {0.0};

    dVec   phases(N);
    double delta = 2.0 * PI / static_cast<double>(N);
    for (size_t i = 0; i < N; ++i)
        phases[i] = i * delta;
    return phases;
}

// Splay phases with random perturbation in [-amplitude, amplitude]
inline dVec splay_perturbed(
    size_t   N,
    double   amplitude,
    unsigned seed
)
{
    // Error Handling
    if (N==0) return {};
    if (std::abs(amplitude) <= 5e-15) return splay(N);

    auto phases = splay(N);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-amplitude, amplitude);
    for (auto& x : phases)
        x += dist(rng);
    return phases;
}

// Generate phases/frequencies for one module based on a given condition (uniform, normal, etc)
inline dVec module_by_condition(
    size_t             module_size,
    InitType           initType,
    double             a,
    double             b,
    unsigned           seed
)
{
    // Error Handling
    if (module_size == 0) return {};
    switch (initType)
    {
        case InitType::Uniform:        return random_uniform(module_size, a, b, seed);break;
        case InitType::Normal:         return random_normal(module_size,a,b,seed);break;
        case InitType::Cauchy:         return random_cauchy(module_size,a,b,seed);break;
        case InitType::Exponential:    return random_exponential(module_size,a,seed);break;
        case InitType::Circle:         return random_circle(module_size,seed);break;
        case InitType::Splay:          return splay(module_size);break;
        case InitType::SplayPerturbed: return splay_perturbed(module_size,a,seed);break;
        default:                       return random_uniform(module_size, a, b, seed);break;
    }
}

// Generate phases/frequencies for one module and copy to all modules (identical modules)
inline dVec identical_modules(
    size_t   N_per_module,
    size_t   num_modules,
    InitType initType,
    double   a,
    double   b,
    unsigned seed
)
{
    // Error Handling
    if (N_per_module == 0) return {};
    if (num_modules == 0) return {};

    dVec base = module_by_condition(N_per_module, initType, a, b, seed);
    dVec result(N_per_module * num_modules);
    for (size_t m = 0; m < num_modules; ++m)
        for (size_t i = 0; i < N_per_module; ++i)
            result[m * N_per_module + i] = base[i];
    return result;
}

inline dVec modules(
    size_t   N_per_module,
    size_t   num_modules,
    InitType initType,
    double   a,
    double   b,
    unsigned seed,
    bool     identical
)
{
    if (N_per_module == 0) return {};
    if (num_modules == 0) return {};
    if (identical)
    {
    	return identical_modules(N_per_module, num_modules, initType, a, b, seed);
    }
    dVec results(N_per_module*num_modules);
    dVec subResults(N_per_module);
    for (size_t mOdule=0; mOdule<num_modules; ++mOdule)
    {
        subResults = module_by_condition(N_per_module,initType,a,b,seed+mOdule);
        for (size_t i=0; i<N_per_module; ++i)
        {
            results[i+mOdule*N_per_module]=subResults[i];
        }
    }
    return results;
}

inline dVec initialize_vector(const InitializerParams& initParams)
{
    switch(initParams.initState)
    {
        case InitState::Uniform:        return random_uniform(initParams.N,initParams.param1,initParams.param2,initParams.seed);break;
        case InitState::Normal:         return random_normal(initParams.N,initParams.param1,initParams.param2,initParams.seed);break;
        case InitState::Cauchy:         return random_cauchy(initParams.N,initParams.param1,initParams.param2,initParams.seed);break;
        case InitState::Exponential:    return random_exponential(initParams.N,initParams.param1,initParams.seed);break;
        case InitState::Circle:         return random_circle(initParams.N,initParams.seed);break;
        case InitState::Splay:          return splay(initParams.N);break;
        case InitState::SplayPerturbed: return splay_perturbed(initParams.N,initParams.param1,initParams.seed);break;
        case InitState::Modules:        return modules(initParams.moduleSize,initParams.numModules,initParams.moduleType,
            									initParams.param1,initParams.param2,initParams.seed,initParams.identical);break;
        default:                        return random_uniform(initParams.N,initParams.param1,initParams.param2,initParams.seed); break;
    }
}
} // End namespace MathEngine

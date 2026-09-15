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

template <FPNumber Num>
struct InitializerParams
{
    InitState initState  = InitState::Uniform;
    InitType  moduleType = InitType::Uniform;
    Num       param1     = static_cast<Num>(0.0);
    Num       param2     = static_cast<Num>(1.0);
    size_t    N          = 50;
    size_t    numModules = 1;
    size_t    moduleSize = 50;
    size_t    seed       = 41;
    bool      identical  = false;
};
// Uniform distribution in [min, max)
template <FPNumber Num>
inline Vec<Num> random_uniform(
    size_t   N,
    Num      min_val,
    Num      max_val,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (min_val > max_val) std::swap(min_val,max_val);
    // Original logic
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Num> dist(min_val, max_val);
    Vec<Num>  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Normal (Gaussian) distribution
template <FPNumber Num>
inline Vec<Num> random_normal(
    size_t   N,
    Num      mean,
    Num      stddev,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (stddev <= 0.0) return Vec<Num>(N,mean);;
    // Original logic
    std::mt19937 rng(seed);
    std::normal_distribution<Num> dist(mean, stddev);
    Vec<Num>  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Cauchy (Lorentzian) distribution
template <FPNumber Num>
inline Vec<Num> random_cauchy(
    size_t   N,
    Num      location,
    Num      scale,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (scale <= 0.0) return Vec<Num>(N,location);
    // Original logic
    std::mt19937 rng(seed);
    std::cauchy_distribution<Num> dist(location, scale);
    Vec<Num>  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Exponential distribution
template <FPNumber Num>
inline Vec<Num> random_exponential(
    size_t   N,
    Num      lambda, // rate parameter
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    if (lambda <= 0.0) return Vec<Num>(N,0.0);
    // Original logic
    std::mt19937 rng(seed);
    std::exponential_distribution<Num> dist(lambda);
    Vec<Num>  phases(N);
    for (auto& x : phases)
        x = dist(rng);
    return phases;
}

// Uniform distribution on the unit circle [-pi, pi)
template <FPNumber Num>
inline Vec<Num> random_circle(
    size_t   N,
    unsigned seed
)
{
    // Error Handling
    if (N == 0) return {};
    return random_uniform(N, -PI, PI, seed);
}

// Splay phases: equidistant around the circle [0, 2*pi)
template <FPNumber Num>
inline Vec<Num> splay(size_t N)
{
    // Error Handling
    if (N == 0) return {};
    if (N == 1) return {Num{}};

    Vec<Num> phases(N);
    Num delta = static_cast<Num>(2.0) * PI / static_cast<Num>(N);
    for (size_t i = 0; i < N; ++i)
        phases[i] = i * delta;
    return phases;
}

// Splay phases with random perturbation in [-amplitude, amplitude]
template <FPNumber Num>
inline Vec<Num> splay_perturbed(
    size_t   N,
    Num      amplitude,
    unsigned seed
)
{
    // Error Handling
    if (N==0) return {};

    auto phases = splay<Num>(N);
    if (std::abs(amplitude) <= std::numeric_limits<Num>::epsilon()) return phases;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Num> dist(-amplitude, amplitude);
    for (auto& x : phases)
        x += dist(rng);
    return phases;
}

// Generate phases/frequencies for one module based on a given condition (uniform, normal, etc)
template <FPNumber Num>
inline Vec<Num> module_by_condition(
    size_t   module_size,
    InitType initType,
    Num      a,
    Num      b,
    unsigned seed
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
        case InitType::Circle:         return random_circle<Num>(module_size,seed);break;
        case InitType::Splay:          return splay<Num>(module_size);break;
        case InitType::SplayPerturbed: return splay_perturbed(module_size,a,seed);break;
        default:                       return random_uniform(module_size, a, b, seed);break;
    }
}

// Generate phases/frequencies for one module and copy to all modules (identical modules)
template <FPNumber Num>
inline Vec<Num> identical_modules(
    size_t   N_per_module,
    size_t   num_modules,
    InitType initType,
    Num      a,
    Num      b,
    unsigned seed
)
{
    // Error Handling
    if (N_per_module == 0) return {};
    if (num_modules == 0) return {};

    Vec<Num> base = module_by_condition(N_per_module, initType, a, b, seed);
    Vec<Num> result(N_per_module * num_modules);
    for (size_t m = 0; m < num_modules; ++m)
        for (size_t i = 0; i < N_per_module; ++i)
            result[m * N_per_module + i] = base[i];
    return result;
}

template <FPNumber Num>
inline Vec<Num> modules(
    size_t   N_per_module,
    size_t   num_modules,
    InitType initType,
    Num      a,
    Num      b,
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
    Vec<Num> results(N_per_module*num_modules);
    Vec<Num> subResults(N_per_module);
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

template <FPNumber Num>
inline Vec<Num> initialize_vector(const InitializerParams<Num>& initParams)
{
    switch(initParams.initState)
    {
        case InitState::Uniform:        return random_uniform(initParams.N,initParams.param1,initParams.param2,initParams.seed);break;
        case InitState::Normal:         return random_normal(initParams.N,initParams.param1,initParams.param2,initParams.seed);break;
        case InitState::Cauchy:         return random_cauchy(initParams.N,initParams.param1,initParams.param2,initParams.seed);break;
        case InitState::Exponential:    return random_exponential(initParams.N,initParams.param1,initParams.seed);break;
        case InitState::Circle:         return random_circle<Num>(initParams.N,initParams.seed);break;
        case InitState::Splay:          return splay<Num>(initParams.N);break;
        case InitState::SplayPerturbed: return splay_perturbed(initParams.N,initParams.param1,initParams.seed);break;
        case InitState::Modules:        return modules(initParams.moduleSize,initParams.numModules,initParams.moduleType,
            									initParams.param1,initParams.param2,initParams.seed,initParams.identical);break;
        default:                        return random_uniform(initParams.N,initParams.param1,initParams.param2,initParams.seed); break;
    }
}
} // End namespace MathEngine

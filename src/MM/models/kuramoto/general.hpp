#pragma once
#include "../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
// General Kuramoto model with phase-lag and flexible dMatrix
inline dVec  kuramoto_general(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const dMatrix&       adj,
    double              alpha
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    // Original logic
    dVec  dtheta_dt = dVec( N, 0.0 );
    auto  k_norm    = K / static_cast<double>( N ); // Renamed k to k_norm for clarity
    for ( size_t i = 0; i < N; ++i )
    {
        double sum = 0.0;
        for ( size_t j = 0; j < N; ++j )
        {
            sum += adj[i][j] * std::sin( theta[j] - theta[i] - alpha );
        }
        dtheta_dt[i] = omega[i] + k_norm * sum;
    }
    return dtheta_dt;
}

// General Kuramoto model with phase-lag and flexible dMatrix (parallel)
inline dVec  kuramoto_general_parallel(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const dMatrix&       adj,
    double              alpha
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    std::vector<std::jthread> threads;
    size_t num_threads = std::min(N, static_cast<size_t>(std::max(1u, std::thread::hardware_concurrency())));
    if ( num_threads == 0 ) num_threads = 1;
    size_t chunk_size  = N / num_threads;
    auto   k_norm      = K / static_cast<double>( N ); // Renamed k = K/N to k_norm for clarity
    dVec   dtheta_dt   = dVec( N, 0.0 );
    
    for (size_t t = 0; t < num_threads; ++t) 
    {
        threads.emplace_back([&, t]() 
        {
            size_t start = t * chunk_size;
            size_t end = (t == num_threads - 1) ? N : (t + 1) * chunk_size;
            
            for (size_t i = start; i < end; ++i) 
            {
                double sum = 0.0;
                for (size_t j = 0; j < N; ++j) 
                {
                    sum += adj[i][j] * std::sin(theta[j] - theta[i] - alpha);
                }
                dtheta_dt[i] = omega[i] + k_norm * sum;
            }
        });
    }
    // No need to join, jthread automatically joins in destructor
    return dtheta_dt;
} 

// ======================================== //
//                                          //
// ------------ Wrappers ------------------ //
//                                          //
// ======================================== //
//
// Usage: 
// KuramotoParams params(N, alpha);
// params.K     = 1.0;
// params.omega = splay(N);
// params.adj   = random(N)
// auto kuramoto_func = kuramoto_general_parallel_wrapper(t, theta, params);
// ======================================== //
//                                          //
// ------------ Wrappers ------------------ //
//                                          //
// ======================================== //
//
// Usage: 
// KuramotoParams params(N, alpha);
// params.K     = 1.0;
// params.omega = splay(N);
// params.adj   = random(N)
// auto kuramoto_func = kuramoto_general_parallel_wrapper(t, theta, params);
//
// -------------------------------------------------------------------------
struct KuramotoParams
{
    // Model specific parameters in Kuramoto model (time and phases/theta are dynamics bounds)
    double K;     // Coupling strength
    int    N;     // Number of oscillators
    dVec   omega; // Natural frequencies
    dMatrix adj;   
    double alpha;

    KuramotoParams(int num_oscillators=50, double phae_lag = 0.0)
        : N(num_oscillators),
          omega(dVec(num_oscillators, 0.0)),
          adj(dMatrix(num_oscillators, dVec(num_oscillators, 0.0))),
          alpha(phae_lag)
    {}
};

inline MyFunc kuramoto_general_wrapper(const KuramotoParams& params)
{
    return [params](double t, const dVec& theta) -> dVec
    {
        return kuramoto_general(t, theta, params.omega, params.K, params.adj, params.alpha);
    };
}

inline MyFunc kuramoto_general_parallel_wrapper(const KuramotoParams& params)
{
    return [params](double t, const dVec& theta) -> dVec
    {
        return kuramoto_general_parallel(t, theta, params.omega, params.K, params.adj, params.alpha);
    };
}
} // End namespace MathEngine

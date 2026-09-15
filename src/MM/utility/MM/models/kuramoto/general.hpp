#pragma once
#include "../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
// General Kuramoto model with phase-lag and flexible Matrix<Num>
template <FPNumber Num>
inline void kuramoto_general(
    Num                time,
    const Vec<Num>&    theta,
    Vec<Num>&          dthetadt,
    const Vec<Num>&    omega,
    Num                K,
    const Matrix<Num>& adj,
    Num                alpha
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    // Original logic
    Num k_norm    = K / static_cast<Num>( N ); // Renamed k to k_norm for clarity
    for ( size_t i = 0; i < N; ++i )
    {
        Num sum = 0.0;
        for ( size_t j = 0; j < N; ++j )
        {
            sum += adj[i,j] * std::sin( theta[j] - theta[i] - alpha );
        }
        dthetadt[i] = omega[i] + k_norm * sum;
    }
    // return dthetadt;
}

// // General Kuramoto model with phase-lag and flexible Matrix<Num> (parallel)
// template <FPNumber Num>
// inline void kuramoto_general_parallel(
//     Num                time,
//     const Vec<Num>&    theta,
//     Vec<Num>&          dthetadt,
//     const Vec<Num>&    omega,
//     Num                K,
//     const Matrix<Num>& adj,
//     Num                alpha
// )
// {
//     const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

//     Vec<std::thread> threads;
//     size_t num_threads = std::min(N, static_cast<size_t>(std::max(1u, std::thread::hardware_concurrency())));
//     if ( num_threads == 0 ) num_threads = 1;
//     size_t chunk_size  = N / num_threads;
//     Num k_norm      = K / static_cast<Num>( N ); // Renamed k = K/N to k_norm for clarity

//     for (size_t t = 0; t < num_threads; ++t) 
//     {
//         threads.emplace_back([&, t]() 
//         {
//             size_t start = t * chunk_size;
//             size_t end = (t == num_threads - 1) ? N : (t + 1) * chunk_size;
            
//             for (size_t i = start; i < end; ++i) 
//             {
//                 Num sum = 0.0;
//                 for (size_t j = 0; j < N; ++j) 
//                 {
//                     sum += adj[i][j] * std::sin(theta[j] - theta[i] - alpha);
//                 }
//                 dthetadt[i] = omega[i] + k_norm * sum;
//             }
//         });
//     }
//     for (auto& t : threads) t.join();
//     // return dthetadt;
// } 

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
template <FPNumber Num>
struct KuramotoParams
{
    // Model specific parameters in Kuramoto model (time and phases/theta are dynamics bounds)
    Matrix<Num> adj;
    Vec<Num>    omega; // Natural frequencies
    Num         alpha;
    Num         K;     // Coupling strength
    int         N;     // Number of oscillators

    // KuramotoParams(int num_oscillators=50, Num phae_lag = 0.0)
    //     : N(num_oscillators),
    //       omega(Vec<Num>(num_oscillators, Num{})),
    //       adj(Matrix<Num>(num_oscillators, num_oscillators)),
    //       alpha(phae_lag)
    // {}
};

template <FPNumber Num>
inline MyFunc<Num> auto kuramoto_general_wrapper(const KuramotoParams<Num>& params)
{
    return [params](Num t, const Vec<Num>& theta, Vec<Num>& dthetadt) -> void
    {
        return kuramoto_general(t, theta, dthetadt, params.omega, params.K, params.adj, params.alpha);
    };
}
// template <FPNumber Num>
// inline MyFunc<Num> auto kuramoto_general_parallel_wrapper(const KuramotoParams<Num>& params)
// {
//     return [params](Num t, const Vec<Num>& theta, Vec<Num>& dthetadt) -> void
//     {
//         return kuramoto_general_parallel(t, theta, dthetadt, params.omega, params.K, params.adj, params.alpha);
//     };
// }
} // End namespace MathEngine

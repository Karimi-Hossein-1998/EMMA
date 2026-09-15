#pragma once
#include "../../typedefs/header.hpp"

namespace MathEngine
{ // namespace MathEngine
// Special modular Kuramoto model (ordinary)
template <FPNumber Num>
inline void kuramoto_special_modular(
    Num             time,
    const Vec<Num>& theta,
    Vec<Num>&       dthetadt,
    const Vec<Num>& omega,
    Num             intra_K,
    Num             inter_K,
    Num             alpha,
    const size_t    module_size
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    Num intra_k = intra_K / static_cast<Num>( N ); // Renamed k = K/N to k_norm for clarity
    Num inter_k = inter_K / static_cast<Num>( N ); // Renamed k = K/N to k_norm for clarity
    for (size_t i = 0; i < N; ++i) 
    {
        size_t module_i = i / module_size;
        Num sum = 0.0;
        for (size_t j = 0; j < N; ++j) 
        {
            size_t module_j = j / module_size;
            if ( i != j )
            {
                if ( module_i == module_j )
                {
                    sum += intra_k * std::sin(theta[j] - theta[i] - alpha);
                }
                else
                {
                    sum += inter_k * std::sin(theta[j] - theta[i] - alpha);
                }
            }
        }
        dthetadt[i] = omega[i] + sum;
    }
    // return dtheta_dt;
}

// // Special modular Kuramoto model (parallel)
// inline void kuramoto_special_modular_parallel(
//     Num              time,
//     const Vec<Num>&         theta,
//     Vec<Num>&               dthetadt,
//     const Vec<Num>&         omega,
//     Num              intra_K,
//     Num              inter_K,
//     Num              alpha,
//     const size_t        module_size
// )
// {
//     const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

//     std::vector<std::thread> threads;
//     size_t num_threads = std::min(N, static_cast<size_t>(std::max(1u, std::thread::hardware_concurrency())));
//     if ( num_threads == 0 ) num_threads = 1;
//     Num intra_k    = intra_K / static_cast<Num>( N ); // Renamed k = K/N to k_norm for clarity
//     Num inter_k    = inter_K / static_cast<Num>( N ); // Renamed k = K/N to k_norm for clarity
//     Vec<Num>   dtheta_dt  = Vec<Num>( N, 0.0 );
//     size_t chunk_size = N / num_threads;

//     for (size_t thread_index = 0; thread_index < num_threads; ++thread_index) 
//     {
//         threads.emplace_back([&, thread_index]() 
//         {
//             size_t start = thread_index * chunk_size;
//             size_t end   = (thread_index == num_threads - 1) ? N : (thread_index + 1) * chunk_size;
//             for (size_t i = start; i < end; ++i) 
//             {
//                 size_t module_i = i / module_size;
//                 Num sum = 0.0;
//                 for (size_t j = 0; j < N; ++j) 
//                 {
//                     size_t module_j = j / module_size;
//                     if ( i != j )
//                     {
//                         if ( module_i == module_j )
//                         {
//                             sum += intra_k * std::sin(theta[j] - theta[i] - alpha);
//                         }
//                         else
//                         {
//                             sum += inter_k * std::sin(theta[j] - theta[i] - alpha);
//                         }
//                     }
//                 }
//                 dthetadt[i] = omega[i] + sum;
//             }
//         });
//     }
//     for (auto& t : threads) t.join();
//     // return dtheta_dt;
// }

// ======================================= //
//                                         //
// ---------------- Wrapper -------------- //
//                                         //
// ======================================= //
template <FPNumber Num>
struct KuramotoModularParams
{
    // Model specific parameters in modular Kuramoto model (time and phases/theta are dynamics bounds)
    Vec<Num> omega;       // Natural frequencies
    Num      intra_K;     // Intra-module coupling strength
    Num      inter_K;     // Inter-module coupling strength
    Num      alpha;
    size_t   module_size; // Size of each module
    size_t   num_modules;
    int      N;           // Number of oscillators

    // KuramotoModularParams(size_t mod_size=50, size_t num_mods=2, Num phae_lag = 0.0)
    //     : N(mod_size * num_mods),
    //       omega(Vec<Num>(mod_size * num_mods, 0.0)),
    //       module_size(mod_size),
    //       alpha(phae_lag)
    // {}
};
template <FPNumber Num>
inline MyFunc<Num> auto kuramoto_special_modular_wrapper(const KuramotoModularParams<Num>& params)
{
    return [params](Num t, const Vec<Num>& theta, Vec<Num>& dthetadt) -> void
    {
        return kuramoto_special_modular(t, theta, dthetadt, params.omega, params.intra_K, params.inter_K, params.alpha, params.module_size);
    };
}
} // End MathEngine namespace

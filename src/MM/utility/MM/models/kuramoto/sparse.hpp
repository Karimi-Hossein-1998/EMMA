#pragma once
#include "../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
// General Kuramoto model with phase-lag and sparse dMatrix
template <FPNumber Num>
inline void kuramoto_sparse(
    Num                       time,
    const Vec<Num>&           theta,
    Vec<Num>&                 dthetadt,
    const Vec<Num>&           omega,
    Num                       K,
    const SparsedMatrix<Num>& sparse_adj,
    Num                       alpha
)
{
    const size_t N = theta.size(); // N > 0 here

    // Original logic
    Num k_norm    = K / static_cast<Num>( N ); // Renamed k to k_norm
    for ( size_t i = 0; i < N; ++i )
    {
        Num sum = 0.0;
        for ( const auto& edge_pair : sparse_adj.rows[i] ) // Renamed [j, weight] to edge_pair for clarity before C++17 structured binding decomposition
        {
            size_t j = edge_pair.first;
            Num weight = edge_pair.second;
            sum += weight * std::sin( theta[j] - theta[i] - alpha );
        }
        dthetadt[i] = omega[i] + k_norm * sum;
    }
    // return dtheta_dt;
}

// // General Kuramoto model with phase-lag and sparse dMatrix (parallel)
// inline void kuramoto_sparse_parallel(
//     Num              time,
//     const Vec<Num>&         theta,
//     Vec<Num>&               dthetadt,
//     const Vec<Num>&         omega,
//     Num              K,
//     const SparsedMatrix& sparse_adj,
//     Num              alpha
// )
// {
//     const size_t N = theta.size(); // N > 0 here

//     // Original logic
//     Num k_norm    = K / static_cast<Num>( N ); // Renamed k to k_norm
//     std::vector<std::thread> threads;
//     size_t num_threads = std::min(N, static_cast<size_t>(std::max(1u, std::thread::hardware_concurrency())));
//     if ( num_threads == 0 ) num_threads = 1;
//     size_t chunk_size  = N / num_threads;
    
//     for (size_t t = 0; t < num_threads; ++t) 
//     {
//         threads.emplace_back([&, t]() 
//         {
//             size_t start = t * chunk_size;
//             size_t end = (t == num_threads - 1) ? N : (t + 1) * chunk_size;
            
//             for (size_t i = start; i < end; ++i) 
//             {
//                 Num sum = 0.0;
//                 for ( const auto& edge_pair : sparse_adj.rows[i] ) // Renamed [j, weight] to edge_pair for clarity before C++17 structured binding decomposition
//                 {
//                     size_t j      = edge_pair.first;
//                     Num weight = edge_pair.second;
//                     sum          += weight * std::sin( theta[j] - theta[i] - alpha );
//                 }
//                 dthetadt[i] = omega[i] + k_norm * sum;
//             }
//         });
//     }
//     for (auto& t : threads) t.join();
//     // return dtheta_dt;
// }

// ======================================== //
//                                          //
// ------------ Wrappers ------------------ //
//                                          //
// ======================================== //
template <FPNumber Num>
struct KuramotoSparseParams
{
    // Model specific parameters in Kuramoto model (time and phases/theta are dynamics bounds)
    SparsedMatrix<Num> sparse_adj;
    Vec<Num>           omega; // Natural frequencies
    Num                alpha;
    Num                K;     // Coupling strength
    int                N;     // Number of oscillators

    // KuramotoSparseParams(int num_oscillators=50, Num phae_lag = 0.0)
    //     : N(num_oscillators),
    //       omega(Vec<Num>(num_oscillators, 0.0)),
    //       sparse_adj(SparsedMatrix(num_oscillators)),
    //       alpha(phae_lag)
    // {}
};

template <FPNumber Num>
inline MyFunc<Num> auto kuramoto_sparse_wrapper(const KuramotoSparseParams<Num>& params)
{
    return [params](Num time, const Vec<Num>& theta, Vec<Num>& dthetadt) -> void
    {
        return kuramoto_sparse(time, theta, dthetadt, params.omega, params.K, params.sparse_adj, params.alpha);
    };
}

// inline MyFunc kuramoto_sparse_parallel_wrapper(const KuramotoSparseParams& params)
// {
//     return [params](Num time, const Vec<Num>& theta, Vec<Num>& dthetadt) -> void
//     {
//         return kuramoto_sparse_parallel(time, theta, dthetadt, params.omega, params.K, params.sparse_adj, params.alpha);
//     };
// }
} // End namespace MathEngine

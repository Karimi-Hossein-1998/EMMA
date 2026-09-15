#pragma once
#include "../../typedefs/header.hpp"

namespace MathEngine
{ // MathEngine namespace
// General Kuramoto model with phase-lag and sparse dMatrix
inline void kuramoto_sparse(
    double               time,
    const Vec<double>&   theta,
    Vec<double>&         dthetadt,
    const Vec<double>&   omega,
    double               K,
    const SparsedMatrix& sparse_adj,
    double               alpha
)
{
    const size_t N = theta.size(); // N > 0 here

    // Original logic
    double k_norm    = K / static_cast<double>( N ); // Renamed k to k_norm
    for ( size_t i = 0; i < N; ++i )
    {
        double sum = 0.0;
        for ( const auto& edge_pair : sparse_adj.rows[i] ) // Renamed [j, weight] to edge_pair for clarity before C++17 structured binding decomposition
        {
            size_t j = edge_pair.first;
            double weight = edge_pair.second;
            sum += weight * std::sin( theta[j] - theta[i] - alpha );
        }
        dthetadt[i] = omega[i] + k_norm * sum;
    }
    // return dtheta_dt;
}

// // General Kuramoto model with phase-lag and sparse dMatrix (parallel)
// inline void kuramoto_sparse_parallel(
//     double              time,
//     const Vec<double>&         theta,
//     Vec<double>&               dthetadt,
//     const Vec<double>&         omega,
//     double              K,
//     const SparsedMatrix& sparse_adj,
//     double              alpha
// )
// {
//     const size_t N = theta.size(); // N > 0 here

//     // Original logic
//     double k_norm    = K / static_cast<double>( N ); // Renamed k to k_norm
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
//                 double sum = 0.0;
//                 for ( const auto& edge_pair : sparse_adj.rows[i] ) // Renamed [j, weight] to edge_pair for clarity before C++17 structured binding decomposition
//                 {
//                     size_t j      = edge_pair.first;
//                     double weight = edge_pair.second;
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
struct KuramotoSparseParams
{
    // Model specific parameters in Kuramoto model (time and phases/theta are dynamics bounds)
    SparsedMatrix sparse_adj = SparsedMatrix();
    Vec<double>   omega; // Natural frequencies
    double        alpha;
    double        K;     // Coupling strength
    int           N;     // doubleber of oscillators
};

inline MyFunc kuramoto_sparse_wrapper(const KuramotoSparseParams& params)
{
    return [params](double time, const Vec<double>& theta, Vec<double>& dthetadt) -> void
    {
        return kuramoto_sparse(time, theta, dthetadt, params.omega, params.K, params.sparse_adj, params.alpha);
    };
}

// inline MyFunc kuramoto_sparse_parallel_wrapper(const KuramotoSparseParams& params)
// {
//     return [params](double time, const Vec<double>& theta, Vec<double>& dthetadt) -> void
//     {
//         return kuramoto_sparse_parallel(time, theta, dthetadt, params.omega, params.K, params.sparse_adj, params.alpha);
//     };
// }
} // End namespace MathEngine

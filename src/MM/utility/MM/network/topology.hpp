#pragma once
#include "../typedefs/header.hpp"
#include <cstddef>

namespace MathEngine
{ // namespace MathEngine
// -----------------------------------------------------------------------------
// Convert dense Matrix<Num> to sparse Matrix<Num>
/*----------------------------------------------------------*/
// Parameters:
//   adj: Dense adjacency matrix (N x N)
// Returns:
//   SparsedMatrix<Num> representation of the input matrix
template <FPNumber Num>
inline SparsedMatrix<Num> dense_to_sparse(const Matrix<Num>& adj)
{
    size_t N = adj.Rows();
    SparsedMatrix<Num> sparse(N);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            if (adj[i,j] != 0.0)
                sparse.rows[i].emplace_back(j, adj[i,j]);
    return sparse;
}

// -----------------------------------------------------------------------------
// Compute density (fraction of nonzero off-diagonal entries) for a dense Matrix<Num>
/*----------------------------------------------------------*/
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider an entry as nonzero (default 1e-12)
// Returns:
//   Density (fraction of nonzero off-diagonal entries)
template <FPNumber Num>
inline Num density(const Matrix<Num>& adj, std::float64_t threshold = 1e-12)
{
    size_t N = adj.size(), nonzero = 0, total = N * (N - 1);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            if (i != j && std::abs(adj[i,j]) > threshold)
                ++nonzero;
    return total ? static_cast<Num>(nonzero*1.0) / static_cast<Num>(total*1.0) : 0.0;
}

// -----------------------------------------------------------------------------
// Convert dense Matrix<Num> to sparse Matrix<Num> if sparse enough, else return density
/*----------------------------------------------------------*/
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   sparse_adj: Output sparse matrix (only filled if density <= density_threshold)
//   density_threshold: Maximum density to allow conversion to sparse (default 0.5)
//   zero_threshold: Minimum absolute value to consider as nonzero (default 1e-12)
// Returns:
//   Pair: (true if sparse_adj is filled, false if only density is set; density value)
template <FPNumber Num>
inline std::pair<bool, Num> dense_to_sparse_conditional(
    const Matrix<Num> &adj,
    SparsedMatrix<Num> &sparse_adj,
    Num density_threshold = static_cast<Num>(0.5),
    Num zero_threshold    = static_cast<Num>(1e-12)
)
{
    size_t N = adj.size();
    auto density_measure = density(adj, zero_threshold);
    if (density_measure <= density_threshold)
    {
        sparse_adj = SparsedMatrix<Num>(N);
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                if (std::abs(adj[i,j]) > zero_threshold)
                {
                    sparse_adj.rows[i].emplace_back(j, adj[i,j]);
                }
            }
        }
        return {true, density_measure};
    }
    else
    {
        return {false, density_measure};
    }
}

// Random Matrix<Num>: each edge weight is random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated dense adjacency matrix (N x N)
template <FPNumber Num>
inline Matrix<Num> random(
    size_t   N          = 100,
    Num      min_weight = static_cast<Num>(0.0),
    Num      max_weight = static_cast<Num>(1.0),
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (min_weight > max_weight) std::swap(min_weight,max_weight);
    // Original logic
    std::mt19937 rng(seed);
    // uniform_real_distribution handles min_weight == max_weight correctly.
    std::uniform_real_distribution<Num> dist(min_weight, max_weight);
    Matrix<Num> adj(N,N);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            adj[i,j] = (i == j) ? 0.0 : dist(rng);
        }
    }
    return adj;
}

// Random Symmetric Matrix<Num>: each edge weight is random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated dense adjacency matrix (N x N)
template <FPNumber Num>
inline Matrix<Num> random_symmetric(
    size_t   N          = 100,
    Num      min_weight = static_cast<Num>(0.0),
    Num      max_weight = static_cast<Num>(1.0),
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (min_weight > max_weight) std::swap(min_weight,max_weight);

    // Original logic
    std::mt19937 rng(seed);
    // uniform_real_distribution handles min_weight == max_weight correctly.
    std::uniform_real_distribution<double> dist(min_weight, max_weight);
    Matrix<Num> adj(N,N);
    for (size_t i = 1; i < N; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            adj[i,j] = dist(rng);
            adj[j,i] = adj[i,j];
        }
    }
    return adj;
}

// Erdos-Renyi Matrix<Num>: each edge exists with probability p, weight random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi adjacency matrix (N x N)
template <FPNumber Num>
inline Matrix<Num> erdos_renyi(
    size_t   N          = 100,
    Num      p          = static_cast<Num>(0.5),
    Num      min_weight = static_cast<Num>(0.0),
    Num      max_weight = static_cast<Num>(1.0),
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (p < 0.0) return Matrix<Num>(N,dVec(N,0.0));
    if (p >= 1.0) return random(N,min_weight,max_weight,seed);
    if (min_weight > max_weight) std::swap(min_weight,max_weight);

    // Original logic
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Num> weight_dist(min_weight, max_weight); // Renamed dist
    std::bernoulli_distribution edge_dist(p);
    Matrix<Num> adj(N,N);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (i != j && edge_dist(rng))
            {
                adj[i,j] = weight_dist(rng);
            }
        }
    }
    return adj;
}

// Erdos-Renyi Matrix<Num> (Uniform): exactly N*(N-1)*p edges with random weights
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi adjacency matrix (N x N) with exactly N*(N-1)*p edges
template <FPNumber Num>
inline Matrix<Num> erdos_renyi_uniform(
    size_t   N          = 100,
    Num      p          = static_cast<Num>(0.5),
    Num      min_weight = static_cast<Num>(0.0),
    Num      max_weight = static_cast<Num>(1.0),
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (p < 0.0) return Matrix<Num>(N,dVec(N,0.0));
    if (p >= 1.0) return random(N,min_weight,max_weight,seed);
    if (min_weight > max_weight) std::swap(min_weight,max_weight);

    std::mt19937 rng(seed);
    std::uniform_real_distribution<Num> weight_dist(min_weight, max_weight);
    Matrix<Num> adj(N,N);
    // Calculate exact number of edges needed
    size_t total_edges = static_cast<size_t>(N * (N - 1) * p);
    // Create vector of all possible edges
    wpairVec possible_edges;
    possible_edges.reserve(N * (N - 1));  // Reserve space for all possible edges
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            if (i != j)  // No self-loops
                possible_edges.emplace_back(i, j);

    // Shuffle the edges
    std::shuffle(possible_edges.begin(), possible_edges.end(), rng);
    // Create exactly total_edges edges
    for (size_t e = 0; e < total_edges; ++e)
    {
        size_t i  = possible_edges[e].first;
        size_t j  = possible_edges[e].second;
        adj[i,j] = weight_dist(rng);
    }

    return adj;
}

// Symmetric Erdos-Renyi Matrix<Num>: each edge exists with probability p, weight random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi (Symmetric) adjacency matrix (N x N)
template <FPNumber Num>
inline Matrix<Num> erdos_renyi_symmetric(
    size_t   N          = 100,
    Num      p          = static_cast<Num>(0.5),
    Num      min_weight = static_cast<Num>(0.0),
    Num      max_weight = static_cast<Num>(1.0),
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (p < 0.0) return Matrix<Num>(N,dVec(N,0.0));
    if (p >= 1.0) return random_symmetric(N,min_weight,max_weight,seed);
    if (min_weight > max_weight) std::swap(min_weight,max_weight);

    // Original logic
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> weight_dist(min_weight, max_weight); // Renamed dist
    std::bernoulli_distribution edge_dist(p);
    Matrix<Num> adj(N,N);
    for (size_t i = 1; i < N; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            if (edge_dist(rng))
            {
                adj[i,j] = weight_dist(rng);
                adj[j,i] = adj[i,j];
            }
        }
    }
    return adj;
}

// Symmetric Erdos-Renyi Matrix<Num> (Uniform): exactly N*(N-1)*p/2 edges with random weights
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi (Symmetric) adjacency matrix (N x N) with exactly N*(N-1)*p/2 edges
template <FPNumber Num>
inline Matrix<Num> erdos_renyi_symmetric_uniform(
    size_t   N          = 100,
    Num      p          = static_cast<Num>(0.5),
    Num      min_weight = static_cast<Num>(0.0),
    Num      max_weight = static_cast<Num>(1.0),
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (p < 0.0) return Matrix<Num>(N,dVec(N,0.0));
    if (p >= 1.0) return random_symmetric(N,min_weight,max_weight,seed);
    if (min_weight > max_weight) std::swap(min_weight,max_weight);

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> weight_dist(min_weight, max_weight);
    Matrix<Num> adj(N,N);

    // Calculate exact number of edges needed
    size_t total_edges = static_cast<size_t>( ( N * (N - 1) * 0.5 ) * p );
    
    // Create vector of all possible edges (only lower triangle)
    wpairVec possible_edges;
    possible_edges.reserve( N * (N - 1) * 0.5 );  // Reserve space for all possible edges
    for (size_t i = 1; i < N; ++i)
        for (size_t j = 0; j < i; ++j)
            possible_edges.emplace_back(i, j);

    // Shuffle the edges
    std::shuffle(possible_edges.begin(), possible_edges.end(), rng);
    // Create exactly total_edges edges
    for (size_t e = 0; e < total_edges; ++e)
    {
        size_t i      = possible_edges[e].first;
        size_t j      = possible_edges[e].second;
        Num    weight = weight_dist(rng);
        adj[i,j]     = weight;
        adj[j,i]     = weight;  // Make it symmetric
    }

    return adj;
}

// Small-world Matrix<Num> (Watts-Strogatz model, ring lattice with rewiring)
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   k: Each node is connected to k nearest neighbors in ring topology (default 4)
//   beta: Rewiring probability (default 0.5)
//   weight: Edge weight for all connections (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Small-world adjacency matrix (N x N)
template <FPNumber Num>
inline Matrix<Num> small_world(
    size_t N      = 100,
    size_t k      = 4,
    Num    beta   = static_cast<Num>(0.5),
    Num    weight = static_cast<Num>(1.0),
    unsigned seed = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (k >= N/2) return random(N,weight,weight,seed);
    if (k == 0) return Matrix<Num>(N,dVec(N,0.0));
    if (beta < 0.0) beta = 0.0;
    else if (beta >1.0) beta = 1.0;
    std::mt19937 rng(seed);
    Matrix<Num> adj(N,N);

    // Initial ring lattice - ensure k/2 neighbors on each side
    size_t half_k = k / 2; // Integer division intentional
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= half_k; ++j)
        {
            size_t right  = (i + j) % N;
            size_t left   = (i + N - j) % N;
            adj[i,right] = weight;
            adj[i,left]  = weight;
        }
        // Handle odd k by adding one more connection
        if (k % 2 == 1 && N > 2*(k/2))
        {
            size_t extra  = (i + (k / 2 + 1)) % N;
            adj[i,extra] = weight;
            adj[extra,i] = weight;
        }
    }
    // Rewiring with improved logic
    std::bernoulli_distribution rewire_dist(beta);
    std::uniform_int_distribution<size_t> node_dist(0, N - 1);

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= half_k; ++j)
        {
            size_t neighbor = (i + j) % N;
            if (rewire_dist(rng))
            {
                size_t attempts     = 0;
                size_t max_attempts = N * static_cast<size_t>(std::sqrt(N)); // Reasonable limit to prevent infinite loops
                bool   found        = false;
                size_t new_neighbor;

                while (!found && attempts < max_attempts)
                {
                    new_neighbor = node_dist(rng);
                    if (new_neighbor        != i &&
                        new_neighbor        != neighbor &&
                        adj[i,new_neighbor] == 0.0)
                    {
                        found = true;
                    }
                    attempts++;
                }

                if (found)
                {
                    adj[i,neighbor]     = 0.0;
                    adj[neighbor,i]     = 0.0;
                    adj[i,new_neighbor] = weight;
                    adj[new_neighbor,i] = weight;
                }
            }
        }

        // Handle rewiring for the extra connection in case of odd k
        if (k % 2 == 1 && N > 2)
        {
            size_t extra = (i + (k / 2 + 1)) % N;
            if (rewire_dist(rng))
            {
                size_t attempts     = 0;
                size_t max_attempts = N * static_cast<size_t>(std::sqrt(N));
                bool   found        = false;
                size_t new_neighbor;

                while (!found && attempts < max_attempts)
                {
                    new_neighbor              = node_dist(rng);
                    if (new_neighbor         != i &&
                        new_neighbor         != extra &&
                        adj[i,new_neighbor] == 0.0)
                    {
                        found = true;
                    }
                    attempts++;
                }

                if (found)
                {
                    adj[i,extra]        = 0.0;
                    adj[extra,i]        = 0.0;
                    adj[i,new_neighbor] = weight;
                    adj[new_neighbor,i] = weight;
                }
            }
        }
    }
    return adj;
}

// Directed Small-world Matrix<Num> (Watts-Strogatz model, ring lattice with rewiring)
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   k: Each node has k outgoing edges in ring topology (default 4)
//   beta: Rewiring probability (default 0.5)
//   weight: Edge weight for all connections (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Directed small-world adjacency matrix (N x N)
template <FPNumber Num>
inline Matrix<Num> small_world_directed(
    size_t N      = 100,
    size_t k      = 4,
    Num    beta   = static_cast<Num>(0.5),
    Num    weight = static_cast<Num>(1.0),
    unsigned seed = 42
)
{
    // Error Handling
    if (N == 0) return {};
    if (k >= N/2) return random(N,weight,weight,seed);
    if (k == 0) return Matrix<Num>(N,dVec(N,0.0));
    if (beta < 0.0) beta = 0.0;
    else if (beta >1.0) beta = 1.0;

    std::mt19937 rng(seed);
    Matrix<Num> adj(N,N);

    // Initial ring lattice - create k outgoing edges for each node
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= k; ++j)
        {
            size_t target  = (i + j) % N;
            adj[target,i] = weight;
        }
    }

    // Rewiring with improved logic
    std::bernoulli_distribution           rewire_dist(beta);
    std::uniform_int_distribution<size_t> node_dist(0, N - 1);

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= k; ++j)
        {
            size_t neighbor = (i + j) % N;
            if (rewire_dist(rng))
            {
                size_t attempts     = 0;
                size_t max_attempts = N * static_cast<size_t>(std::sqrt(N)); // Reasonable limit to prevent infinite loops
                bool   found        = false;
                size_t new_neighbor;

                while (!found && attempts < max_attempts)
                {
                    new_neighbor = node_dist(rng);
                    if (new_neighbor != i && 
                        new_neighbor != neighbor && 
                        adj[i,new_neighbor] == 0.0)
                    {
                        found = true;
                    }
                    attempts++;
                }

                if (found)
                {
                    adj[neighbor,i]     = 0.0;
                    adj[new_neighbor,i] = weight;
                }
            }
        }
    }
    return adj;
}

// Modular Matrix<Num>: nodes are divided into modules, dense within, sparse between
/*----------------------------------------------------------*/
// Parameters:
//   module_size: Number of nodes per module (default 100)
//   num_modules: Number of modules (default 10)
//   p_in: Probability of within-module connection (default 0.9)
//   p_out: Probability of between-module connection (default 0.1)
//   in_weight: Weight for within-module connections (default 1.0)
//   out_weight: Weight for between-module connections (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Modular adjacency matrix (N x N), where N = module_size * num_modules
template <FPNumber Num>
inline Matrix<Num> modular(
    size_t   module_size = 100,
    size_t   num_modules = 10,
    Num      p_in        = static_cast<Num>(0.9),
    Num      p_out       = static_cast<Num>(0.1),
    Num      in_weight   = static_cast<Num>(1.0),
    Num      out_weight  = static_cast<Num>(1.0),
    unsigned seed        = 42
)
{
    // Error handling
    if (module_size == 0) return {};
    if (num_modules == 0) return {};
    if (p_in < 0.0) p_in = 0.0;
    else if (p_in > 1.0) p_in = 1.0;
    if (p_out < 0.0) p_out = 0.0;
    else if (p_out > 1.0) p_out = 1.0;
    // Calculate total network size
    size_t N = module_size*num_modules;
	Matrix<Num> adj(N,N);
    std::mt19937 rng(seed);
    std::bernoulli_distribution in_dist(p_in), out_dist(p_out);
    for (size_t m = 0; m < num_modules; ++m)
    {
        size_t start = m * module_size;
        size_t end   = (m == num_modules - 1) ? N : (m + 1) * module_size;
        // Within module
        for (size_t i = start; i < end; ++i)
        {
            for (size_t j = start; j < end; ++j)
            {
                if (i != j && in_dist(rng))
                {
                    adj[i,j] = in_weight;
                }
            }
        }
        // Between modules
        for (size_t n = m + 1; n < num_modules; ++n)
        {
            size_t n_start = n * module_size;
            size_t n_end   = (n == num_modules - 1) ? N : (n + 1) * module_size;
            for (size_t i = start; i < end; ++i)
            {
                for (size_t j = n_start; j < n_end; ++j)
                {
                    if (out_dist(rng))
                    {
                        adj[i,j] = out_weight;
                        adj[j,i] = out_weight;
                    }
                }
            }
        }
    }
    return adj;
}

// Hierarchical Matrix<Num>: recursively nested modules with level-dependent connection probabilities
/*----------------------------------------------------------*/
// Parameters:
//   N: Base module size at the lowest level (default 100)
//   levels: Number of hierarchical levels (0 returns Erdos-Renyi N x N matrix) (default 2)
//   p_in: Base within-module connection probability (decays with level) (default 0.9)
//   p_out: Base between-module connection probability (decays with level) (default 0.1)
//   in_weight: Weight for within-module connections (default 1.0)
//   out_weight: Weight for between-module connections (default 1.0)
//   level_decay: Decay factor for connection probabilities at each level (0 < level_decay <= 1) (default 0.5)
//   seed: Random seed (default 42)
//   base_module_num: Number of modules in the base level of the hierarchy (default 2)
// Returns:
//   Hierarchical adjacency matrix (size: N * 2^(levels-1) * base_module_num)
template <FPNumber Num>
inline Matrix<Num> hierarchical(
    size_t   N               = 100,
    size_t   levels          = 2,
    Num      p_in            = static_cast<Num>(0.9),
    Num      p_out           = static_cast<Num>(0.1),
    Num      in_weight       = static_cast<Num>(1.0),
    Num      out_weight      = static_cast<Num>(1.0),
    Num      level_decay     = static_cast<Num>(0.5),
    unsigned seed            = 42,
    size_t   base_module_num = 2
)
{
    // Parameter validation and base cases
    if (N == 0) return {};
    if (base_module_num == 0) return {};
    p_in = std::clamp(p_in, 0.0, 1.0);
    p_out = std::clamp(p_out, 0.0, 1.0);
    level_decay = std::clamp(level_decay,0.0,1.0);
    if (levels == 0) return erdos_renyi(N, p_in, in_weight, out_weight, seed);
    // Check for potential size overflow
    size_t level_multiplier=(static_cast<size_t>(1)<<(levels-1)),module_num = base_module_num*level_multiplier;
    size_t max_nodes = N * module_num;
    if (max_nodes / (N * base_module_num) != level_multiplier) return erdos_renyi(module_num,p_in,in_weight,out_weight,seed);
    Num this_p_in  = p_in * std::pow(level_decay, levels - 1);
    Num this_p_out = p_out * std::pow(level_decay, levels - 1);
    if (levels == 1)
    {
        return modular(N, base_module_num, this_p_in, this_p_out, in_weight, out_weight, seed);
    }
    Matrix<Num> adj(max_nodes,max_nodes);
    // Recursively build left and right submodules
    auto left  = hierarchical(N, levels - 1, p_in, p_out, in_weight, out_weight, level_decay, seed + 1, base_module_num);
    auto right = hierarchical(N, levels - 1, p_in, p_out, in_weight, out_weight, level_decay, seed + level_multiplier, base_module_num);
    size_t half  = N * (static_cast<size_t>(1)<<(levels-2)) * base_module_num;
    for (size_t i = 0; i < half; ++i)
    {
        for (size_t j = 0; j < half; ++j)
        {
            adj[i,j]               = left[i,j];
            adj[half + i,half + j] = right[i,j];
        }
    }
    std::mt19937                rng(seed);
    std::bernoulli_distribution out_dist(this_p_out);
    for (size_t i = 0; i < half; ++i)
    {
        for (size_t j = half; j < max_nodes; ++j)
        {
            if (out_dist(rng))
            {
                adj[i,j] = out_weight;
                adj[j,i] = out_weight;
            }
        }
    }
    return adj;
}

enum class NetworkTopology
{
    Uniform=0,
    UniformSymmetric,
    ErdosRenyi,
    ErdosRenyiUniform,
    ErdosRenyiSymmetric,
    ErdosRenyiSymmetricUniform,
    SmallWorld,
    SmallWorldDirected,
    Modular,
    Hierarchical
};

template <FPNumber Num>
struct NetworkParams
{
    NetworkTopology topology=NetworkTopology::Uniform;
    Num    p1      = static_cast<Num>(1.0);
    Num    p2      = static_cast<Num>(0.5);
    Num    weight1 = static_cast<Num>(1.0);
    Num    weight2 = static_cast<Num>(0.5);
    Num    ld      = static_cast<Num>(0.1);
    size_t degree  = 2;
    size_t seed    = 41;
    size_t baseN   = 50;
    size_t N       = 50;
    size_t sModule = 50;
    size_t nModule = 1;
    size_t hLevel  = 1;
};

template <FPNumber Num>
inline Matrix<Num> initialize_network(const NetworkParams<Num>& networkParams)
{
	switch(networkParams.topology)
    {
        case NetworkTopology::Uniform: return random(networkParams.N,networkParams.weight1,networkParams.weight2,networkParams.seed);break;
        case NetworkTopology::UniformSymmetric:
            return random_symmetric(networkParams.N,networkParams.weight1,networkParams.weight2,networkParams.seed);break;
        case NetworkTopology::ErdosRenyi:
            return erdos_renyi(networkParams.N,networkParams.p1,networkParams.weight1,networkParams.weight2,networkParams.seed);break;
        case NetworkTopology::ErdosRenyiUniform:
            return erdos_renyi_uniform(networkParams.N,networkParams.p1,networkParams.weight1,networkParams.weight2,networkParams.seed);break;
        case NetworkTopology::ErdosRenyiSymmetric:
            return erdos_renyi_symmetric(networkParams.N,networkParams.p1,networkParams.weight1,networkParams.weight2,networkParams.seed);break;
        case NetworkTopology::ErdosRenyiSymmetricUniform:
            return erdos_renyi_symmetric_uniform(networkParams.N,networkParams.p1,networkParams.weight1,networkParams.weight2,networkParams.seed);
            break;
        case NetworkTopology::SmallWorld:
            return small_world(networkParams.N,networkParams.degree,networkParams.p1,networkParams.weight1,networkParams.seed);break;
        case NetworkTopology::SmallWorldDirected:
            return small_world_directed(networkParams.N,networkParams.degree,networkParams.p1,networkParams.weight1,networkParams.seed);break;
        case NetworkTopology::Modular:
            return modular(networkParams.sModule,networkParams.nModule,networkParams.p1,networkParams.p2,
                           networkParams.weight1,networkParams.weight2,networkParams.seed);break;
        case NetworkTopology::Hierarchical:
            return hierarchical(networkParams.sModule,networkParams.hLevel,networkParams.p1,networkParams.p2,
                                networkParams.weight1,networkParams.weight2,networkParams.ld,networkParams.seed,networkParams.nModule);break;
        default: return random_symmetric(networkParams.N,networkParams.weight1,networkParams.weight2,networkParams.seed);break;
    }
}

// -----------------------------------------------------------------------------
// Sparse small-world Matrix<Num> (Watts-Strogatz, returns SparsedMatrix<Num>)
template <FPNumber Num>
inline SparsedMatrix<Num> small_world_sparse(
    size_t   N,
    size_t   k,
    Num      beta,
    Num      weight,
    unsigned seed = 42
)
{
    auto   dense  = small_world(N, k, beta, weight, seed);
    auto   sparse = dense_to_sparse(dense);
    return sparse;
}


// Multilayered Matrix<Num>: block diagonal matrix, each block is a layer
template <FPNumber Num>
inline Matrix<Num> multilayered(const Vec<Matrix<Num>> &layers)
{
    size_t offset = 0;
    size_t total_N = 0;
    for (const auto &layer : layers) total_N += layer.size();
    Matrix<Num> adj(total_N, dVec(total_N, 0.0));
    for (const auto &layer : layers)
    {
        size_t n = layer.size();
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                adj[offset+i,offset + j] = layer[i][j];
            }
        }
        offset += n;
    }
    return adj;
}

// effective_multiplex: Weighted sum of multiple network layers
/*----------------------------------------------------------*/
// Parameters:
//   layers: Vector of adjacency matrices (all must be N x N)
//   layer_weights: Vector of weights (one per layer)
// Returns:
//   Weighted sum of layers as a single adjacency matrix
template <FPNumber Num>
inline Matrix<Num> effective_multiplex(
    const Vec<Matrix<Num>> &layers,
    const Vec<Num>         &layer_weights
)
{
    if (layers.empty()) return {};
    size_t num_layers = layers.size();
    if (num_layers > layer_weights.size()) {};
    const size_t N = layers[0].size();
    if (N == 0) return {};
    for (size_t l = 0; l < num_layers; ++l)
    {
        if (layers[l].size() != N) return {};
        for (size_t i = 0; i < N; ++i)
        {
            if (layers[l][i].size() != N) return {};
        }
    }
    Matrix<Num> adj_out(N, dVec(N, 0.0));
    for (size_t l = 0; l < num_layers; ++l)
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
                adj_out[i][j] += layer_weights[l] * layers[l][i][j];
    return adj_out;
}

// generate_multiplex_network_layers: Generate multiple network layers for a multiplex network
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (consistent across all layers) (default 100)
//   num_layers: Number of layers to generate (default 2)
//   layer_generation_types: Vector of topology types for each layer (default "erdos_renyi")
//   layer_params: Matrix<Num> of parameters for each layer (one row per layer)
//   layer_seeds: Vector of random seeds for each layer (default 42)
// Returns:
//   Vector of adjacency matrices, one per layer
template <FPNumber Num>
inline Vec<Matrix<Num>> generate_multiplex_network_layers(const Vec<NetworkParams<Num>>& layer_configs)
{
    if (layer_configs.empty()) return {};
    Vec<Matrix<Num>> layers; layers.reserve(layer_configs.size());
    for (const auto& config : layer_configs) layers.push_back(initialize_network(config));
	return layers;
}

// -----------------------------------------------------------------------------
// Calculate in-degrees of a given adjacency matrix
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of in-degrees (number of incoming edges for each node)
template <FPNumber Num>
inline wVec in_degrees(const Matrix<Num>& adj, Num threshold = 1e-12)
{
    size_t N = adj.size();
    wVec   indeg(N, 0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            if (std::abs(adj[i,j]) > threshold)
                ++indeg[i];
    return indeg;
}

// -----------------------------------------------------------------------------
// Calculate out-degrees of a given adjacency matrix
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of out-degrees (number of outgoing edges for each node)
template <FPNumber Num>
inline wVec out_degrees(const Matrix<Num>& adj, Num threshold = 1e-12)
{
    size_t N = adj.size();
    wVec   outdeg(N, 0);
    for (size_t j = 0; j < N; ++j)
        for (size_t i = 0; i < N; ++i)
            if (std::abs(adj[i,j]) > threshold)
                ++outdeg[j];
    return outdeg;
}

// -----------------------------------------------------------------------------
// Calculate in-degrees of a given sparse adjacency matrix
// Parameters:
//   sparse: Sparse adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of in-degrees (number of incoming edges for each node)
template <FPNumber Num>
inline wVec in_degrees(const SparsedMatrix<Num>& sparse, Num threshold = 1e-12)
{
    size_t N = sparse.rows.size();
    wVec   indeg(N, 0);
    for (size_t i = 0; i < N; ++i)
    {
        for (const auto& entry : sparse.rows[i])
        {
            size_t j     = entry.first;
            double value = entry.second;
            if (std::abs(value) > threshold)
                ++indeg[i];
        }
    }
    return indeg;
}

// -----------------------------------------------------------------------------
// Calculate out-degrees of a given sparse adjacency matrix
// Parameters:
//   sparse: Sparse adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of out-degrees (number of outgoing edges for each node)
template <FPNumber Num>
inline wVec out_degrees(const SparsedMatrix<Num>& sparse, Num threshold = 1e-12)
{
    size_t N = sparse.rows.size();
    wVec   outdeg(N, 0);
    for (size_t i = 0; i < N; ++i)
    {
        for (const auto& entry : sparse.rows[i])
        {
            size_t j     = entry.first;
            double value = entry.second;
            if (std::abs(value) > threshold)
                ++outdeg[j];
        }
    }
    return outdeg;
}
} // End MathEngine namespace

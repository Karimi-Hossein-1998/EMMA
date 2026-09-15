#pragma once
#include "../../typedefs/header.hpp"
#include <ranges>
namespace MathEngine
{ // namespace MathEngine
inline dVec calculate_order(const dVec& phases)
{
    const auto N = phases.size();
    dVec  once_order(3, 0.0);
    for ( auto theta : phases )
    {
        once_order[0] += std::sin(theta);
        once_order[1] += std::cos(theta);
    }
    once_order[0] /= N;
    once_order[1] /= N;
    once_order[2] = std::sqrt( once_order[0] * once_order[0] + once_order[1] * once_order[1] );
    return once_order;
}

inline dMatrix calculate_order(const SolverResults& results)
{
    const auto&  sol      = results.solution;
    const auto&  time     = results.timePoints;
    const auto   N        = sol[0].size();
    const double NinVerse = 1.0 / static_cast<double>(N);
    dMatrix order;

    if (time.size() != sol.size())
    {
        throw std::invalid_argument("[calculate_order] Dimensions mismatch: timePoints.size() (" +
              std::to_string(time.size()) + ") != solution.size() (" + 
              std::to_string(sol.size()) + ")");
    }

    dVec  tempOrder(3, 0.0);
    for (auto [sol_iter, time_iter] : std::views::zip(sol, time))
    {
        tempOrder = calculate_order(sol_iter);
        tempOrder.push_back(time_iter);
        order.push_back(tempOrder);
    }
    
    return order;
}

inline dVec  calculate_order_per_node( 
    const  dVec& phases, 
    size_t node_number, 
    const  dMatrix& adj,
    double tolerance = 1e-12
)
{
    auto   N       = phases.size();
    double weight  = 0.0;
    size_t counter = 0;
    dVec   local_order(3,0.0);
    for ( size_t i = 0; i<N; ++i )
    {
        if ( std::abs(adj[node_number][i]) > tolerance )
        {
            local_order[0] += adj[node_number][i] * std::sin(phases[i]);
            local_order[1] += adj[node_number][i] * std::cos(phases[i]);
            weight         += adj[node_number][i];
            ++counter;
        }
    }
    if ( weight == 0.0 )
        weight = counter;
    local_order[0] /= weight;
    local_order[1] /= weight;
    local_order[2]  = std::sqrt( local_order[0] * local_order[0] + local_order[1] * local_order[1] );
    return local_order;
}

inline dMatrix calculate_order_per_node( 
    const  dVec&    phases, 
    const  dMatrix& adj,
    double tolerance = 1e-12
)
{
    auto   N       = phases.size();
    double weight  = 0.0;
    size_t counter = 0;
    dMatrix Local_Order(N,{0.0,0.0,0.0});
    for ( size_t j = 0; j<N; ++j )
    {
        Local_Order[j] = calculate_order_per_node(phases,j,adj);
    }
    return Local_Order;
}

// Calculate the order parameter for a single node using a sparse adjacency matrix
inline dVec  calculate_order_per_node(
    const dVec& phases,
    size_t node,
    const SparsedMatrix& sparse_adj
)
{
    double total_weight = 0.0;
    dVec    local_order(3,0.0);
    for (const auto& entry : sparse_adj.rows[node])
    {
        size_t j      = entry.first;
        double weight = entry.second;
        local_order[0] += weight * std::sin(phases[j]);
        local_order[1] += weight * std::cos(phases[j]);
        total_weight   += weight;
    }
    if ( total_weight == 0.0 )
        total_weight = sparse_adj.rows[node].size();
    local_order[0] /= total_weight;
    local_order[1] /= total_weight;
    local_order[2]  = std::sqrt( local_order[0] * local_order[0] + local_order[1] * local_order[1] );
    return local_order;
}

// Calculate the order parameter for a single node using a sparse adjacency matrix
inline dMatrix calculate_order_per_node(
    const dVec& phases,
    const SparsedMatrix& sparse_adj
)
{
    auto   N            = phases.size();
    double total_weight = 0.0;
    size_t counter      = 0;
    dMatrix Local_Order(N,{0.0,0.0,0.0});
    for ( size_t j = 0; j<N; ++j )
    {
        Local_Order[j] = calculate_order_per_node(phases,j,sparse_adj);
    }
    return Local_Order;
}

// Calculate the order parameter for each module in a modular network
inline dVec  calculate_order_per_module( 
    const  dVec& phases, 
    const  wVec& module_assignments, 
    size_t       num_modules
)
{
    dVec   module_orders(num_modules*3 + 3, 0.0);
    wVec   counts(num_modules, 0);
    size_t total_count = 0;
    for ( size_t i = 0; i < phases.size(); ++i )
    {
        size_t m                = module_assignments[i];
        module_orders[3*m]     += std::sin(phases[i]);
        module_orders[3*m + 1] += std::cos(phases[i]);
        counts[m]++;
    }
    for ( size_t m = 0; m < num_modules; ++m )
    {
        if ( counts[m] > 0 )
        {
            module_orders[3*num_modules]     += module_orders[3*m];
            module_orders[3*num_modules + 1] += module_orders[3*m + 1];
            module_orders[3*m]               /= counts[m];
            module_orders[3*m + 1]           /= counts[m];
            module_orders[3*m + 2]            = std::sqrt( module_orders[3*m] * module_orders[3*m] + module_orders[3*m + 1] * module_orders[3*m + 1] );
            total_count                      += counts[m];
        }
    }
    module_orders[3*num_modules]     /= total_count;
    module_orders[3*num_modules + 1] /= total_count;
    module_orders[3*num_modules + 2] += std::sqrt( module_orders[3*num_modules] * module_orders[3*num_modules] + module_orders[3*num_modules + 1] * module_orders[3*num_modules + 1] );
    return module_orders;
}

// Calculate the per-module order parameter evolution for a simulation
// Returns: Vec<std::pair<Triple, double>> where each entry is (order parameters for all modules, time)
inline dMatrix calculate_order_per_module(
    const SolverResults& results,
    const wVec&          module_assignments,
    size_t               num_modules
)
{
    const auto& sol  = results.solution;
    const auto& time = results.timePoints;
    const auto  N    = sol[0].size();
    dVec         order_per_time(num_modules*3 + 3, 0.0);
    dMatrix      all_orders;
    for ( size_t t = 0; t < sol.size(); ++t )
    {
        order_per_time = calculate_order_per_module(sol[t], module_assignments, num_modules);
        order_per_time.push_back(time[t]);
        all_orders.push_back(order_per_time);
    }
    return all_orders;
}

// Optimized hierarchical order parameter calculation using aggregation from lower levels
// Returns: Vec<Triple> where result[l][g] is the order parameter for group g at level l (0=bottom/fine, last=top/global)
inline dMatrix calculate_order_hierarchical(
    const dVec&      phases,
    const Vec<wVec>& level_assignments,
    const wVec&      num_groups_per_level
)
{
    size_t       N          = phases.size();
    size_t       num_levels = level_assignments.size();
    Vec<dMatrix> components;
    Vec<wVec>    counts;
    // Isolating this block!
    {
    dMatrix component(num_groups_per_level[0], {0.0, 0.0}); 
    wVec    count(num_groups_per_level[0], 0);
    for ( size_t i = 0; i < N; ++i )
    {
        size_t g         = level_assignments[0][i];
        component[g][0] += std::sin(phases[i]);
        component[g][1] += std::cos(phases[i]);
        count[g]        += 1;
    }
    components.push_back(component);
    counts.push_back(count);
    } // End of the isolated block.
    // Higher levels: aggregate from lower level
    for ( size_t L = 1; L < num_levels; ++L )
    {
        dMatrix component(num_groups_per_level[L], {0.0, 0.0}); 
        wVec    count(num_groups_per_level[L], 0);
        auto    pre_components = components[L-1];
        auto    pre_counts     = counts[L-1];
        for ( size_t i = 0; i < N; ++i )
        {
            size_t g_lo          = level_assignments[L-1][i];
            size_t g_hi          = level_assignments[L][i];
            component[g_hi][0]  += pre_components[g_lo][0];
            component[g_hi][1]  += pre_components[g_lo][1];
            count[g_hi]          = pre_counts[g_lo];
            pre_components[g_lo] = {0.0, 0.0};            
            pre_counts[g_lo]     = 0;
        }
        components.push_back(component);
        counts.push_back(count);
    }

    // Isolating this block!
    {
    dMatrix component(1, {0.0, 0.0});
    wVec    count(1, 0);
    size_t L = num_levels - 1;
    for ( size_t i = 0; i < num_groups_per_level[L]; ++i )
    {
        component[0][0] += components[L][i][0];
        component[0][1] += components[L][i][1];
        count[0]        += counts[L][i];
    }
    components.push_back(component);
    counts.push_back(count);
    } // End of the isolated block.

    // Now, for each level, for each group, compute the normalized order parameter
    dMatrix level_orders;
    for ( size_t L = 0; L < num_levels + 1; ++L )
    {
        dVec    level_order(num_groups_per_level[L]*3, 0.0);
        for ( size_t g = 0; g < num_groups_per_level[L]; ++g )
        {
            double sum_sin = components[L][g][0];
            double sum_cos = components[L][g][1];
            size_t count   = counts[L][g];
            if ( count > 0 )
            {
                double mean_sin      = sum_sin / count;
                double mean_cos      = sum_cos / count;
                level_order[g*3]     = mean_sin;
                level_order[g*3 + 1] = mean_cos;
                level_order[g*3 + 2] = std::sqrt( mean_sin * mean_sin + mean_cos * mean_cos );
            }
        }
        level_orders.push_back(level_order);
    }
    return level_orders;
}
} // End MathEngine namespace

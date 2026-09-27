#pragma once
#include "../../typedefs/header.hpp"

namespace MathEngine
{ // namespace MathEngine

// -----------------------------------------------------------------------------
// Order-parameter helpers for the Kuramoto model.
//
// Every helper returns the (mean) sine/cosine of the phase ensemble plus the
// resulting magnitude (the Kuramoto order parameter), i.e. a 3-tuple
//   [ <sin theta>, <cos theta>, rho ],   rho = sqrt(<sin>^2 + <cos>^2).
// -----------------------------------------------------------------------------

// Order parameter of a single phase ensemble.
inline dVec calculate_order(const dVec& phases)
{
    const size_t N = phases.size();
    dVec once_order(3, 0.0);
    for (double theta : phases)
    {
        once_order[0] += std::sin(theta);
        once_order[1] += std::cos(theta);
    }
    once_order[0] /= static_cast<double>(N);
    once_order[1] /= static_cast<double>(N);
    once_order[2] = std::sqrt(once_order[0] * once_order[0] + once_order[1] * once_order[1]);
    return once_order;
}

// Order-parameter time series from solver results.
// Returns a matrix with columns [<sin>, <cos>, rho, time] (one row per time step).
inline dMatrix calculate_order(const SolverResults& results)
{
    const auto&  sol  = results.solution;
    const auto&  time = results.timePoints;
    const size_t rows = sol.Rows();

    if (rows != time.size())
    {
        throw std::invalid_argument(
            "[calculate_order] Dimensions mismatch: timePoints.size() (" +
            std::to_string(time.size()) + ") != solution.Rows() (" +
            std::to_string(rows) + ")");
    }

    dMatrix order(rows, 4, 0.0);
    for (size_t t = 0; t < rows; ++t)
    {
        const dVec row(sol[t].begin(), sol[t].end());
        const dVec o = calculate_order(row);
        order[t, 0] = o[0];
        order[t, 1] = o[1];
        order[t, 2] = o[2];
        order[t, 3] = time[t];
    }
    return order;
}

// Local (weighted) order parameter around a single node using a dense adjacency.
inline dVec calculate_order_per_node(
    const dVec&    phases,
    size_t         node_number,
    const dMatrix& adj,
    double         tolerance = 1e-12
)
{
    const size_t N = phases.size();
    double   weight  = 0.0;
    size_t   counter = 0;
    dVec     local_order(3, 0.0);

    for (size_t i = 0; i < N; ++i)
    {
        if (std::abs(adj[node_number, i]) > tolerance)
        {
            local_order[0] += adj[node_number, i] * std::sin(phases[i]);
            local_order[1] += adj[node_number, i] * std::cos(phases[i]);
            weight         += adj[node_number, i];
            ++counter;
        }
    }
    if (weight == 0.0) weight = counter;

    local_order[0] /= weight;
    local_order[1] /= weight;
    local_order[2]  = std::sqrt(local_order[0] * local_order[0] + local_order[1] * local_order[1]);
    return local_order;
}

// Per-node local order parameters using a dense adjacency.
// Returns an N x 3 matrix (rows are [<sin>, <cos>, rho] per node).
inline dMatrix calculate_order_per_node(
    const dVec&    phases,
    const dMatrix& adj,
    double         tolerance = 1e-12
)
{
    const size_t N = phases.size();
    dMatrix Local_Order(N, 3, 0.0);
    for (size_t j = 0; j < N; ++j)
    {
        const dVec local = calculate_order_per_node(phases, j, adj, tolerance);
        Local_Order[j, 0] = local[0];
        Local_Order[j, 1] = local[1];
        Local_Order[j, 2] = local[2];
    }
    return Local_Order;
}

// Local (weighted) order parameter around a single node using a sparse adjacency.
inline dVec calculate_order_per_node(
    const dVec&          phases,
    size_t               node,
    const SparsedMatrix& sparse_adj
)
{
    double total_weight = 0.0;
    dVec   local_order(3, 0.0);

    for (const auto& entry : sparse_adj.rows[node])
    {
        size_t j      = entry.first;
        double weight = entry.second;
        local_order[0] += weight * std::sin(phases[j]);
        local_order[1] += weight * std::cos(phases[j]);
        total_weight   += weight;
    }
    if (total_weight == 0.0) total_weight = sparse_adj.rows[node].size();

    local_order[0] /= total_weight;
    local_order[1] /= total_weight;
    local_order[2]  = std::sqrt(local_order[0] * local_order[0] + local_order[1] * local_order[1]);
    return local_order;
}

// Per-node local order parameters using a sparse adjacency.
// Returns an N x 3 matrix.
inline dMatrix calculate_order_per_node(
    const dVec&          phases,
    const SparsedMatrix& sparse_adj
)
{
    const size_t N = phases.size();
    dMatrix Local_Order(N, 3, 0.0);
    for (size_t j = 0; j < N; ++j)
    {
        const dVec local = calculate_order_per_node(phases, j, sparse_adj);
        Local_Order[j, 0] = local[0];
        Local_Order[j, 1] = local[1];
        Local_Order[j, 2] = local[2];
    }
    return Local_Order;
}

// Order parameter for each module in a modular network.
// Layout: per-module [<sin>, <cos>, rho] for modules 0..num_modules-1, then a
// final global [<sin>, <cos>, rho]; size = 3*num_modules + 3.
inline dVec calculate_order_per_module(
    const dVec&  phases,
    const wVec&  module_assignments,
    size_t       num_modules
)
{
    dVec   module_orders(num_modules * 3 + 3, 0.0);
    wVec   counts(num_modules, 0);
    size_t total_count = 0;

    for (size_t i = 0; i < phases.size(); ++i)
    {
        const size_t m = module_assignments[i];
        module_orders[3*m]     += std::sin(phases[i]);
        module_orders[3*m + 1] += std::cos(phases[i]);
        counts[m]++;
    }

    for (size_t m = 0; m < num_modules; ++m)
    {
        if (counts[m] > 0)
        {
            module_orders[3*num_modules]     += module_orders[3*m];
            module_orders[3*num_modules + 1] += module_orders[3*m + 1];

            module_orders[3*m]     /= static_cast<double>(counts[m]);
            module_orders[3*m + 1] /= static_cast<double>(counts[m]);
            module_orders[3*m + 2]  = std::sqrt(module_orders[3*m] * module_orders[3*m]
                                               + module_orders[3*m + 1] * module_orders[3*m + 1]);
            total_count += counts[m];
        }
    }

    module_orders[3*num_modules]     /= static_cast<double>(total_count);
    module_orders[3*num_modules + 1] /= static_cast<double>(total_count);
    module_orders[3*num_modules + 2]  = std::sqrt(module_orders[3*num_modules] * module_orders[3*num_modules]
                                                + module_orders[3*num_modules + 1] * module_orders[3*num_modules + 1]);
    return module_orders;
}

// Per-module order-parameter time series from solver results.
// Returns a matrix with columns [per-module triples..., global triple, time].
inline dMatrix calculate_order_per_module(
    const SolverResults& results,
    const wVec&          module_assignments,
    size_t               num_modules
)
{
    const auto&  sol  = results.solution;
    const auto&  time = results.timePoints;
    const size_t rows = sol.Rows();

    const size_t width = 3 * num_modules + 4;   // per-module triples + global triple + time
    dMatrix all_orders(rows, width, 0.0);

    for (size_t t = 0; t < rows; ++t)
    {
        const dVec row(sol[t].begin(), sol[t].end());
        dVec order_per_time = calculate_order_per_module(row, module_assignments, num_modules);
        order_per_time.push_back(time[t]);
        for (size_t k = 0; k < width; ++k)
            all_orders[t, k] = order_per_time[k];
    }
    return all_orders;
}

// Order parameter for each group at each level of a hierarchical network.
// result[L] (L in [0, num_levels)) holds the per-group [<sin>, <cos>, rho] triples
// at level L; result[num_levels] holds the single global [<sin>, <cos>, rho].
inline Vec<dVec> calculate_order_hierarchical(
    const dVec&      phases,
    const Vec<wVec>& level_assignments,
    const wVec&      num_groups_per_level
)
{
    const size_t N          = phases.size();
    const size_t num_levels = level_assignments.size();

    Vec<dVec> level_orders;
    level_orders.reserve(num_levels + 1);

    for (size_t L = 0; L < num_levels; ++L)
    {
        dVec level_order(num_groups_per_level[L] * 3, 0.0);
        wVec counts(num_groups_per_level[L], 0);

        for (size_t i = 0; i < N; ++i)
        {
            const size_t g = level_assignments[L][i];
            level_order[3*g]     += std::sin(phases[i]);
            level_order[3*g + 1] += std::cos(phases[i]);
            counts[g]++;
        }

        for (size_t g = 0; g < num_groups_per_level[L]; ++g)
        {
            if (counts[g] > 0)
            {
                const double s = level_order[3*g]     / static_cast<double>(counts[g]);
                const double c = level_order[3*g + 1] / static_cast<double>(counts[g]);
                level_order[3*g]     = s;
                level_order[3*g + 1] = c;
                level_order[3*g + 2] = std::sqrt(s * s + c * c);
            }
        }
        level_orders.push_back(std::move(level_order));
    }

    // Global order parameter across the whole ensemble.
    {
        double s = 0.0, c = 0.0;
        for (size_t i = 0; i < N; ++i)
        {
            s += std::sin(phases[i]);
            c += std::cos(phases[i]);
        }
        s /= static_cast<double>(N);
        c /= static_cast<double>(N);
        level_orders.push_back(dVec{ s, c, std::sqrt(s * s + c * c) });
    }

    return level_orders;
}

} // End MathEngine namespace

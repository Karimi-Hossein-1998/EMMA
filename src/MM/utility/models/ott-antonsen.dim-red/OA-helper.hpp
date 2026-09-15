#pragma once
#include "../../typedefs/header.hpp"

// Helpers for converting solver state to (rho,phi) pairs and building hierarchical
// layer representations. Uses project aliases (Vec, dVec) for consistency.

// Convert a solver state that stores interleaved Re/Im into a vector of (rho,phi)
inline dpairVec state_to_rhophi_from_reim(const dVec& state)
{
    size_t n = state.size() / 2;
    dpairVec out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        double re = state[2*i + 0];
        double im = state[2*i + 1];
        double rho = std::hypot(re, im);
        double phi = std::atan2(im, re);
        out.emplace_back(rho, phi);
    }
    return out;
}

// Convert a solver state that stores interleaved rho/phi into vector of (rho,phi)
inline dpairVec state_to_rhophi_from_rhophi(const dVec& state)
{
    size_t n = state.size() / 2;
    dpairVec out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        double rho = state[2*i + 0];
        double phi = state[2*i + 1];
        out.emplace_back(rho, phi);
    }
    return out;
}

// Normalize eta: if empty or sum == 0, return uniform; otherwise normalize to sum==1
inline dVec normalize_eta(const dVec& eta, size_t C)
{
    if (eta.size() != C) return dVec(C, 1.0 / static_cast<double>(C));
    double s = std::accumulate(eta.begin(), eta.end(), 0.0);
    if (s == 0.0) return dVec(C, 1.0 / static_cast<double>(C));
    dVec out = eta;
    for (double &v : out) v /= s;
    return out;
}

// (groups_at_level removed) contiguous-block folding is used instead.
// Build hierarchical layers of rho/phi pairs.
// - base_pairs: vector of (rho,phi) for the base (finest) level (size C)
// - num_modules: base number of modules at level 1
// - hierarchy_level: number of levels (>=1)
// - eta: community fractions for the base level (size C)
// - level (optional): if >=1 and <= hierarchy_level, return only that level (1-based)
// - grouping: CONTIGUOUS (default), ROUND_ROBIN, or EXPLICIT
// - explicit_map: if grouping==EXPLICIT, must be size C mapping each base index to a group index
// Contiguous-block hierarchical aggregation.
// At each step, consecutive blocks of size `num_modules` are folded together.
inline Vec<dpairVec> build_hierarchical_layers(
    const dpairVec& base_pairs,
    size_t          num_modules,
    size_t          hierarchy_level,
    const dVec&     eta,
    int             level = -1
)
{
    Vec<dpairVec> layers;
    if (hierarchy_level == 0) return layers;

    const size_t C         = base_pairs.size();
    dVec         eta_local = normalize_eta(eta, C);

    // Level 0 is the base level
    layers.emplace_back(base_pairs);

    // Iterative bottom-to-top aggregation: each level is built from the previous level
    // This mirrors the behavior in order.hpp where lower-level components are aggregated
    // to produce higher-level components.
    // Prepare base complex and weight vectors
    dcVec base_order_vecs; base_order_vecs.reserve(C);
    Vec<double> base_weights; base_weights.reserve(C);
    for (size_t c = 0; c < C; ++c)
    {
        base_order_vecs.emplace_back(std::polar(base_pairs[c].first, base_pairs[c].second));
        base_weights.push_back(eta_local[c]);
    }

    // Folding pattern: first fold groups of size `num_modules` (level 2), then
    // for higher levels fold consecutive pairs (size 2) repeatedly. This
    // mirrors order.hpp where after the first aggregation higher levels are
    // obtained by aggregating pairs of previous-level components.
    for (size_t L = 2; L <= hierarchy_level; ++L)
    {
        size_t prev_count = base_order_vecs.size();
        size_t block_size = (L == 2) ? num_modules : 2u;
        size_t groups = (prev_count + block_size - 1) / block_size; // ceil(prev_count/block_size)

        dpairVec this_level; this_level.reserve(groups);
        dcVec     this_cplx(groups, std::complex<double>(0.0, 0.0));
        dVec                          this_weights(groups, 0.0);

        for (size_t i = 0; i < prev_count; ++i)
        {
            size_t g = std::min(i / block_size, groups - 1);
            this_cplx[g] += base_order_vecs[i] * base_weights[i];
            this_weights[g] += base_weights[i];
        }

        for (size_t g = 0; g < groups; ++g)
        {
            if (this_weights[g] > 0.0) this_cplx[g] /= this_weights[g];
            double rho_agg = std::abs(this_cplx[g]);
            double phi_agg = std::atan2(this_cplx[g].imag(), this_cplx[g].real());
            this_level.emplace_back(rho_agg, phi_agg);
        }

        layers.emplace_back(std::move(this_level));
        base_order_vecs.swap(this_cplx);
        base_weights.swap(this_weights);
    }

    // Append single global aggregation (coarsest / total order)
    {
        std::complex<double> global_c(0.0, 0.0);
        double global_w = 0.0;
        for (size_t i = 0; i < base_order_vecs.size(); ++i) {
            global_c += base_order_vecs[i] * base_weights[i];
            global_w += base_weights[i];
        }
        if (global_w > 0.0) global_c /= global_w;
        double rho_g = std::abs(global_c);
        double phi_g = std::atan2(global_c.imag(), global_c.real());
        layers.emplace_back(dpairVec(1, {rho_g, phi_g}));
    }

    if (level >= 1 && static_cast<size_t>(level) <= layers.size()) {
        return { layers[level - 1] };
    }

    return layers;
}

// Overload: return hierarchical layers as complex r instead of (rho,phi)
inline dcMatrix build_hierarchical_layers_complex(
    const dpairVec& base_pairs,
    size_t          num_modules,
    size_t          hierarchy_level,
    const dVec&     eta,
    int             level = -1
)
{
    dcMatrix layers;
    if (hierarchy_level == 0) return layers;

    const size_t C = base_pairs.size();
    dVec eta_local = normalize_eta(eta, C);

    // base level
    dcVec base_cplx;
    base_cplx.reserve(C);
    Vec<double> base_weights; base_weights.reserve(C);
    for (size_t c = 0; c < C; ++c)
    {
        base_cplx.emplace_back(std::polar(base_pairs[c].first, base_pairs[c].second));
        base_weights.push_back(eta_local[c]);
    }
    layers.emplace_back(base_cplx);

    for (size_t L = 2; L <= hierarchy_level; ++L)
    {
        size_t prev_count = base_cplx.size();
        size_t block_size = (L == 2) ? num_modules : 2u;
        size_t groups = (prev_count + block_size - 1) / block_size; // ceil(prev_count/block_size)

        dcVec this_level(groups, std::complex<double>(0.0, 0.0));
        Vec<double> this_weights(groups, 0.0);

        for (size_t i = 0; i < prev_count; ++i)
        {
            size_t g = std::min(i / block_size, groups - 1);
            this_level[g] += base_cplx[i] * base_weights[i];
            this_weights[g] += base_weights[i];
        }

        for (size_t g = 0; g < groups; ++g)
        {
            if (this_weights[g] > 0.0) this_level[g] /= this_weights[g];
        }

        layers.emplace_back(this_level);
        base_cplx.swap(this_level);
        base_weights.swap(this_weights);
    }

    // Append final global mean-field (single complex) as top-most level
    {
        std::complex<double> global_c(0.0, 0.0);
        double global_w = 0.0;
        for (size_t i = 0; i < base_cplx.size(); ++i)
        {
            global_c += base_cplx[i] * base_weights[i];
            global_w += base_weights[i];
        }
        if (global_w > 0.0) global_c /= global_w;
        layers.emplace_back(dcVec(1, global_c));
    }

    if (level > 0 && static_cast<size_t>(level) <= layers.size())
    {
        return { layers[level - 1] };
    }

    return layers;
}

#pragma once
#include "../../typedefs/header.hpp"

// Optimized real-form RHS for Ott-Antonsen communities.
// State layout: interleaved [rho_0, phi_0, rho_1, phi_1, ...]
// This implementation uses pure real arithmetic (no std::complex or std::polar
// in the hot inner loop) and precomputes sin/cos of phases.

inline dVec ott_antonsen_expanded(
    double               time,
    const  dVec&         state,
    const  ott_antonsen& parameters,
    size_t               hierarchy_level = 0,
    size_t               num_modules     = 0
)
{
    const size_t C = parameters.gammas.size();
    dVec deriv(state.size(), 0.0);

    // Precompute Re/Im (a,b) and sin/cos for each community
    std::vector<double> a(C), b(C), cos_phi(C), sin_phi(C);
    for (size_t i = 0; i < C; ++i)
    {
        double rho = state[2*i + 0];
        double phi = state[2*i + 1];
        cos_phi[i] = std::cos(phi);
        sin_phi[i] = std::sin(phi);
        a[i] = rho * cos_phi[i];
        b[i] = rho * sin_phi[i];
    }

    auto K = [&](size_t i, size_t j) -> double
    {
        return (i < parameters.connections.size() && j < parameters.connections[i].size()) ? parameters.connections[i][j] : 0.0;
    };

    const double eps = 1e-14;

    for (size_t c = 0; c < C; ++c)
    {
        // compute coupling_sum = sum_{cp} K_{c,cp} ( r_c^2 conj(r_cp) - r_cp )
        // expand everything with real parts a,b and c,d for r_cp = c+id
        double ac = a[c], bc = b[c];

        // r_c^2 = (ac+ibc)^2 = (ac^2 - bc^2) + i(2*ac*bc)
        double rc2_re = ac*ac - bc*bc;
        double rc2_im = 2.0 * ac * bc;

        double sum_re = 0.0;
        double sum_im = 0.0;

        for (size_t cp = 0; cp < C; ++cp) {
            double weight = K(c, cp);
            double c_re = a[cp];
            double c_im = b[cp];

            // conj(r_cp) = c_re - i c_im
            // rc2 * conj(r_cp) = (rc2_re + i rc2_im) * (c_re - i c_im)
            // = (rc2_re*c_re + rc2_im*c_im) + i( rc2_im*c_re - rc2_re*c_im )
            double term1_re = rc2_re * c_re + rc2_im * c_im;
            double term1_im = rc2_im * c_re - rc2_re * c_im;

            // r_c^2 * conj(r_cp) - r_cp
            double pair_re = term1_re - c_re;
            double pair_im = term1_im - c_im;

            sum_re += weight * pair_re;
            sum_im += weight * pair_im;
        }

        // linear term: (-gamma + i mu) * r_c = (-gamma*ac - mu*bc) + i(-gamma*bc + mu*ac)
        double lin_re = -parameters.gammas[c] * ac - parameters.mus[c] * bc;
        double lin_im = -parameters.gammas[c] * bc + parameters.mus[c] * ac;

        double dRe = lin_re - sum_re;
        double dIm = lin_im - sum_im;

        // project to rho, phi
        double rho = std::sqrt(ac*ac + bc*bc);
        double inv_rho = 1.0 / std::max(rho, eps);

        // z = conj(r_c) * dr_c = (ac - i bc) * (dRe + i dIm)
        // Re(z) = ac*dRe + bc*dIm  = rho * d_rho
        // Im(z) = -bc*dRe + ac*dIm = rho^2 * d_phi
        double Re_z = ac * dRe + bc * dIm;
        double Im_z = -bc * dRe + ac * dIm;

        double d_rho = Re_z * inv_rho;
        double d_phi;
        if (rho > eps)
        {
            double inv_rho2 = inv_rho * inv_rho;
            d_phi = Im_z * inv_rho2; // Im_z / rho^2
        }
        else
        {
            // when rho is effectively zero, phase is undefined — keep phi static
            d_phi = 0.0;
        }

        deriv[2*c + 0] = d_rho;
        deriv[2*c + 1] = d_phi;
    }

    (void) hierarchy_level;
    (void) num_modules;

    return deriv;
}

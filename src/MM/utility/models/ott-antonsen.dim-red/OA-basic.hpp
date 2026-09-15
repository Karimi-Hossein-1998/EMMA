#pragma once
#include "../../typedefs/header.hpp"
#include <complex>

struct ott_antonsen
{
    dVec    gammas;
    dVec    mus;
    dMatrix connections;
};

inline dVec ott_antonsen_basic(
    double              time,
    const dVec&         state,
    const ott_antonsen& parameters
)
{
    const size_t C = parameters.gammas.size();
    dVec deriv( state.size(), 0.0 );

    // Helper lambda to access connection weight K_{c,c'} safely
    auto K = [&](size_t i, size_t j) -> double 
    {
        if ( i < parameters.connections.size() && j < parameters.connections[i].size() ) return parameters.connections[i][j];
        return 0.0;
    };

    for ( size_t c = 0; c < C; ++c )
    {
        // build complex r_c
        std::complex<double> r_c( state[2*c], state[2*c + 1] );

        // accumulate coupling term: sum_{c'} K_{c,c'} ( r_c^2 conj(r_{c'}) - r_{c'} )
        std::complex<double> coupling_sum(0.0, 0.0);
        for ( size_t cp = 0; cp < C; ++cp )
        {
            std::complex<double> r_cp( state[2*cp], state[2*cp + 1] );
            coupling_sum += K(c,cp) * ( r_c * r_c * std::conj(r_cp) - r_cp );
        }

        // local linear term: (-gamma_c + i mu_c) * r_c
        std::complex<double> linear_term( -parameters.gammas[c], parameters.mus[c] );
        std::complex<double> dr_c = linear_term * r_c - coupling_sum;

        deriv[2*c]     = std::real( dr_c );
        deriv[2*c + 1] = std::imag( dr_c );
    }

    return deriv;
}

inline MyFunc ott_antonsen_basic_wrapper(const ott_antonsen& params)
{
    return [params](double time, const dVec& state) -> dVec
    {
        return ott_antonsen_basic(time, state, params);
    };
}

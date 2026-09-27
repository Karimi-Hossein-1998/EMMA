#pragma once
#include "../typedefs/header.hpp"

namespace MathEngine
{ // namespace MathEngine

// ============================================================ //
//                                                              //
//  Ott-Antonsen (OA) Ansatz: dimensionality-reduced Kuramoto   //
//                                                              //
//  The full all-to-all Kuramoto model (an infinite-dimensional //
//  distribution of phases) collapses onto a single complex ODE //
//  for the order parameter r(t) when the intrinsic frequencies //
//  are Lorentzian-distributed. The complex state is stored in  //
//  Cartesian form (Re, Im) so the equations have no coordinate //
//  singularity (unlike the polar rho/phi form, whose phase     //
//  equation diverges as 1/rho as rho -> 0).                    //
//                                                              //
//  See OA-Ansatz.md for the derivation.                        //
//                                                              //
// ============================================================ //

// -----------------------------------------------------------------------------
// Single community (maps to the reduction of kuramoto-general with all-to-all
// adjacency).
//
//   dr/dt = (-gamma + i*mu) r - (K/2) (r^2 r* - r)
//         = (-gamma + i*mu) r + (K/2) (1 - |r|^2) r
//
// State layout: [x, y] = [Re(r), Im(r)].
// -----------------------------------------------------------------------------
inline void OA(
    double      time,
    const dVec& state,
    dVec&       dstate,
    double      gamma,   // Lorentzian half-width (decay rate)
    double      mu,      // Lorentzian center (mean intrinsic frequency)
    double      K        // Coupling strength
)
{
    const double x    = state[0];
    const double y    = state[1];
    const double rho2 = x*x + y*y;              // |r|^2
    const double drive = 0.5 * K * (1.0 - rho2);

    dstate[0] = -gamma*x - mu*y + drive*x;
    dstate[1] =  mu*x - gamma*y + drive*y;
    (void) time;
}

struct OAParams
{
    double gamma = 1.0; // Lorentzian half-width
    double mu    = 0.0; // Lorentzian center
    double K     = 1.0; // Coupling strength
};

inline MyFunc OA_wrapper(const OAParams& params)
{
    return [params](double t, const dVec& state, dVec& dstate) -> void
    {
        OA(t, state, dstate, params.gamma, params.mu, params.K);
    };
}

// -----------------------------------------------------------------------------
// Multi-community case.
//
//   dr_c/dt = (-gamma_c + i*mu_c) r_c
//           - (1/2) sum_{c'} [ K_{c,c'} * eta_{c'} ] ( r_c^2 conj(r_{c'}) - r_{c'} )
//
//   eta_{c'} = fraction of the population in community c'.
//
// State layout: interleaved [Re(r_0), Im(r_0), Re(r_1), Im(r_1), ...].
// -----------------------------------------------------------------------------
inline void OAGeneral(
    double        time,
    const dVec&   state,
    dVec&         dstate,
    const dVec&   gammas,   // per-community gamma (size C)
    const dVec&   mus,      // per-community mu    (size C)
    const dVec&   eta,      // per-community population fractions (size C)
    const dMatrix& K        // C x C coupling strengths K_{c,c'}
)
{
    const size_t C = gammas.size();

    for (size_t c = 0; c < C; ++c)
    {
        const double x = state[2*c + 0];
        const double y = state[2*c + 1];

        // r_c^2 = (x + i y)^2 = (x^2 - y^2) + i (2 x y)
        const double rc2_re = x*x - y*y;
        const double rc2_im = 2.0*x*y;

        // coupling_sum = sum_{c'} [ 0.5 * K_{c,c'} * eta_{c'} ] * ( r_c^2 conj(r_{c'}) - r_{c'} )
        double sum_re = 0.0;
        double sum_im = 0.0;
        for (size_t cp = 0; cp < C; ++cp)
        {
            const double w  = 0.5 * K[c, cp] * eta[cp];
            const double xp = state[2*cp + 0];
            const double yp = state[2*cp + 1];

            // Re( r_c^2 conj(r_{c'}) ) = rc2_re*xp + rc2_im*yp
            // Im( r_c^2 conj(r_{c'}) ) = rc2_im*xp - rc2_re*yp
            sum_re += w * (rc2_re*xp + rc2_im*yp - xp);
            sum_im += w * (rc2_im*xp - rc2_re*yp - yp);
        }

        // linear term: (-gamma + i mu) r_c = (-gamma*x - mu*y) + i (mu*x - gamma*y)
        dstate[2*c + 0] = -gammas[c]*x - mus[c]*y - sum_re;
        dstate[2*c + 1] =  mus[c]*x - gammas[c]*y - sum_im;
    }
    (void) time;
}

struct OAGeneralParams
{
    dVec    gammas; // per-community gamma (size C)
    dVec    mus;    // per-community mu    (size C)
    dVec    eta;    // per-community population fractions (size C)
    dMatrix K;      // C x C coupling strengths K_{c,c'}
    int     C = 0;  // number of communities
};

inline MyFunc OAGeneral_wrapper(const OAGeneralParams& params)
{
    return [params](double t, const dVec& state, dVec& dstate) -> void
    {
        OAGeneral(t, state, dstate, params.gammas, params.mus, params.eta, params.K);
    };
}

// -----------------------------------------------------------------------------
// Order parameter.
//
// The per-community order parameter is the state itself:
//   r_c = state[2c] + i state[2c+1],  rho_c = |r_c|.
// The global order parameter is the eta-weighted mean field:
//   R = sum_c eta_c r_c,  rho = |R|.
// -----------------------------------------------------------------------------

// Global order parameter, returns [Re(R), Im(R), rho].
// state: interleaved [Re(r_0), Im(r_0), ...]; eta: population fractions (size C).
inline dVec OAOrder(const dVec& state, const dVec& eta)
{
    const size_t C = state.size() / 2;
    double sum_re = 0.0;
    double sum_im = 0.0;
    for (size_t c = 0; c < C; ++c)
    {
        sum_re += eta[c] * state[2*c + 0];
        sum_im += eta[c] * state[2*c + 1];
    }

    dVec out(3, 0.0);
    out[0] = sum_re;
    out[1] = sum_im;
    out[2] = std::hypot(sum_re, sum_im);
    return out;
}

// Per-community order-parameter magnitudes rho_c, size C.
inline dVec OAOrderPerCommunity(const dVec& state)
{
    const size_t C = state.size() / 2;
    dVec rho(C, 0.0);
    for (size_t c = 0; c < C; ++c)
        rho[c] = std::hypot(state[2*c + 0], state[2*c + 1]);
    return rho;
}

// Global order parameter time series from solver results.
// Returns a matrix with columns [time, Re(R), Im(R), rho].
inline dMatrix OAOrder(const SolverResults& results, const dVec& eta)
{
    const auto&  sol  = results.solution;
    const auto&  time = results.timePoints;
    const size_t rows = sol.Rows();
    if (rows != time.size())
    {
        throw std::invalid_argument(
            "[OAOrder] Dimensions mismatch: timePoints.size() (" +
            std::to_string(time.size()) + ") != solution.Rows() (" +
            std::to_string(rows) + ")");
    }

    dMatrix out(rows, 4, 0.0);
    for (size_t t = 0; t < rows; ++t)
    {
        double re = 0.0;
        double im = 0.0;
        for (size_t c = 0; c < eta.size(); ++c)
        {
            re += eta[c] * sol[t][2*c + 0];
            im += eta[c] * sol[t][2*c + 1];
        }
        out[t, 0] = time[t];
        out[t, 1] = re;
        out[t, 2] = im;
        out[t, 3] = std::hypot(re, im);
    }
    return out;
}

} // namespace MathEngine

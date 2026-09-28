#pragma once
// -----------------------------------------------------------------------------
// Analysis & data gathering for the molecular-dynamics model.
//
// Provides the structural observables that play the role of the Kuramoto "order
// parameter" for a condensed-matter system:
//   - bond-orientational order psi4 / psi6 (solid / liquid / hexatic detection)
//   - radial distribution function g(r)
//   - mean-squared displacement (diffusion)
// plus the thermodynamic observables (energies, temperature, virial pressure).
// -----------------------------------------------------------------------------
#include <cmath>
#include <vector>

#include "engine.hpp"

namespace MathEngine
{

inline constexpr double kMDPi = 3.14159265358979323846;

// ----------------------------------------------------------------------------
// Bond-orientational order parameter.
//   psi_n = | (1/N) sum_i [ (1/N_i) sum_{j in nbr(i)} exp(i n theta_ij) ] |
// psi6 ~ 1 for a perfect triangular/hexatic crystal and ~ 0 for a liquid.
// ----------------------------------------------------------------------------
struct OrderParams { double psi4 = 0.0, psi6 = 0.0; };

inline OrderParams ComputeBondOrientationalOrder(const MolecularDynamics& md, double cutoff)
{
    const double rc  = cutoff * md.potentialParams.sigma;
    const double rc2 = rc * rc;
    const double W = md.width, H = md.height;
    const bool pbc = md.periodicBoundaryCondition;
    const size_t N = md.numParticles;

    double sumSin4 = 0.0, sumCos4 = 0.0, sumSin6 = 0.0, sumCos6 = 0.0;

    for (size_t i = 0; i < N; ++i)
    {
        double s4 = 0.0, c4 = 0.0, s6 = 0.0, c6 = 0.0;
        size_t nb = 0;
        for (size_t j = 0; j < N; ++j)
        {
            if (i == j) continue;
            double dx = md.posX[j] - md.posX[i];
            double dy = md.posY[j] - md.posY[i];
            if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); }
            const double r2 = dx * dx + dy * dy;
            if (r2 > rc2) continue;
            const double theta = std::atan2(dy, dx);
            s4 += std::sin(4.0 * theta); c4 += std::cos(4.0 * theta);
            s6 += std::sin(6.0 * theta); c6 += std::cos(6.0 * theta);
            ++nb;
        }
        if (nb > 0) { s4 /= static_cast<double>(nb); c4 /= static_cast<double>(nb);
                      s6 /= static_cast<double>(nb); c6 /= static_cast<double>(nb); }
        sumSin4 += s4; sumCos4 += c4;
        sumSin6 += s6; sumCos6 += c6;
    }

    const double invN = 1.0 / static_cast<double>(N);
    sumSin4 *= invN; sumCos4 *= invN; sumSin6 *= invN; sumCos6 *= invN;

    OrderParams op;
    op.psi4 = std::sqrt(sumSin4 * sumSin4 + sumCos4 * sumCos4);
    op.psi6 = std::sqrt(sumSin6 * sumSin6 + sumCos6 * sumCos6);
    return op;
}

// ----------------------------------------------------------------------------
// Radial distribution function g(r) over [0, rmax] with `nbins` bins.
// ----------------------------------------------------------------------------
struct RDFResult { std::vector<double> r; std::vector<double> g; };

inline RDFResult ComputeRDF(const MolecularDynamics& md, size_t nbins, double rmax)
{
    const double W = md.width, H = md.height;
    const bool pbc = md.periodicBoundaryCondition;
    const size_t N = md.numParticles;
    const double dr = rmax / static_cast<double>(nbins);

    std::vector<double> hist(nbins, 0.0);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = i + 1; j < N; ++j)
        {
            double dx = md.posX[j] - md.posX[i];
            double dy = md.posY[j] - md.posY[i];
            if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); }
            const double r = std::sqrt(dx * dx + dy * dy);
            if (r >= rmax) continue;
            const size_t b = static_cast<size_t>(r / dr);
            if (b < nbins) hist[b] += 1.0;
        }
    }

    const double rho = static_cast<double>(N) / md.Volume();
    RDFResult out;
    out.r.resize(nbins);
    out.g.resize(nbins);
    for (size_t k = 0; k < nbins; ++k)
    {
        const double r     = (static_cast<double>(k) + 0.5) * dr;
        const double shell = 2.0 * kMDPi * r * dr;
        const double ideal = 0.5 * static_cast<double>(N) * rho * shell;
        out.r[k] = r;
        out.g[k] = (ideal > 0.0) ? hist[k] / ideal : 0.0;
    }
    return out;
}

// ----------------------------------------------------------------------------
// Mean-squared displacement relative to the reference (t=0) configuration.
// Uses the minimum-image convention, so it is accurate while |dr| < L/2.
// ----------------------------------------------------------------------------
inline double ComputeMSD(const MolecularDynamics& md)
{
    const double W = md.width, H = md.height;
    const bool pbc = md.periodicBoundaryCondition;
    double sum = 0.0;
    for (size_t i = 0; i < md.numParticles; ++i)
    {
        double dx = md.posX[i] - md.initPosX[i];
        double dy = md.posY[i] - md.initPosY[i];
        if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); }
        sum += dx * dx + dy * dy;
    }
    return sum / static_cast<double>(md.numParticles);
}

// ----------------------------------------------------------------------------
// Convenience aggregate of everything at one instant.
// ----------------------------------------------------------------------------
struct Observables
{
    double time = 0.0;
    double temperature = 0.0;
    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;
    double totalEnergy = 0.0;
    double pressure = 0.0;
    double psi4 = 0.0;
    double psi6 = 0.0;
    double msd = 0.0;
};

inline Observables CollectObservables(MolecularDynamics& md, double time, double orderCutoff = 1.4)
{
    md.ComputeKineticEnergy();
    md.ComputePotentialEnergy();

    Observables o;
    o.time            = time;
    o.temperature     = md.currentTemperature;
    o.kineticEnergy   = md.kineticEnergy;
    o.potentialEnergy = md.potentialEnergy;
    o.totalEnergy     = md.kineticEnergy + md.potentialEnergy;
    o.pressure        = ComputePressure(md);
    o.msd             = ComputeMSD(md);

    const OrderParams op = ComputeBondOrientationalOrder(md, orderCutoff);
    o.psi4 = op.psi4;
    o.psi6 = op.psi6;
    return o;
}

} // namespace MathEngine

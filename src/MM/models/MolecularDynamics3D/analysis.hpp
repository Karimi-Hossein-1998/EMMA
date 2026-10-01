#pragma once
// -----------------------------------------------------------------------------
// Analysis for the 3D MD engine.
//
//   - Steinhardt bond-orientational order Q4 / Q6 (local, Lechner-Dellago)
//     computed via the bond-angle (Legendre polynomial) identity, which avoids
//     spherical harmonics: q_l(i) = sqrt( (1/N_i^2) sum_{j,k} P_l( r_j . r_k ) ).
//   - radial distribution function g(r) (3D normalisation, 4 pi r^2 dr shell)
//   - mean-squared displacement
//   - thermodynamic observables (energies, temperature, 3D virial pressure)
// -----------------------------------------------------------------------------
#include <cmath>
#include <vector>

#include "engine.hpp"
#include "thermostats.hpp"

namespace MathEngine
{

inline constexpr double kMDPi3D = 3.14159265358979323846;

inline double Legendre4(double x)
{
    const double x2 = x * x;
    return 0.125 * (35.0 * x2 * x2 - 30.0 * x2 + 3.0);
}
inline double Legendre6(double x)
{
    const double x2 = x * x;
    return 0.0625 * (231.0 * x2 * x2 * x2 - 315.0 * x2 * x2 + 105.0 * x2 - 5.0);
}

struct SteinhardtOrder { double q4 = 0.0, q6 = 0.0; };

inline SteinhardtOrder ComputeSteinhardtOrder3D(const MolecularDynamics3D& md, double cutoff)
{
    const double rc  = cutoff * md.potentialParams.sigma;
    const double rc2 = rc * rc;
    const double W = md.width, H = md.height, D = md.depth;
    const bool pbc = md.periodicBoundaryCondition;
    const size_t N = md.numParticles;

    double sumQ4 = 0.0, sumQ6 = 0.0;
    std::vector<double> ux, uy, uz;
    ux.reserve(32); uy.reserve(32); uz.reserve(32);

    for (size_t i = 0; i < N; ++i)
    {
        ux.clear(); uy.clear(); uz.clear();
        for (size_t j = 0; j < N; ++j)
        {
            if (i == j) continue;
            double dx = md.posX[j] - md.posX[i];
            double dy = md.posY[j] - md.posY[i];
            double dz = md.posZ[j] - md.posZ[i];
            if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); dz -= D * std::round(dz / D); }
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > rc2) continue;
            const double inv = 1.0 / std::sqrt(r2);
            ux.push_back(dx * inv); uy.push_back(dy * inv); uz.push_back(dz * inv);
        }

        const size_t nb = ux.size();
        if (nb == 0) continue;

        double s4 = 0.0, s6 = 0.0;
        for (size_t a = 0; a < nb; ++a)
            for (size_t b = 0; b < nb; ++b)
            {
                const double c = ux[a] * ux[b] + uy[a] * uy[b] + uz[a] * uz[b];
                s4 += Legendre4(c);
                s6 += Legendre6(c);
            }
        const double invNb2 = 1.0 / static_cast<double>(nb * nb);
        sumQ4 += std::sqrt(s4 * invNb2);
        sumQ6 += std::sqrt(s6 * invNb2);
    }

    SteinhardtOrder o;
    o.q4 = sumQ4 / static_cast<double>(N);
    o.q6 = sumQ6 / static_cast<double>(N);
    return o;
}

struct RDFResult3D { std::vector<double> r, g; };

inline RDFResult3D ComputeRDF3D(const MolecularDynamics3D& md, size_t nbins, double rmax)
{
    const double W = md.width, H = md.height, D = md.depth;
    const bool pbc = md.periodicBoundaryCondition;
    const size_t N = md.numParticles;
    const double dr = rmax / static_cast<double>(nbins);

    std::vector<double> hist(nbins, 0.0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = i + 1; j < N; ++j)
        {
            double dx = md.posX[j] - md.posX[i], dy = md.posY[j] - md.posY[i], dz = md.posZ[j] - md.posZ[i];
            if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); dz -= D * std::round(dz / D); }
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r >= rmax) continue;
            const size_t b = static_cast<size_t>(r / dr);
            if (b < nbins) hist[b] += 1.0;
        }

    const double rho = static_cast<double>(N) / md.Volume();
    RDFResult3D out;
    out.r.resize(nbins);
    out.g.resize(nbins);
    for (size_t k = 0; k < nbins; ++k)
    {
        const double r     = (static_cast<double>(k) + 0.5) * dr;
        const double shell = 4.0 * kMDPi3D * r * r * dr;
        const double ideal = 0.5 * static_cast<double>(N) * rho * shell;
        out.r[k] = r;
        out.g[k] = (ideal > 0.0) ? hist[k] / ideal : 0.0;
    }
    return out;
}

inline double ComputeMSD3D(const MolecularDynamics3D& md)
{
    const double W = md.width, H = md.height, D = md.depth;
    const bool pbc = md.periodicBoundaryCondition;
    double sum = 0.0;
    for (size_t i = 0; i < md.numParticles; ++i)
    {
        double dx = md.posX[i] - md.initPosX[i];
        double dy = md.posY[i] - md.initPosY[i];
        double dz = md.posZ[i] - md.initPosZ[i];
        if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); dz -= D * std::round(dz / D); }
        sum += dx * dx + dy * dy + dz * dz;
    }
    return sum / static_cast<double>(md.numParticles);
}

struct Observables3D
{
    double time = 0.0;
    double temperature = 0.0;
    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;
    double totalEnergy = 0.0;
    double pressure = 0.0;
    double q4 = 0.0;
    double q6 = 0.0;
    double msd = 0.0;
};

inline Observables3D CollectObservables3D(MolecularDynamics3D& md, double time, double orderCutoff = 1.4)
{
    md.ComputeKineticEnergy();
    md.ComputePotentialEnergy();

    Observables3D o;
    o.time = time;
    o.temperature = md.currentTemperature;
    o.kineticEnergy = md.kineticEnergy;
    o.potentialEnergy = md.potentialEnergy;
    o.totalEnergy = md.kineticEnergy + md.potentialEnergy;
    o.pressure = ComputePressure3D(md);
    o.msd = ComputeMSD3D(md);

    const SteinhardtOrder q = ComputeSteinhardtOrder3D(md, orderCutoff);
    o.q4 = q.q4;
    o.q6 = q.q6;
    return o;
}

} // namespace MathEngine

#pragma once
// -----------------------------------------------------------------------------
// Thermostats, a barostat, and pressure for the 3D MD engine.
// Same processes as the 2D version, extended to three velocity components and
// three translational DOF. The virial pressure coefficient is 1/3 (not 1/2).
// -----------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <random>

#include "engine.hpp"

namespace MathEngine
{

enum class ThermostatType3D
{
    None = 0,
    Rescale,
    Berendsen,
    Andersen,
    Langevin,
    NoseHoover
};

// 3D virial pressure: P = (N T - (1/3) virial) / V.
inline double ComputePressure3D(MolecularDynamics3D& md)
{
    md.ComputeKineticEnergy();
    const double virial = md.ComputeVirialSum();
    return (static_cast<double>(md.numParticles) * md.currentTemperature - (1.0 / 3.0) * virial)
           / md.Volume();
}

inline void ApplyVelocityRescale3D(MolecularDynamics3D& md, double T)
{
    md.ComputeKineticEnergy();
    if (md.currentTemperature > 1e-14) md.ScaleVelocities(std::sqrt(T / md.currentTemperature));
}

inline void ApplyBerendsen3D(MolecularDynamics3D& md, double T, double dt, double tau)
{
    md.ComputeKineticEnergy();
    const double lam = std::sqrt(1.0 + (dt / tau) * (T / md.currentTemperature - 1.0));
    md.ScaleVelocities(lam);
}

inline void ApplyAndersen3D(MolecularDynamics3D& md, double T, double dt, double nu, size_t& seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    const double p = std::min(1.0, nu * dt);
    for (size_t i = 0; i < md.numParticles; ++i)
    {
        if (ud(rng) < p)
        {
            const double s = std::sqrt(T / md.mass[i]);
            md.velX[i] = nd(rng) * s; md.velY[i] = nd(rng) * s; md.velZ[i] = nd(rng) * s;
        }
    }
    md.ComputeKineticEnergy();
    seed = rng();
}

inline void ApplyLangevin3D(MolecularDynamics3D& md, double T, double dt, double gamma, size_t& seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    const double f = std::exp(-gamma * dt);
    const double onef = 1.0 - f * f;
    for (size_t i = 0; i < md.numParticles; ++i)
    {
        const double ns = std::sqrt(onef * T / md.mass[i]);
        md.velX[i] = md.velX[i] * f + ns * nd(rng);
        md.velY[i] = md.velY[i] * f + ns * nd(rng);
        md.velZ[i] = md.velZ[i] * f + ns * nd(rng);
    }
    md.ComputeKineticEnergy();
    seed = rng();
}

inline void StepNoseHoover3D(MolecularDynamics3D& md, double dt, double T, double tau)
{
    const size_t N = md.numParticles;
    const double Ndf = 3.0 * static_cast<double>(N);
    md.nhQ = Ndf * T * tau * tau;
    const double dtHalf = 0.5 * dt;

    md.ComputeKineticEnergy();
    md.nhXi += dtHalf * (2.0 * md.kineticEnergy - Ndf * T) / md.nhQ;

    for (size_t i = 0; i < N; ++i)
    {
        md.velX[i] += dtHalf * (md.accX[i] - md.nhXi * md.velX[i]); md.posX[i] += dt * md.velX[i];
        md.velY[i] += dtHalf * (md.accY[i] - md.nhXi * md.velY[i]); md.posY[i] += dt * md.velY[i];
        md.velZ[i] += dtHalf * (md.accZ[i] - md.nhXi * md.velZ[i]); md.posZ[i] += dt * md.velZ[i];
    }
    md.ApplyBoundaries();
    md.CalculateAccelerations();
    for (size_t i = 0; i < N; ++i)
    {
        md.velX[i] += dtHalf * (md.accX[i] - md.nhXi * md.velX[i]);
        md.velY[i] += dtHalf * (md.accY[i] - md.nhXi * md.velY[i]);
        md.velZ[i] += dtHalf * (md.accZ[i] - md.nhXi * md.velZ[i]);
    }

    md.ComputeKineticEnergy();
    md.nhXi += dtHalf * (2.0 * md.kineticEnergy - Ndf * T) / md.nhQ;
}

inline void ApplyBerendsenBarostat3D(MolecularDynamics3D& md, double P_target, double dt, double tauP)
{
    const double P = ComputePressure3D(md);
    const double mu = std::pow(1.0 - (dt / tauP) * (P_target - P), 1.0 / 3.0);
    md.ScalePositionsAndBox(mu);
}

} // namespace MathEngine

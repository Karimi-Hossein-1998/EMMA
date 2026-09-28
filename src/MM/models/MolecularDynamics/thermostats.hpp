#pragma once
// -----------------------------------------------------------------------------
// Thermostats, a barostat, and pressure evaluation for the MD engine.
//
// These are "processes" applied on top of (or in place of) a plain NVE step:
//   - Rescale   : instantaneous velocity rescaling (also available as
//                 MolecularDynamics::SetTemperature).
//   - Berendsen : weak velocity coupling toward a target temperature.
//   - Andersen  : stochastic velocity reassignment (collision bath).
//   - Langevin  : friction + fluctuating force (FDT-consistent OU update).
//   - NoseHoover: extended-system thermostat integrated inside the step.
//   - BerendsenBarostat : weak pressure coupling (NPT-like).
// -----------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <random>

#include "engine.hpp"

namespace MathEngine
{

enum class ThermostatType
{
    None = 0,
    Rescale,
    Berendsen,
    Andersen,
    Langevin,
    NoseHoover
};

// ---- pressure ---------------------------------------------------------------
// Virial pressure in 2D: P = (N T - 0.5 * virialSum) / V.
// `virialSum` accumulates  r_ij . F_ij  with r_ij = r_j - r_i and F_ij the force
// ON i, which is NEGATIVE for repulsion; the virial theorem therefore enters
// with a minus sign so that repulsion raises the pressure.
inline double ComputePressure(MolecularDynamics& md)
{
    md.ComputeKineticEnergy();
    const double virial = md.ComputeVirialSum();
    return (static_cast<double>(md.numParticles) * md.currentTemperature - 0.5 * virial)
           / md.Volume();
}

// ---- thermostats ------------------------------------------------------------
inline void ApplyVelocityRescale(MolecularDynamics& md, double T)
{
    md.ComputeKineticEnergy();
    if (md.currentTemperature > 1e-14)
        md.ScaleVelocities(std::sqrt(T / md.currentTemperature));
}

inline void ApplyBerendsen(MolecularDynamics& md, double T, double dt, double tau)
{
    md.ComputeKineticEnergy();
    const double lam = std::sqrt(1.0 + (dt / tau) * (T / md.currentTemperature - 1.0));
    md.ScaleVelocities(lam);
}

inline void ApplyAndersen(MolecularDynamics& md, double T, double dt, double nu, size_t& seed)
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
            md.velX[i] = nd(rng) * s;
            md.velY[i] = nd(rng) * s;
        }
    }
    md.ComputeKineticEnergy();
    seed = rng();
}

// Langevin thermostat using the exact Ornstein-Uhlenbeck update for the random
// (free-particle) part, which is unconditionally stable for any friction.
inline void ApplyLangevin(MolecularDynamics& md, double T, double dt, double gamma, size_t& seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    const double f    = std::exp(-gamma * dt);
    const double onef = 1.0 - f * f;

    for (size_t i = 0; i < md.numParticles; ++i)
    {
        const double noiseScale = std::sqrt(onef * T / md.mass[i]);
        md.velX[i] = md.velX[i] * f + noiseScale * nd(rng);
        md.velY[i] = md.velY[i] * f + noiseScale * nd(rng);
    }
    md.ComputeKineticEnergy();
    seed = rng();
}

// Nosé-Hoover thermostat integrated with a velocity-Verlet scheme. The extended
// variable `xi` (with inertia Q = N_dof T tau^2) drives the kinetic temperature
// toward the target. Use this INSTEAD of MolecularDynamics::Step.
inline void StepNoseHoover(MolecularDynamics& md, double dt, double T, double tau)
{
    const size_t N = md.numParticles;
    const double Ndf = 2.0 * static_cast<double>(N);   // 2 translational DOF / particle
    md.nhQ = Ndf * T * tau * tau;
    const double dtHalf = 0.5 * dt;

    // half step of the thermostat variable
    md.ComputeKineticEnergy();
    md.nhXi += dtHalf * (2.0 * md.kineticEnergy - Ndf * T) / md.nhQ;

    // first half kick + drift
    for (size_t i = 0; i < N; ++i)
    {
        md.velX[i] += dtHalf * (md.accX[i] - md.nhXi * md.velX[i]);
        md.posX[i] += dt * md.velX[i];
        md.velY[i] += dtHalf * (md.accY[i] - md.nhXi * md.velY[i]);
        md.posY[i] += dt * md.velY[i];
    }
    md.ApplyBoundaries();
    md.CalculateAccelerations();

    // second half kick
    for (size_t i = 0; i < N; ++i)
    {
        md.velX[i] += dtHalf * (md.accX[i] - md.nhXi * md.velX[i]);
        md.velY[i] += dtHalf * (md.accY[i] - md.nhXi * md.velY[i]);
    }

    // full step of the thermostat variable
    md.ComputeKineticEnergy();
    md.nhXi += dtHalf * (2.0 * md.kineticEnergy - Ndf * T) / md.nhQ;
}

// ---- barostat ---------------------------------------------------------------
// Berendsen (weak) pressure coupling: scales positions and the box so the
// instantaneous pressure relaxes toward P_target over a time scale tauP.
inline void ApplyBerendsenBarostat(MolecularDynamics& md, double P_target, double dt, double tauP)
{
    const double P  = ComputePressure(md);
    const double mu = std::pow(1.0 - (dt / tauP) * (P_target - P), 0.5); // 2D linear scale
    md.ScalePositionsAndBox(mu);
}

} // namespace MathEngine

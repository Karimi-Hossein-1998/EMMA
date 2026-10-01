#pragma once
// -----------------------------------------------------------------------------
// Core 3D molecular-dynamics engine (Structure-of-Arrays, backend-independent).
//
// A separate mirror of the 2D `MolecularDynamics`: 3D positions/velocities/
// accelerations, a cuboid box, and symplectic velocity-Verlet / position-Verlet
// integrators. Temperature uses 3 translational DOF (T = 2 KE / (3 N)).
// Wall reflection uses the same folded ("periodic-reflection") map as the 2D
// engine, which handles arbitrarily large overshoots in O(1).
// -----------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "potential.hpp"
#include "initial-conditions.hpp"

namespace MathEngine
{

enum class IntegratorType3D
{
    VelocityVerlet = 0,
    Leapfrog
};

struct MDConfig3D
{
    size_t numParticles = 216;
    double width  = 11.0;
    double height = 11.0;
    double depth  = 11.0;

    double mass        = 1.0;
    double radius      = 0.1;
    double massRatio   = 2.0;
    double radiusRatio = 1.5;

    double sigma       = 1.0;
    double epsilon     = 1.0;
    double cutoffCoeff = 2.5;
    double morseAlpha  = 1.0;
    double powerN      = 9.0;
    double yukawaKappa = 1.0;

    double temperature   = 1.0;
    double restitution   = 1.0;
    double minSeparation = 0.8;
    size_t seed          = 41;

    bool periodicBoundaryCondition = true;
    bool bounce         = false;
    bool hardSphereCollisions = false;

    PotentialType3D        potential        = PotentialType3D::LennardJones;
    InitialConditionType3D initialCondition = InitialConditionType3D::SC;
};

class MolecularDynamics3D
{
public:
    std::vector<double> mass, massInv;
    std::vector<double> posX, posY, posZ;
    std::vector<double> velX, velY, velZ;
    std::vector<double> accX, accY, accZ;
    std::vector<double> radius;
    std::vector<int>    species;

    double width, height, depth;
    double restitution;
    bool   periodicBoundaryCondition, bounce, hardSphereCollisions;

    PotentialType3D   potentialType;
    PotentialParams3D potentialParams;

    size_t numParticles;
    double currentTemperature = 0.0;
    double kineticEnergy      = 0.0;
    double potentialEnergy    = 0.0;

    std::vector<double> initPosX, initPosY, initPosZ;

    double nhXi = 0.0, nhQ = 1.0;

    MolecularDynamics3D() = default;

    explicit MolecularDynamics3D(const MDConfig3D& cfg)
    {
        numParticles = cfg.numParticles;
        width = cfg.width; height = cfg.height; depth = cfg.depth;
        restitution = std::clamp(cfg.restitution, 0.0, 1.0);
        periodicBoundaryCondition = cfg.periodicBoundaryCondition;
        bounce = cfg.bounce;
        hardSphereCollisions = cfg.hardSphereCollisions;
        potentialType = cfg.potential;

        potentialParams.sigma = cfg.sigma;
        potentialParams.epsilon = cfg.epsilon;
        potentialParams.cutoffCoeff = cfg.cutoffCoeff;
        potentialParams.morseAlpha = cfg.morseAlpha;
        potentialParams.powerN = cfg.powerN;
        potentialParams.yukawaKappa = cfg.yukawaKappa;
        FinalizePotential3D(cfg.potential, potentialParams);

        InitialConditionParams3D icp;
        icp.N = cfg.numParticles;
        icp.width = cfg.width; icp.height = cfg.height; icp.depth = cfg.depth;
        icp.sigma = cfg.sigma;
        icp.mass = cfg.mass; icp.radius = cfg.radius;
        icp.massRatio = cfg.massRatio; icp.radiusRatio = cfg.radiusRatio;
        icp.temperature = cfg.temperature;
        icp.minSeparation = cfg.minSeparation;
        icp.seed = cfg.seed;
        icp.type = cfg.initialCondition;

        InitialState3D s = MakeInitialState3D(icp);
        posX = std::move(s.posX); posY = std::move(s.posY); posZ = std::move(s.posZ);
        velX = std::move(s.velX); velY = std::move(s.velY); velZ = std::move(s.velZ);
        mass = std::move(s.mass); radius = std::move(s.radius); species = std::move(s.species);

        massInv.resize(numParticles);
        accX.assign(numParticles, 0.0);
        accY.assign(numParticles, 0.0);
        accZ.assign(numParticles, 0.0);
        for (size_t i = 0; i < numParticles; ++i) massInv[i] = 1.0 / mass[i];

        initPosX = posX; initPosY = posY; initPosZ = posZ;

        nhQ = 3.0 * static_cast<double>(numParticles) * cfg.temperature;
        ComputeKineticEnergy();
        CalculateAccelerations();
    }

    void CalculateAccelerations()
    {
        const double W = width, H = height, D = depth;
        const bool pbc = periodicBoundaryCondition;
        const PotentialType3D ptype = potentialType;
        const PotentialParams3D pp = potentialParams;
        const double invW = 1.0 / W, invH = 1.0 / H, invD = 1.0 / D;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            const double px = posX[i], py = posY[i], pz = posZ[i];
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (size_t j = 0; j < numParticles; ++j)
            {
                if (i == j) continue;
                double dx = posX[j] - px, dy = posY[j] - py, dz = posZ[j] - pz;
                if (pbc)
                {
                    dx -= W * std::nearbyint(dx * invW);
                    dy -= H * std::nearbyint(dy * invH);
                    dz -= D * std::nearbyint(dz * invD);
                }
                const double r2 = dx * dx + dy * dy + dz * dz;
                double fFactor, pe;
                PairForceEnergy3D(r2, ptype, pp, fFactor, pe);
                ax += dx * fFactor; ay += dy * fFactor; az += dz * fFactor;
            }
            accX[i] = ax * massInv[i]; accY[i] = ay * massInv[i]; accZ[i] = az * massInv[i];
        }
    }

    void Step(double dt)
    {
        const double dtHalf = 0.5 * dt;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            velX[i] += dtHalf * accX[i]; posX[i] += dt * velX[i];
            velY[i] += dtHalf * accY[i]; posY[i] += dt * velY[i];
            velZ[i] += dtHalf * accZ[i]; posZ[i] += dt * velZ[i];
        }
        ApplyBoundaries();
        CalculateAccelerations();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            velX[i] += dtHalf * accX[i];
            velY[i] += dtHalf * accY[i];
            velZ[i] += dtHalf * accZ[i];
        }
    }

    void StepLeapfrog(double dt)
    {
        std::vector<double> oldAx(accX), oldAy(accY), oldAz(accZ);
        const double dt2 = dt * dt;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            posX[i] += velX[i] * dt + 0.5 * accX[i] * dt2;
            posY[i] += velY[i] * dt + 0.5 * accY[i] * dt2;
            posZ[i] += velZ[i] * dt + 0.5 * accZ[i] * dt2;
        }
        ApplyBoundaries();
        CalculateAccelerations();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            velX[i] += 0.5 * (oldAx[i] + accX[i]) * dt;
            velY[i] += 0.5 * (oldAy[i] + accY[i]) * dt;
            velZ[i] += 0.5 * (oldAz[i] + accZ[i]) * dt;
        }
    }

    void ApplyBoundaries()
    {
        if (periodicBoundaryCondition)
        {
            for (size_t i = 0; i < numParticles; ++i)
            {
                posX[i] -= width  * std::floor(posX[i] / width);
                posY[i] -= height * std::floor(posY[i] / height);
                posZ[i] -= depth  * std::floor(posZ[i] / depth);
            }
        }
        if (bounce) ResolveWallCollisions();
        if (hardSphereCollisions) ResolveParticleCollisions();
    }

    // Folded ("periodic-reflection") wall map: robust to large overshoots.
    void ResolveWallCollisions()
    {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            const double r = radius[i];

            const double loX = r, hiX = width - r, spanX = hiX - loX;
            if (spanX > 0.0)
            {
                const double period = 2.0 * spanX;
                double f = std::fmod(posX[i] - loX, period);
                if (f < 0.0) f += period;
                if (f > spanX) { f = period - f; velX[i] = -velX[i] * restitution; }
                posX[i] = loX + f;
            }
            else posX[i] = loX;

            const double loY = r, hiY = height - r, spanY = hiY - loY;
            if (spanY > 0.0)
            {
                const double period = 2.0 * spanY;
                double f = std::fmod(posY[i] - loY, period);
                if (f < 0.0) f += period;
                if (f > spanY) { f = period - f; velY[i] = -velY[i] * restitution; }
                posY[i] = loY + f;
            }
            else posY[i] = loY;

            const double loZ = r, hiZ = depth - r, spanZ = hiZ - loZ;
            if (spanZ > 0.0)
            {
                const double period = 2.0 * spanZ;
                double f = std::fmod(posZ[i] - loZ, period);
                if (f < 0.0) f += period;
                if (f > spanZ) { f = period - f; velZ[i] = -velZ[i] * restitution; }
                posZ[i] = loZ + f;
            }
            else posZ[i] = loZ;
        }
    }

    void ResolveParticleCollisions()
    {
        for (size_t i = 0; i < numParticles; ++i)
        {
            for (size_t j = i + 1; j < numParticles; ++j)
            {
                double dx = posX[j] - posX[i], dy = posY[j] - posY[i], dz = posZ[j] - posZ[i];
                if (periodicBoundaryCondition)
                {
                    dx -= width  * std::round(dx / width);
                    dy -= height * std::round(dy / height);
                    dz -= depth  * std::round(dz / depth);
                }
                const double r2 = dx * dx + dy * dy + dz * dz;
                const double minDist = radius[i] + radius[j];
                if (r2 > minDist * minDist) continue;
                const double dist = std::sqrt(r2);
                const double nx = dx / dist, ny = dy / dist, nz = dz / dist;
                const double overlap = minDist - dist;
                const double massInvSum = massInv[i] + massInv[j];

                posX[i] -= nx * overlap * massInv[i] / massInvSum;
                posY[i] -= ny * overlap * massInv[i] / massInvSum;
                posZ[i] -= nz * overlap * massInv[i] / massInvSum;
                posX[j] += nx * overlap * massInv[j] / massInvSum;
                posY[j] += ny * overlap * massInv[j] / massInvSum;
                posZ[j] += nz * overlap * massInv[j] / massInvSum;

                const double rvx = velX[j] - velX[i], rvy = velY[j] - velY[i], rvz = velZ[j] - velZ[i];
                const double nv = rvx * nx + rvy * ny + rvz * nz;
                if (nv >= 0.0) continue;
                const double impulse = -(1.0 + restitution) * nv / massInvSum;
                velX[i] -= impulse * massInv[i] * nx;
                velY[i] -= impulse * massInv[i] * ny;
                velZ[i] -= impulse * massInv[i] * nz;
                velX[j] += impulse * massInv[j] * nx;
                velY[j] += impulse * massInv[j] * ny;
                velZ[j] += impulse * massInv[j] * nz;
            }
        }
    }

    void ScaleVelocities(double s)
    {
        for (size_t i = 0; i < numParticles; ++i) { velX[i] *= s; velY[i] *= s; velZ[i] *= s; }
        ComputeKineticEnergy();
    }

    void ScalePositionsAndBox(double s)
    {
        for (size_t i = 0; i < numParticles; ++i) { posX[i] *= s; posY[i] *= s; posZ[i] *= s; }
        width *= s; height *= s; depth *= s;
    }

    void SetTemperature(double temperature)
    {
        ComputeKineticEnergy();
        if (currentTemperature > std::numeric_limits<double>::epsilon())
            ScaleVelocities(std::sqrt(temperature / currentTemperature));
    }

    double ComputeKineticEnergy()
    {
        kineticEnergy = 0.0;
        for (size_t i = 0; i < numParticles; ++i)
            kineticEnergy += 0.5 * mass[i] * (velX[i] * velX[i] + velY[i] * velY[i] + velZ[i] * velZ[i]);
        currentTemperature = 2.0 * kineticEnergy / (3.0 * static_cast<double>(numParticles));
        return kineticEnergy;
    }

    double ComputePotentialEnergy()
    {
        const double W = width, H = height, D = depth;
        const bool pbc = periodicBoundaryCondition;
        double pe = 0.0;
        for (size_t i = 0; i < numParticles; ++i)
            for (size_t j = i + 1; j < numParticles; ++j)
            {
                double dx = posX[j] - posX[i], dy = posY[j] - posY[i], dz = posZ[j] - posZ[i];
                if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); dz -= D * std::round(dz / D); }
                double fFactor, peij;
                PairForceEnergy3D(dx * dx + dy * dy + dz * dz, potentialType, potentialParams, fFactor, peij);
                pe += peij; (void)fFactor;
            }
        potentialEnergy = pe;
        return pe;
    }

    double ComputeVirialSum()
    {
        const double W = width, H = height, D = depth;
        const bool pbc = periodicBoundaryCondition;
        double virial = 0.0;
        for (size_t i = 0; i < numParticles; ++i)
            for (size_t j = i + 1; j < numParticles; ++j)
            {
                double dx = posX[j] - posX[i], dy = posY[j] - posY[i], dz = posZ[j] - posZ[i];
                if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); dz -= D * std::round(dz / D); }
                const double r2 = dx * dx + dy * dy + dz * dz;
                double fFactor, peij;
                PairForceEnergy3D(r2, potentialType, potentialParams, fFactor, peij);
                virial += r2 * fFactor; (void)peij;
            }
        return virial;
    }

    double Volume() const { return width * height * depth; }
};

} // namespace MathEngine

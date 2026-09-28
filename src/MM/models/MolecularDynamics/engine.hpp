#pragma once
// -----------------------------------------------------------------------------
// Core 2D molecular-dynamics engine.
//
// The engine owns all particle state in Structure-of-Arrays layout and is
// completely independent of any rendering backend. Interactions are dispatched
// through `potential.hpp`, initial conditions through `initial-conditions.hpp`,
// and the integration / thermostat / analysis pieces live in `thermostats.hpp`
// and `analysis.hpp` respectively.
//
// Integrators are symplectic (velocity-Verlet / position-Verlet), which keeps
// the energy approximately conserved in the NVE ensemble — this is why the
// engine is NOT fed through a generic RK/AB ODE solver.
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

enum class IntegratorType
{
    VelocityVerlet = 0,   // default, symplectic
    Leapfrog              // implemented as the equivalent position-Verlet
};

// All the knobs required to build a simulation. This is the library-level
// analogue of the `MDParams` struct the app exposes in its UI.
struct MDConfig
{
    size_t numParticles = 100;
    double width  = 11.0;   // reduced units: ~11 sigma -> ~0.8 density liquid
    double height = 11.0;

    double mass         = 1.0;
    double radius       = 0.1;   // hard-sphere / visual radius (species 0)
    double massRatio    = 2.0;   // binary species 1 = mass   * massRatio
    double radiusRatio  = 1.5;   // binary species 1 = radius * radiusRatio

    double sigma         = 1.0;
    double epsilon       = 1.0;
    double cutoffCoeff   = 2.5;
    double morseAlpha    = 1.0;
    double powerN        = 9.0;   // SoftSphere exponent
    double yukawaKappa   = 1.0;   // Yukawa inverse screening length

    double temperature   = 1.0;  // initial (and reference) temperature
    double restitution   = 1.0;
    double minSeparation = 0.8;  // Random IC (x sigma)
    size_t seed          = 41;

    bool periodicBoundaryCondition = true;
    bool bounce         = false;
    bool hardSphereCollisions = false;

    PotentialType        potential        = PotentialType::LennardJones;
    InitialConditionType initialCondition = InitialConditionType::SquareLattice;
};

class MolecularDynamics
{
public:
    // ---- Structure-of-Arrays state -----------------------------------------
    std::vector<double> mass;
    std::vector<double> massInv;
    std::vector<double> posX, posY;
    std::vector<double> velX, velY;
    std::vector<double> accX, accY;
    std::vector<double> radius;
    std::vector<int>    species;

    // ---- geometry ----------------------------------------------------------
    double width, height;
    double restitution;
    bool   periodicBoundaryCondition;
    bool   bounce;
    bool   hardSphereCollisions;

    // ---- potential ---------------------------------------------------------
    PotentialType   potentialType;
    PotentialParams potentialParams;

    // ---- bookkeeping -------------------------------------------------------
    size_t numParticles;
    double currentTemperature = 0.0;
    double kineticEnergy      = 0.0;
    double potentialEnergy    = 0.0;

    // ---- reference positions for MSD --------------------------------------
    std::vector<double> initPosX, initPosY;

    // ---- Nosé-Hoover thermostat state -------------------------------------
    double nhXi    = 0.0;   // thermostat "velocity" xi
    double nhQ     = 1.0;   // thermostat inertia Q = N_dof * T * tau^2

    MolecularDynamics() = default;

    explicit MolecularDynamics(const MDConfig& cfg)
    {
        numParticles           = cfg.numParticles;
        width                  = cfg.width;
        height                 = cfg.height;
        restitution            = std::clamp(cfg.restitution, 0.0, 1.0);
        periodicBoundaryCondition = cfg.periodicBoundaryCondition;
        bounce                 = cfg.bounce;
        hardSphereCollisions   = cfg.hardSphereCollisions;
        potentialType          = cfg.potential;

        potentialParams.sigma       = cfg.sigma;
        potentialParams.epsilon     = cfg.epsilon;
        potentialParams.cutoffCoeff = cfg.cutoffCoeff;
        potentialParams.morseAlpha  = cfg.morseAlpha;
        potentialParams.powerN      = cfg.powerN;
        potentialParams.yukawaKappa = cfg.yukawaKappa;
        FinalizePotential(cfg.potential, potentialParams);

        InitialConditionParams icp;
        icp.N             = cfg.numParticles;
        icp.width         = cfg.width;
        icp.height        = cfg.height;
        icp.sigma         = cfg.sigma;
        icp.mass          = cfg.mass;
        icp.radius        = cfg.radius;
        icp.massRatio     = cfg.massRatio;
        icp.radiusRatio   = cfg.radiusRatio;
        icp.temperature   = cfg.temperature;
        icp.minSeparation = cfg.minSeparation;
        icp.seed          = cfg.seed;
        icp.type          = cfg.initialCondition;

        InitialState s = MakeInitialState(icp);
        posX    = std::move(s.posX);
        posY    = std::move(s.posY);
        velX    = std::move(s.velX);
        velY    = std::move(s.velY);
        mass    = std::move(s.mass);
        radius  = std::move(s.radius);
        species = std::move(s.species);

        massInv.resize(numParticles);
        accX.assign(numParticles, 0.0);
        accY.assign(numParticles, 0.0);
        for (size_t i = 0; i < numParticles; ++i) massInv[i] = 1.0 / mass[i];

        initPosX = posX;
        initPosY = posY;

        nhQ = 2.0 * static_cast<double>(numParticles) * cfg.temperature;
        ComputeKineticEnergy();
        CalculateAccelerations();
    }

    // ---- force computation ------------------------------------------------
    void CalculateAccelerations()
    {
        const double W = width, H = height;
        const bool pbc = periodicBoundaryCondition;
        const PotentialType ptype = potentialType;
        const PotentialParams pp   = potentialParams;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            const double px = posX[i];
            const double py = posY[i];
            double ax = 0.0, ay = 0.0;

            for (size_t j = 0; j < numParticles; ++j)
            {
                if (i == j) continue;
                double dx = posX[j] - px;
                double dy = posY[j] - py;
                if (pbc)
                {
                    dx -= W * std::round(dx / W);
                    dy -= H * std::round(dy / H);
                }
                const double r2 = dx * dx + dy * dy;
                double fFactor, pe;
                PairForceEnergy(r2, ptype, pp, fFactor, pe);
                ax += dx * fFactor;
                ay += dy * fFactor;
            }
            accX[i] = ax * massInv[i];
            accY[i] = ay * massInv[i];
        }
    }

    // ---- integration -------------------------------------------------------
    void Step(double dt)
    {
        const double dtHalf = 0.5 * dt;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            velX[i] += dtHalf * accX[i];
            posX[i] += dt * velX[i];
            velY[i] += dtHalf * accY[i];
            posY[i] += dt * velY[i];
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
        }
    }

    void StepLeapfrog(double dt)
    {
        // Position-Verlet, mathematically identical to leapfrog.
        std::vector<double> oldAx(accX), oldAy(accY);

        const double dt2 = dt * dt;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            posX[i] += velX[i] * dt + 0.5 * accX[i] * dt2;
            posY[i] += velY[i] * dt + 0.5 * accY[i] * dt2;
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
        }
    }

    // ---- boundaries --------------------------------------------------------
    void ApplyBoundaries()
    {
        if (periodicBoundaryCondition)
        {
            for (size_t i = 0; i < numParticles; ++i)
            {
                posX[i] -= width  * std::floor(posX[i] / width);
                posY[i] -= height * std::floor(posY[i] / height);
            }
        }
        if (bounce) ResolveWallCollisions();
        if (hardSphereCollisions) ResolveParticleCollisions();
    }

    void ResolveWallCollisions()
    {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < numParticles; ++i)
        {
            const double r = radius[i];
            if (posX[i] - r < 0.0)      { posX[i] = r;        velX[i] = -velX[i] * restitution; }
            else if (posX[i] + r > width) { posX[i] = width - r; velX[i] = -velX[i] * restitution; }
            if (posY[i] - r < 0.0)      { posY[i] = r;         velY[i] = -velY[i] * restitution; }
            else if (posY[i] + r > height) { posY[i] = height - r; velY[i] = -velY[i] * restitution; }
        }
    }

    void ResolveParticleCollisions()
    {
        for (size_t i = 0; i < numParticles; ++i)
        {
            for (size_t j = i + 1; j < numParticles; ++j)
            {
                double dx = posX[j] - posX[i];
                double dy = posY[j] - posY[i];
                if (periodicBoundaryCondition)
                {
                    dx -= width  * std::round(dx / width);
                    dy -= height * std::round(dy / height);
                }
                const double r2 = dx * dx + dy * dy;
                const double minDist = radius[i] + radius[j];
                if (r2 > minDist * minDist) continue;

                const double dist = std::sqrt(r2);
                const double nx = dx / dist;
                const double ny = dy / dist;
                const double overlap = minDist - dist;
                const double massInvSum = massInv[i] + massInv[j];

                posX[i] -= nx * overlap * massInv[i] / massInvSum;
                posY[i] -= ny * overlap * massInv[i] / massInvSum;
                posX[j] += nx * overlap * massInv[j] / massInvSum;
                posY[j] += ny * overlap * massInv[j] / massInvSum;

                const double relVX = velX[j] - velX[i];
                const double relVY = velY[j] - velY[i];
                const double normalVel = relVX * nx + relVY * ny;
                if (normalVel >= 0.0) continue;

                const double impulse = -(1.0 + restitution) * normalVel / massInvSum;
                velX[i] -= impulse * massInv[i] * nx;
                velY[i] -= impulse * massInv[i] * ny;
                velX[j] += impulse * massInv[j] * nx;
                velY[j] += impulse * massInv[j] * ny;
            }
        }
    }

    // ---- velocity / position rescaling ------------------------------------
    void ScaleVelocities(double s)
    {
        for (size_t i = 0; i < numParticles; ++i) { velX[i] *= s; velY[i] *= s; }
        ComputeKineticEnergy();
    }

    void ScalePositionsAndBox(double s)
    {
        for (size_t i = 0; i < numParticles; ++i) { posX[i] *= s; posY[i] *= s; }
        width  *= s;
        height *= s;
    }

    void SetTemperature(double temperature)
    {
        ComputeKineticEnergy();
        if (currentTemperature > std::numeric_limits<double>::epsilon())
            ScaleVelocities(std::sqrt(temperature / currentTemperature));
    }

    // ---- energy / temperature ---------------------------------------------
    double ComputeKineticEnergy()
    {
        kineticEnergy = 0.0;
        for (size_t i = 0; i < numParticles; ++i)
            kineticEnergy += 0.5 * mass[i] * (velX[i] * velX[i] + velY[i] * velY[i]);
        currentTemperature = kineticEnergy / static_cast<double>(numParticles); // 2D: T = <m v^2>/2
        return kineticEnergy;
    }

    // Total pair potential energy (sum over unordered pairs). O(N^2).
    double ComputePotentialEnergy()
    {
        const double W = width, H = height;
        const bool pbc = periodicBoundaryCondition;
        double pe = 0.0;
        for (size_t i = 0; i < numParticles; ++i)
        {
            for (size_t j = i + 1; j < numParticles; ++j)
            {
                double dx = posX[j] - posX[i];
                double dy = posY[j] - posY[i];
                if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); }
                double fFactor, peij;
                PairForceEnergy(dx * dx + dy * dy, potentialType, potentialParams, fFactor, peij);
                pe += peij;
                (void)fFactor;
            }
        }
        potentialEnergy = pe;
        return pe;
    }

    // Virial sum  Sum_{i<j} r_ij . f_ij  (needed for the pressure).
    double ComputeVirialSum()
    {
        const double W = width, H = height;
        const bool pbc = periodicBoundaryCondition;
        double virial = 0.0;
        for (size_t i = 0; i < numParticles; ++i)
        {
            for (size_t j = i + 1; j < numParticles; ++j)
            {
                double dx = posX[j] - posX[i];
                double dy = posY[j] - posY[i];
                if (pbc) { dx -= W * std::round(dx / W); dy -= H * std::round(dy / H); }
                const double r2 = dx * dx + dy * dy;
                double fFactor, peij;
                PairForceEnergy(r2, potentialType, potentialParams, fFactor, peij);
                virial += r2 * fFactor;
                (void)peij;
            }
        }
        return virial;
    }

    double Volume() const { return width * height; }
};

} // namespace MathEngine

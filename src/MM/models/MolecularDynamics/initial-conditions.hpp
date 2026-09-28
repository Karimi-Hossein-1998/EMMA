#pragma once
// -----------------------------------------------------------------------------
// Initial-condition generators for the 2D molecular-dynamics model.
//
// Produces an `InitialState` (positions, velocities, masses, radii) for a box of
// size (width, height). Velocities are drawn from a Maxwell-Boltzmann
// distribution at the requested temperature, the centre-of-mass velocity is
// subtracted, and the velocities are rescaled so the kinetic temperature matches
// the target exactly.
// -----------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace MathEngine
{

enum class InitialConditionType
{
    SquareLattice    = 0,  // regular square lattice (single species)
    HexagonalLattice,      // regular triangular lattice (single species)
    Random,                // uniform random with a minimum separation
    TwoPhaseSlab,          // dense slab in the middle, rarefied vapour around it
    BinaryMixture          // two species on a lattice (different mass/radius)
};

struct InitialState
{
    std::vector<double> posX, posY;
    std::vector<double> velX, velY;
    std::vector<double> mass;
    std::vector<double> radius;
    std::vector<int>    species;   // 0-based species index (== 0 unless binary)
};

struct InitialConditionParams
{
    size_t N             = 100;
    double width         = 800.0;
    double height        = 600.0;
    double sigma         = 1.0;    // length scale used for lattice spacing
    double mass          = 1.0;    // species-0 mass
    double radius        = 0.1;    // species-0 hard-sphere radius (visual + collision)
    double massRatio     = 2.0;    // species-1 mass   = mass   * massRatio
    double radiusRatio   = 1.5;    // species-1 radius = radius * radiusRatio
    double temperature   = 1.0;    // target kinetic temperature
    double minSeparation = 0.8;    // Random IC: minimum pair separation (x sigma)
    size_t seed          = 41;
    InitialConditionType type = InitialConditionType::SquareLattice;
};

// Chooses a grid (cols x rows) that holds at least `N` cells and matches the
// box aspect ratio, so the lattice always fills the box without wrapping onto
// itself (which would otherwise create exact particle overlaps under PBC).
inline void GridDims(size_t N, double W, double H, size_t& cols, size_t& rows)
{
    cols = std::max<size_t>(1, static_cast<size_t>(std::llround(std::sqrt(static_cast<double>(N) * W / H))));
    rows = (N + cols - 1) / cols;
    if (rows == 0) rows = 1;
}

// Generates positions / masses / radii for the requested initial condition.
inline void GenerateConfiguration(const InitialConditionParams& p, InitialState& s)
{
    s.posX.assign(p.N, 0.0);
    s.posY.assign(p.N, 0.0);
    s.mass.assign(p.N, p.mass);
    s.radius.assign(p.N, p.radius);
    s.species.assign(p.N, 0);

    switch (p.type)
    {
        case InitialConditionType::HexagonalLattice:
        {
            // Triangular (hexagonal) lattice: alternate rows offset by half a cell.
            size_t cols, rows;
            GridDims(p.N, p.width, p.height, cols, rows);
            const double sx = p.width  / static_cast<double>(cols + 0.5); // room for the offset
            const double sy = p.height / static_cast<double>(rows);
            for (size_t i = 0; i < p.N; ++i)
            {
                const size_t row = i / cols;
                const size_t col = i % cols;
                s.posX[i] = (static_cast<double>(col) + 0.5) * sx + (row % 2 == 1 ? 0.5 * sx : 0.0);
                s.posY[i] = (static_cast<double>(row) + 0.5) * sy;
            }
            break;
        }
        case InitialConditionType::Random:
        {
            // Rejection-sampled uniform positions with a minimum separation.
            std::mt19937_64 rng(p.seed);
            std::uniform_real_distribution<double> dx(0.0, p.width);
            std::uniform_real_distribution<double> dy(0.0, p.height);
            const double min2 = p.minSeparation * p.minSeparation * p.sigma * p.sigma;
            for (size_t i = 0; i < p.N; ++i)
            {
                bool ok = false;
                for (size_t attempt = 0; attempt < 10000 && !ok; ++attempt)
                {
                    const double x = dx(rng);
                    const double y = dy(rng);
                    ok = true;
                    for (size_t j = 0; j < i; ++j)
                    {
                        const double ddx = s.posX[j] - x;
                        const double ddy = s.posY[j] - y;
                        if (ddx * ddx + ddy * ddy < min2) { ok = false; break; }
                    }
                    if (ok) { s.posX[i] = x; s.posY[i] = y; }
                }
                // If we could not place it after many attempts, just accept.
                if (!ok) { s.posX[i] = dx(rng); s.posY[i] = dy(rng); }
            }
            break;
        }
        case InitialConditionType::TwoPhaseSlab:
        {
            // A dense band in the middle (30% of the width) + sparse surroundings.
            std::mt19937_64 rng(p.seed);
            std::uniform_real_distribution<double> dx(0.0, p.width);
            std::uniform_real_distribution<double> dy(0.0, p.height);
            const double slabL = 0.3 * p.width;
            const double slabR = 0.7 * p.width;
            const size_t nSlab = static_cast<size_t>(0.7 * p.N);
            size_t placed = 0;
            while (placed < nSlab)
            {
                const double x = slabL + (slabR - slabL) * (double(rng()) / rng.max());
                const double y = dy(rng);
                s.posX[placed] = x; s.posY[placed] = y; ++placed;
            }
            while (placed < p.N)
            {
                s.posX[placed] = dx(rng); s.posY[placed] = dy(rng); ++placed;
            }
            break;
        }
        case InitialConditionType::BinaryMixture:
        {
            // Two interleaved species on a square lattice.
            size_t cols, rows;
            GridDims(p.N, p.width, p.height, cols, rows);
            const double sx = p.width  / static_cast<double>(cols);
            const double sy = p.height / static_cast<double>(rows);
            for (size_t i = 0; i < p.N; ++i)
            {
                const size_t row = i / cols;
                const size_t col = i % cols;
                s.posX[i] = (static_cast<double>(col) + 0.5) * sx;
                s.posY[i] = (static_cast<double>(row) + 0.5) * sy;
                const int sp = static_cast<int>(i % 2);
                s.species[i] = sp;
                if (sp == 1)
                {
                    s.mass[i]   = p.mass * p.massRatio;
                    s.radius[i] = p.radius * p.radiusRatio;
                }
            }
            break;
        }
        case InitialConditionType::SquareLattice:
        default:
        {
            size_t cols, rows;
            GridDims(p.N, p.width, p.height, cols, rows);
            const double sx = p.width  / static_cast<double>(cols);
            const double sy = p.height / static_cast<double>(rows);
            for (size_t i = 0; i < p.N; ++i)
            {
                const size_t row = i / cols;
                const size_t col = i % cols;
                s.posX[i] = (static_cast<double>(col) + 0.5) * sx;
                s.posY[i] = (static_cast<double>(row) + 0.5) * sy;
            }
            break;
        }
    }
}

// Assigns Maxwell-Boltzmann velocities at `temperature`, subtracts the
// centre-of-mass velocity and rescales to the target temperature exactly.
inline void AssignVelocities(InitialState& s, double temperature, size_t seed)
{
    const size_t N = s.posX.size();
    s.velX.assign(N, 0.0);
    s.velY.assign(N, 0.0);

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);

    double totalPX = 0.0, totalPY = 0.0, totalMass = 0.0, ke = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        const double sigmaV = std::sqrt(temperature / s.mass[i]);
        s.velX[i] = nd(rng) * sigmaV;
        s.velY[i] = nd(rng) * sigmaV;
        totalPX += s.mass[i] * s.velX[i];
        totalPY += s.mass[i] * s.velY[i];
        totalMass += s.mass[i];
    }

    // Remove the centre-of-mass drift.
    const double comVX = totalPX / totalMass;
    const double comVY = totalPY / totalMass;
    for (size_t i = 0; i < N; ++i)
    {
        s.velX[i] -= comVX;
        s.velY[i] -= comVY;
        ke += 0.5 * s.mass[i] * (s.velX[i] * s.velX[i] + s.velY[i] * s.velY[i]);
    }

    // Rescale to the requested temperature: T = <m v^2>/2 in 2D, so T = KE/N.
    const double currentT = ke / static_cast<double>(N);
    if (currentT > std::numeric_limits<double>::epsilon())
    {
        const double scale = std::sqrt(temperature / currentT);
        for (size_t i = 0; i < N; ++i) { s.velX[i] *= scale; s.velY[i] *= scale; }
    }
}

inline InitialState MakeInitialState(const InitialConditionParams& p)
{
    InitialState s;
    GenerateConfiguration(p, s);
    AssignVelocities(s, p.temperature, p.seed);
    return s;
}

} // namespace MathEngine

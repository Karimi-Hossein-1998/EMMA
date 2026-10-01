#pragma once
// -----------------------------------------------------------------------------
// Initial-condition generators for the 3D molecular-dynamics model.
//
// Lattice types: SC (simple cubic), BCC (body-centred), FCC (face-centred),
// Random (rejection-sampled with minimum separation), Slab (dense region) and
// Binary (two species on a simple-cubic lattice).
// Velocities are Maxwell-Boltzmann at the requested temperature, the COM is
// removed and they are rescaled so T = 2 KE / (3 N) exactly (3 translational DOF).
// -----------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace MathEngine
{

enum class InitialConditionType3D
{
    SC     = 0,
    BCC    = 1,
    FCC    = 2,
    Random = 3,
    Slab   = 4,
    Binary = 5
};

struct InitialState3D
{
    std::vector<double> posX, posY, posZ;
    std::vector<double> velX, velY, velZ;
    std::vector<double> mass;
    std::vector<double> radius;
    std::vector<int>    species;
};

struct InitialConditionParams3D
{
    size_t N             = 216;
    double width         = 11.0;
    double height        = 11.0;
    double depth         = 11.0;
    double sigma         = 1.0;
    double mass          = 1.0;
    double radius        = 0.1;
    double massRatio     = 2.0;
    double radiusRatio   = 1.5;
    double temperature   = 1.0;
    double minSeparation = 0.8;
    size_t seed          = 41;
    InitialConditionType3D type = InitialConditionType3D::SC;
};

inline void GridDims3D(size_t ncells, double W, double H, double D,
                       size_t& nx, size_t& ny, size_t& nz)
{
    nx = std::max<size_t>(1, static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(ncells)))));
    ny = nx; nz = nx;
    while (nx * ny * nz < ncells)
    {
        const double rx = static_cast<double>(nx) / W;
        const double ry = static_cast<double>(ny) / H;
        const double rz = static_cast<double>(nz) / D;
        if (rx <= ry && rx <= rz) ++nx;
        else if (ry <= rz) ++ny;
        else ++nz;
    }
}

inline size_t AtomsPerCell3D(InitialConditionType3D t)
{
    return (t == InitialConditionType3D::BCC) ? 2 : (t == InitialConditionType3D::FCC ? 4 : 1);
}

inline void CellOffsets3D(InitialConditionType3D t, size_t a, double& fx, double& fy, double& fz)
{
    switch (t)
    {
        case InitialConditionType3D::BCC:
            if (a == 0) { fx = 0.0; fy = 0.0; fz = 0.0; }
            else        { fx = 0.5; fy = 0.5; fz = 0.5; }
            break;
        case InitialConditionType3D::FCC:
            if      (a == 0) { fx = 0.0; fy = 0.0; fz = 0.0; }
            else if (a == 1) { fx = 0.5; fy = 0.5; fz = 0.0; }
            else if (a == 2) { fx = 0.5; fy = 0.0; fz = 0.5; }
            else             { fx = 0.0; fy = 0.5; fz = 0.5; }
            break;
        case InitialConditionType3D::SC:
        case InitialConditionType3D::Binary:
        default:
            fx = 0.5; fy = 0.5; fz = 0.5;
            break;
    }
}

inline void GenerateConfiguration3D(const InitialConditionParams3D& p, InitialState3D& s)
{
    s.posX.assign(p.N, 0.0);
    s.posY.assign(p.N, 0.0);
    s.posZ.assign(p.N, 0.0);
    s.mass.assign(p.N, p.mass);
    s.radius.assign(p.N, p.radius);
    s.species.assign(p.N, 0);

    switch (p.type)
    {
        case InitialConditionType3D::SC:
        case InitialConditionType3D::BCC:
        case InitialConditionType3D::FCC:
        {
            const size_t apc = AtomsPerCell3D(p.type);
            size_t nx, ny, nz;
            GridDims3D((p.N + apc - 1) / apc, p.width, p.height, p.depth, nx, ny, nz);
            const double sx = p.width / static_cast<double>(nx);
            const double sy = p.height / static_cast<double>(ny);
            const double sz = p.depth / static_cast<double>(nz);
            size_t idx = 0;
            for (size_t k = 0; k < nz && idx < p.N; ++k)
            for (size_t j = 0; j < ny && idx < p.N; ++j)
            for (size_t i = 0; i < nx && idx < p.N; ++i)
            for (size_t a = 0; a < apc && idx < p.N; ++a)
            {
                double fx, fy, fz;
                CellOffsets3D(p.type, a, fx, fy, fz);
                s.posX[idx] = (static_cast<double>(i) + fx) * sx;
                s.posY[idx] = (static_cast<double>(j) + fy) * sy;
                s.posZ[idx] = (static_cast<double>(k) + fz) * sz;
                ++idx;
            }
            break;
        }
        case InitialConditionType3D::Binary:
        {
            size_t nx, ny, nz;
            GridDims3D(p.N, p.width, p.height, p.depth, nx, ny, nz);
            const double sx = p.width / static_cast<double>(nx);
            const double sy = p.height / static_cast<double>(ny);
            const double sz = p.depth / static_cast<double>(nz);
            for (size_t i = 0; i < p.N; ++i)
            {
                const size_t k = i / (nx * ny);
                const size_t rem = i % (nx * ny);
                const size_t j = rem / nx;
                const size_t ii = rem % nx;
                s.posX[i] = (static_cast<double>(ii) + 0.5) * sx;
                s.posY[i] = (static_cast<double>(j)  + 0.5) * sy;
                s.posZ[i] = (static_cast<double>(k)  + 0.5) * sz;
                const int sp = static_cast<int>(i % 2);
                s.species[i] = sp;
                if (sp == 1) { s.mass[i] = p.mass * p.massRatio; s.radius[i] = p.radius * p.radiusRatio; }
            }
            break;
        }
        case InitialConditionType3D::Random:
        {
            std::mt19937_64 rng(p.seed);
            std::uniform_real_distribution<double> dx(0.0, p.width);
            std::uniform_real_distribution<double> dy(0.0, p.height);
            std::uniform_real_distribution<double> dz(0.0, p.depth);
            const double min2 = p.minSeparation * p.minSeparation * p.sigma * p.sigma;
            for (size_t i = 0; i < p.N; ++i)
            {
                bool ok = false;
                for (size_t attempt = 0; attempt < 10000 && !ok; ++attempt)
                {
                    const double x = dx(rng), y = dy(rng), z = dz(rng);
                    ok = true;
                    for (size_t j = 0; j < i; ++j)
                    {
                        const double ddx = s.posX[j] - x, ddy = s.posY[j] - y, ddz = s.posZ[j] - z;
                        if (ddx * ddx + ddy * ddy + ddz * ddz < min2) { ok = false; break; }
                    }
                    if (ok) { s.posX[i] = x; s.posY[i] = y; s.posZ[i] = z; }
                }
                if (!ok) { s.posX[i] = dx(rng); s.posY[i] = dy(rng); s.posZ[i] = dz(rng); }
            }
            break;
        }
        case InitialConditionType3D::Slab:
        {
            std::mt19937_64 rng(p.seed);
            std::uniform_real_distribution<double> dx(0.0, p.width);
            std::uniform_real_distribution<double> dy(0.0, p.height);
            std::uniform_real_distribution<double> dz(0.0, p.depth);
            const double slabL = 0.3 * p.width, slabR = 0.7 * p.width;
            const size_t nSlab = static_cast<size_t>(0.7 * p.N);
            const double urd = static_cast<double>(rng.max() - rng.min());
            for (size_t i = 0; i < p.N; ++i)
            {
                if (i < nSlab)
                {
                    s.posX[i] = slabL + (slabR - slabL) * (static_cast<double>(rng() - rng.min()) / urd);
                    s.posY[i] = dy(rng); s.posZ[i] = dz(rng);
                }
                else { s.posX[i] = dx(rng); s.posY[i] = dy(rng); s.posZ[i] = dz(rng); }
            }
            break;
        }
    }
}

inline void AssignVelocities3D(InitialState3D& s, double temperature, size_t seed)
{
    const size_t N = s.posX.size();
    s.velX.assign(N, 0.0);
    s.velY.assign(N, 0.0);
    s.velZ.assign(N, 0.0);

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);

    double totalPX = 0.0, totalPY = 0.0, totalPZ = 0.0, totalMass = 0.0, ke = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        const double sv = std::sqrt(temperature / s.mass[i]);
        s.velX[i] = nd(rng) * sv; s.velY[i] = nd(rng) * sv; s.velZ[i] = nd(rng) * sv;
        totalPX += s.mass[i] * s.velX[i];
        totalPY += s.mass[i] * s.velY[i];
        totalPZ += s.mass[i] * s.velZ[i];
        totalMass += s.mass[i];
    }
    const double comX = totalPX / totalMass, comY = totalPY / totalMass, comZ = totalPZ / totalMass;
    for (size_t i = 0; i < N; ++i)
    {
        s.velX[i] -= comX; s.velY[i] -= comY; s.velZ[i] -= comZ;
        ke += 0.5 * s.mass[i] * (s.velX[i] * s.velX[i] + s.velY[i] * s.velY[i] + s.velZ[i] * s.velZ[i]);
    }

    const double currentT = 2.0 * ke / (3.0 * static_cast<double>(N));
    if (currentT > std::numeric_limits<double>::epsilon())
    {
        const double scale = std::sqrt(temperature / currentT);
        for (size_t i = 0; i < N; ++i) { s.velX[i] *= scale; s.velY[i] *= scale; s.velZ[i] *= scale; }
    }
}

inline InitialState3D MakeInitialState3D(const InitialConditionParams3D& p)
{
    InitialState3D s;
    GenerateConfiguration3D(p, s);
    AssignVelocities3D(s, p.temperature, p.seed);
    return s;
}

} // namespace MathEngine

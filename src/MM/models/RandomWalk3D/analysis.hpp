#pragma once
// -----------------------------------------------------------------------------
// Measurements for the 3D random-walk model.
// -----------------------------------------------------------------------------
#include <cmath>
#include <cstddef>

#include "engine.hpp"

namespace MathEngine
{

struct WalkerObservables3D
{
    std::uint64_t step = 0;
    double time = 0.0;

    double meanX  = 0.0, meanY  = 0.0, meanZ  = 0.0;
    double meanX2 = 0.0, meanY2 = 0.0, meanZ2 = 0.0;
    double meanXY = 0.0, meanXZ = 0.0, meanYZ = 0.0;
    double varX   = 0.0, varY   = 0.0, varZ   = 0.0;
    double covXY  = 0.0, covXZ  = 0.0, covYZ  = 0.0;

    double meanR     = 0.0;
    double rmsR      = 0.0;
    double msd       = 0.0;
    double diffusion = 0.0; // D = MSD / (6 t)  (3D)
    double radiusOfGyration = 0.0;
};

inline WalkerObservables3D CollectObservables3D(const RandomWalk3D& rw, double time)
{
    WalkerObservables3D o;
    o.step = rw.stepCount;
    o.time = time;

    const size_t N = rw.numWalkers;
    if (N == 0) return o;

    double sx = 0, sy = 0, sz = 0;
    double sx2 = 0, sy2 = 0, sz2 = 0;
    double sxy = 0, sxz = 0, syz = 0;
    double sr = 0, sd = 0;

    for (size_t i = 0; i < N; ++i)
    {
        const double x = rw.posX[i], y = rw.posY[i], z = rw.posZ[i];
        const double dx = x - rw.initX[i], dy = y - rw.initY[i], dz = z - rw.initZ[i];

        sx += x;  sy += y;  sz += z;
        sx2 += x * x; sy2 += y * y; sz2 += z * z;
        sxy += x * y; sxz += x * z; syz += y * z;
        sr += std::sqrt(dx * dx + dy * dy + dz * dz);
        sd += dx * dx + dy * dy + dz * dz;
    }

    const double invN = 1.0 / static_cast<double>(N);
    o.meanX  = sx * invN;  o.meanY  = sy * invN;  o.meanZ  = sz * invN;
    o.meanX2 = sx2 * invN; o.meanY2 = sy2 * invN; o.meanZ2 = sz2 * invN;
    o.meanXY = sxy * invN; o.meanXZ = sxz * invN; o.meanYZ = syz * invN;
    o.varX   = o.meanX2 - o.meanX * o.meanX;
    o.varY   = o.meanY2 - o.meanY * o.meanY;
    o.varZ   = o.meanZ2 - o.meanZ * o.meanZ;
    o.covXY  = o.meanXY - o.meanX * o.meanY;
    o.covXZ  = o.meanXZ - o.meanX * o.meanZ;
    o.covYZ  = o.meanYZ - o.meanY * o.meanZ;
    o.meanR  = sr * invN;
    o.msd    = sd * invN;
    o.rmsR   = std::sqrt(o.msd);
    o.diffusion = (time > 0.0) ? o.msd / (6.0 * time) : 0.0;
    o.radiusOfGyration = o.meanX2 + o.meanY2 + o.meanZ2
                       - o.meanX * o.meanX - o.meanY * o.meanY - o.meanZ * o.meanZ;
    return o;
}

} // namespace MathEngine

#pragma once
// -----------------------------------------------------------------------------
// Measurements for the random-walk model.
//
// Ensemble averages over the N walkers at a given time. Position moments are in
// canvas coordinates (walkers start at the canvas centre); the displacement
// based quantities (meanR, rmsR, MSD, D, R_g) are origin-independent and are the
// physically meaningful diffusion observables.
// -----------------------------------------------------------------------------
#include <cmath>
#include <cstddef>

#include "engine.hpp"

namespace MathEngine
{

struct WalkerObservables
{
    std::uint64_t step = 0;
    double time = 0.0;

    // First and second moments of the position (canvas coordinates).
    double meanX  = 0.0;   // <x>
    double meanY  = 0.0;   // <y>
    double meanX2 = 0.0;   // <x^2>
    double meanY2 = 0.0;   // <y^2>
    double meanXY = 0.0;   // <xy>
    double varX   = 0.0;   // <x^2> - <x>^2
    double varY   = 0.0;   // <y^2> - <y>^2
    double covXY  = 0.0;   // <xy> - <x><y>

    // Distance / diffusion observables (relative to the starting point).
    double meanR     = 0.0; // <|r|>,   r = displacement from start
    double rmsR      = 0.0; // sqrt(<|r|^2>) == sqrt(MSD)
    double msd       = 0.0; // <|r|^2>
    double diffusion = 0.0; // D = MSD / (4 t)   (2D)
    double radiusOfGyration = 0.0; // <x^2+y^2> - <x>^2 - <y>^2
};

inline WalkerObservables CollectObservables(const RandomWalk& rw, double time)
{
    WalkerObservables o;
    o.step = rw.stepCount;
    o.time = time;

    const size_t N = rw.numWalkers;
    if (N == 0) return o;

    double sx = 0.0, sy = 0.0, sx2 = 0.0, sy2 = 0.0, sxy = 0.0;
    double sr = 0.0, sd = 0.0;

    for (size_t i = 0; i < N; ++i)
    {
        const double x = rw.posX[i];
        const double y = rw.posY[i];
        const double dx = x - rw.initX[i];
        const double dy = y - rw.initY[i];

        sx  += x;  sy  += y;
        sx2 += x * x; sy2 += y * y; sxy += x * y;
        sr  += std::sqrt(dx * dx + dy * dy);
        sd  += dx * dx + dy * dy;
    }

    const double invN = 1.0 / static_cast<double>(N);
    o.meanX  = sx * invN;
    o.meanY  = sy * invN;
    o.meanX2 = sx2 * invN;
    o.meanY2 = sy2 * invN;
    o.meanXY = sxy * invN;
    o.varX   = o.meanX2 - o.meanX * o.meanX;
    o.varY   = o.meanY2 - o.meanY * o.meanY;
    o.covXY  = o.meanXY - o.meanX * o.meanY;
    o.meanR  = sr * invN;
    o.msd    = sd * invN;
    o.rmsR   = std::sqrt(o.msd);
    o.diffusion = (time > 0.0) ? o.msd / (4.0 * time) : 0.0;
    o.radiusOfGyration = o.meanX2 + o.meanY2 - o.meanX * o.meanX - o.meanY * o.meanY;
    return o;
}

} // namespace MathEngine

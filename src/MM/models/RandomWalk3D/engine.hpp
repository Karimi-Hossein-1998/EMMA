#pragma once
// -----------------------------------------------------------------------------
// Random-walk 3D engine (Structure-of-Arrays, backend-independent).
//
// A separate mirror of the 2D `RandomWalk`: 3D positions, a 3D direction set
// (straight / plane-diagonal / body-diagonal / full and their combinations,
// plus "with center" and continuous variants) and a cuboid boundary. Positions
// are tracked *unwrapped*; the wrapped coordinate is produced via `WrapX/Y/Z`.
// -----------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "../RandomWalk/engine.hpp" // for MathEngine::BoundaryMode

namespace MathEngine
{

// 3D move styles. "Diagonal" = body diagonals (8), "PlaneDiagonal" = face
// diagonals (12), "FullDiagonal" = both (20). Combined styles add the straight
// (6) set; "StraightFullDiagonal" is the full 26-neighbourhood.
enum class WalkerMoveStyle3D : std::uint8_t
{
    Straight                     = 0,  // 6
    PlaneDiagonal                = 1,  // 12
    Diagonal                     = 2,  // 8
    FullDiagonal                 = 3,  // 20
    StraightPlaneDiagonal        = 4,  // 18
    StraightDiagonal             = 5,  // 14
    StraightFullDiagonal         = 6,  // 26
    // ... with a "stay" option
    StraightWCenter              = 7,
    PlaneDiagonalWCenter         = 8,
    DiagonalWCenter              = 9,
    FullDiagonalWCenter          = 10,
    StraightPlaneDiagonalWCenter = 11,
    StraightDiagonalWCenter      = 12,
    StraightFullDiagonalWCenter  = 13,
    // ... continuous step magnitude in [0,1)
    StraightContinuous           = 14,
    PlaneDiagonalContinuous      = 15,
    DiagonalContinuous           = 16,
    FullDiagonalContinuous       = 17,
    StraightPlaneDiagonalContinuous = 18,
    StraightDiagonalContinuous   = 19,
    StraightFullDiagonalContinuous  = 20,
};

struct RandomWalk3DConfig
{
    size_t          numWalkers = 200;
    double          width      = 900.0;
    double          height     = 600.0;
    double          depth      = 600.0;
    double          size       = 5.0;
    double          startX     = 0.0;
    double          startY     = 0.0;
    double          startZ     = 0.0;
    double          stepSize   = 1.0;
    WalkerMoveStyle3D moveStyle = WalkerMoveStyle3D::Straight;
    BoundaryMode    boundary   = BoundaryMode::Free;
    std::uint64_t   seed       = 0;
};

class RandomWalk3D
{
public:
    // ---- Structure-of-Arrays state -----------------------------------------
    std::vector<double> posX, posY, posZ;    // current position (unwrapped)
    std::vector<double> initX, initY, initZ; // t = 0 reference
    std::vector<double> size;                // per-walker draw size

    // ---- geometry / bookkeeping -------------------------------------------
    double width  = 900.0;
    double height = 600.0;
    double depth  = 600.0;
    WalkerMoveStyle3D moveStyle = WalkerMoveStyle3D::Straight;
    BoundaryMode    boundary   = BoundaryMode::Free;
    double          stepSize   = 1.0;
    size_t          numWalkers = 0;
    std::uint64_t   stepCount  = 0;
    std::uint64_t   seed       = 0;

    std::mt19937_64 rng;

    // Direction set built from `moveStyle`.
    std::vector<std::array<double, 3>> dirs;
    bool continuous = false;

    RandomWalk3D() = default;

    explicit RandomWalk3D(const RandomWalk3DConfig& cfg)
        : width(cfg.width), height(cfg.height), depth(cfg.depth),
          moveStyle(cfg.moveStyle), boundary(cfg.boundary),
          stepSize(cfg.stepSize),
          numWalkers(cfg.numWalkers), seed(cfg.seed),
          rng(cfg.seed)
    {
        posX.assign(numWalkers, cfg.startX);
        posY.assign(numWalkers, cfg.startY);
        posZ.assign(numWalkers, cfg.startZ);
        initX.assign(numWalkers, cfg.startX);
        initY.assign(numWalkers, cfg.startY);
        initZ.assign(numWalkers, cfg.startZ);
        size.assign(numWalkers, cfg.size);
        BuildDirections(cfg.moveStyle, dirs, continuous);
    }

    void Step()
    {
        if (dirs.empty()) { ++stepCount; return; }
        std::uniform_int_distribution<int> uid(0, static_cast<int>(dirs.size()) - 1);
        std::uniform_real_distribution<double> umag(0.0, 1.0);

        for (size_t i = 0; i < numWalkers; ++i)
        {
            const std::array<double, 3>& d = dirs[uid(rng)];
            const double mag = continuous ? umag(rng) : 1.0;
            posX[i] += stepSize * mag * d[0];
            posY[i] += stepSize * mag * d[1];
            posZ[i] += stepSize * mag * d[2];

            if (boundary == BoundaryMode::Reflective)
            {
                posX[i] = Reflect(posX[i], 0.0, width);
                posY[i] = Reflect(posY[i], 0.0, height);
                posZ[i] = Reflect(posZ[i], 0.0, depth);
            }
        }
        ++stepCount;
    }

    double WrapX(double x) const { return x - width  * std::floor(x / width); }
    double WrapY(double y) const { return y - height * std::floor(y / height); }
    double WrapZ(double z) const { return z - depth  * std::floor(z / depth); }

    void Display(double& x, double& y, double& z, size_t i) const
    {
        x = posX[i]; y = posY[i]; z = posZ[i];
        if (boundary == BoundaryMode::Periodic)
        {
            x = WrapX(x); y = WrapY(y); z = WrapZ(z);
        }
    }

private:
    static double Reflect(double v, double lo, double hi)
    {
        const double span = hi - lo;
        if (span <= 0.0) return v;
        while (v < lo || v > hi)
        {
            if (v < lo) v = 2.0 * lo - v;
            if (v > hi) v = 2.0 * hi - v;
        }
        return v;
    }

    static void BuildDirections(WalkerMoveStyle3D ms, std::vector<std::array<double, 3>>& out, bool& cont)
    {
        out.clear();
        cont = false;
        bool s = false, p = false, d = false, c = false;

        switch (ms)
        {
            case WalkerMoveStyle3D::Straight:                     s = true; break;
            case WalkerMoveStyle3D::PlaneDiagonal:                p = true; break;
            case WalkerMoveStyle3D::Diagonal:                     d = true; break;
            case WalkerMoveStyle3D::FullDiagonal:                 p = d = true; break;
            case WalkerMoveStyle3D::StraightPlaneDiagonal:        s = p = true; break;
            case WalkerMoveStyle3D::StraightDiagonal:             s = d = true; break;
            case WalkerMoveStyle3D::StraightFullDiagonal:         s = p = d = true; break;
            case WalkerMoveStyle3D::StraightWCenter:              s = c = true; break;
            case WalkerMoveStyle3D::PlaneDiagonalWCenter:         p = c = true; break;
            case WalkerMoveStyle3D::DiagonalWCenter:              d = c = true; break;
            case WalkerMoveStyle3D::FullDiagonalWCenter:          p = d = c = true; break;
            case WalkerMoveStyle3D::StraightPlaneDiagonalWCenter: s = p = c = true; break;
            case WalkerMoveStyle3D::StraightDiagonalWCenter:      s = d = c = true; break;
            case WalkerMoveStyle3D::StraightFullDiagonalWCenter:  s = p = d = c = true; break;
            case WalkerMoveStyle3D::StraightContinuous:           s = true; cont = true; break;
            case WalkerMoveStyle3D::PlaneDiagonalContinuous:      p = true; cont = true; break;
            case WalkerMoveStyle3D::DiagonalContinuous:           d = true; cont = true; break;
            case WalkerMoveStyle3D::FullDiagonalContinuous:       p = d = true; cont = true; break;
            case WalkerMoveStyle3D::StraightPlaneDiagonalContinuous: s = p = true; cont = true; break;
            case WalkerMoveStyle3D::StraightDiagonalContinuous:   s = d = true; cont = true; break;
            case WalkerMoveStyle3D::StraightFullDiagonalContinuous: s = p = d = true; cont = true; break;
        }

        for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz)
        {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            const int nz = (dx != 0) + (dy != 0) + (dz != 0);
            bool keep = false;
            if (nz == 1 && s) keep = true;
            else if (nz == 2 && p) keep = true;
            else if (nz == 3 && d) keep = true;
            if (keep) out.push_back({{static_cast<double>(dx), static_cast<double>(dy), static_cast<double>(dz)}});
        }
        if (c) out.push_back({{0.0, 0.0, 0.0}});
    }
};

} // namespace MathEngine

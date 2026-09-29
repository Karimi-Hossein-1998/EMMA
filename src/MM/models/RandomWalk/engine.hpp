#pragma once
// -----------------------------------------------------------------------------
// Random-walk engine (Structure-of-Arrays, backend-independent).
//
// Owns all walker state in parallel arrays and is completely independent of any
// rendering backend — the direct analogue of `MolecularDynamics` for stochastic
// lattice/continuous walks. Positions are tracked *unwrapped* (raw accumulation
// of step increments); the on-canvas (wrapped/reflected) coordinate is produced
// on demand via `WrapX/Y` so that MSD, moments and diffusion stay meaningful.
// -----------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace MathEngine
{

// Discrete / continuous move styles.
enum class WalkerMoveStyle : std::uint8_t
{
    Straight                   = 0, // one of the 4 cardinal directions
    Diagonal                   = 1, // one of the 4 diagonal directions
    StraightDiagonal           = 2, // one of the 8 directions
    StraightWCenter            = 3, // cardinal + stay
    DiagonalWCenter            = 4, // diagonal + stay
    StraightDiagonalWCenter    = 5, // 8 directions + stay
    StraightContinuous         = 6, // cardinal axis, continuous step in [-1,1)
    DiagonalContinuous         = 7, // diagonal, continuous step in [-1,1)
    StraightDiagonalContinuous = 8, // independent continuous steps in x and y
};

// How walkers interact with the canvas edges.
enum class BoundaryMode : std::uint8_t
{
    Periodic   = 0, // wrap around (torus); displacement kept unwrapped
    Reflective = 1, // bounce off the walls
    Free       = 2, // no boundary; view auto-scales to fit
};

// All the knobs required to build a run (library-level analogue of RWParams).
struct RandomWalkConfig
{
    size_t          numWalkers = 200;
    double          width      = 900.0;
    double          height     = 600.0;
    double          size       = 5.0;
    double          startX     = 0.0;   // common starting position
    double          startY     = 0.0;
    double          stepSize   = 1.0;   // distance moved per step
    WalkerMoveStyle moveStyle  = WalkerMoveStyle::Straight;
    BoundaryMode    boundary   = BoundaryMode::Free;
    std::uint64_t   seed       = 0;
};

// Discrete direction tables (step increments).
static constexpr std::array<std::array<double, 2>, 4> kStraight{{
    {{1.0, 0.0}}, {{0.0, 1.0}}, {{-1.0, 0.0}}, {{0.0, -1.0}}
}};
static constexpr std::array<std::array<double, 2>, 4> kDiagonal{{
    {{1.0, 1.0}}, {{-1.0, 1.0}}, {{-1.0, -1.0}}, {{1.0, -1.0}}
}};
static constexpr std::array<std::array<double, 2>, 8> kStraightDiagonal{{
    {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}}, {{-1.0, 1.0}},
    {{-1.0, 0.0}}, {{-1.0, -1.0}}, {{0.0, -1.0}}, {{1.0, -1.0}}
}};
static constexpr std::array<std::array<double, 2>, 5> kStraightWCenter{{
    {{0.0, 0.0}}, {{1.0, 0.0}}, {{0.0, 1.0}}, {{-1.0, 0.0}}, {{0.0, -1.0}}
}};
static constexpr std::array<std::array<double, 2>, 5> kDiagonalWCenter{{
    {{0.0, 0.0}}, {{1.0, 1.0}}, {{-1.0, 1.0}}, {{-1.0, -1.0}}, {{1.0, -1.0}}
}};
static constexpr std::array<std::array<double, 2>, 9> kStraightDiagonalWCenter{{
    {{0.0, 0.0}}, {{1.0, 0.0}}, {{1.0, 1.0}}, {{0.0, 1.0}}, {{-1.0, 1.0}},
    {{-1.0, 0.0}}, {{-1.0, -1.0}}, {{0.0, -1.0}}, {{1.0, -1.0}}
}};

class RandomWalk
{
public:
    // ---- Structure-of-Arrays state -----------------------------------------
    std::vector<double> posX, posY;   // current position (unwrapped)
    std::vector<double> initX, initY; // t = 0 reference (for displacement/MSD)
    std::vector<double> size;         // per-walker draw size

    // ---- geometry / bookkeeping -------------------------------------------
    double          width      = 900.0;
    double          height     = 600.0;
    WalkerMoveStyle moveStyle  = WalkerMoveStyle::Straight;
    BoundaryMode    boundary   = BoundaryMode::Free;
    double          stepSize   = 1.0;
    size_t          numWalkers = 0;
    std::uint64_t   stepCount  = 0;
    std::uint64_t   seed       = 0;

    std::mt19937_64 rng;

    RandomWalk() = default;

    explicit RandomWalk(const RandomWalkConfig& cfg)
        : width(cfg.width), height(cfg.height),
          moveStyle(cfg.moveStyle), boundary(cfg.boundary),
          stepSize(cfg.stepSize),
          numWalkers(cfg.numWalkers), seed(cfg.seed),
          rng(cfg.seed)
    {
        posX.assign(numWalkers, cfg.startX);
        posY.assign(numWalkers, cfg.startY);
        initX.assign(numWalkers, cfg.startX);
        initY.assign(numWalkers, cfg.startY);
        size.assign(numWalkers, cfg.size);
    }

    // Advance every walker by exactly one step.
    void Step()
    {
        std::uniform_int_distribution<int> uid4(0, 3);
        std::uniform_int_distribution<int> uid5(0, 4);
        std::uniform_int_distribution<int> uid8(0, 7);
        std::uniform_int_distribution<int> uid9(0, 8);
        std::uniform_real_distribution<double> urd1(-1.0, 1.0);
        std::bernoulli_distribution bd(0.5);

        for (size_t i = 0; i < numWalkers; ++i)
        {
            double dx = 0.0, dy = 0.0;
            switch (moveStyle)
            {
                case WalkerMoveStyle::Straight: {
                    const int r = uid4(rng);
                    dx = kStraight[r][0]; dy = kStraight[r][1]; break;
                }
                case WalkerMoveStyle::Diagonal: {
                    const int r = uid4(rng);
                    dx = kDiagonal[r][0]; dy = kDiagonal[r][1]; break;
                }
                case WalkerMoveStyle::StraightDiagonal: {
                    const int r = uid8(rng);
                    dx = kStraightDiagonal[r][0]; dy = kStraightDiagonal[r][1]; break;
                }
                case WalkerMoveStyle::StraightWCenter: {
                    const int r = uid5(rng);
                    dx = kStraightWCenter[r][0]; dy = kStraightWCenter[r][1]; break;
                }
                case WalkerMoveStyle::DiagonalWCenter: {
                    const int r = uid5(rng);
                    dx = kDiagonalWCenter[r][0]; dy = kDiagonalWCenter[r][1]; break;
                }
                case WalkerMoveStyle::StraightDiagonalWCenter: {
                    const int r = uid9(rng);
                    dx = kStraightDiagonalWCenter[r][0]; dy = kStraightDiagonalWCenter[r][1]; break;
                }
                case WalkerMoveStyle::StraightContinuous: {
                    if (bd(rng)) dx = urd1(rng); else dy = urd1(rng); break;
                }
                case WalkerMoveStyle::DiagonalContinuous: {
                    const double s = urd1(rng);
                    dx = s; dy = (bd(rng) ? s : -s); break;
                }
                case WalkerMoveStyle::StraightDiagonalContinuous: {
                    dx = urd1(rng); dy = urd1(rng); break;
                }
                default: {
                    const int r = uid4(rng);
                    dx = kStraight[r][0]; dy = kStraight[r][1]; break;
                }
            }

            posX[i] += stepSize * dx;
            posY[i] += stepSize * dy;

            if (boundary == BoundaryMode::Reflective)
            {
                posX[i] = Reflect(posX[i], 0.0, width);
                posY[i] = Reflect(posY[i], 0.0, height);
            }
            // Periodic: leave unwrapped (wrap only for display).
            // Free: leave unbounded.
        }
        ++stepCount;
    }

    double WrapX(double x) const { return x - width  * std::floor(x / width); }
    double WrapY(double y) const { return y - height * std::floor(y / height); }

    // Display-space coordinate of walker i (respecting the boundary mode).
    void Display(double& x, double& y, size_t i) const
    {
        x = posX[i];
        y = posY[i];
        if (boundary == BoundaryMode::Periodic)
        {
            x = WrapX(x);
            y = WrapY(y);
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
};

} // namespace MathEngine

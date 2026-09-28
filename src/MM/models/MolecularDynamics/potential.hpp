#pragma once
// -----------------------------------------------------------------------------
// Pair potentials for the 2D molecular-dynamics model.
//
// Every potential is expressed through a single helper `PairForceEnergy` that
// returns, for a pair with squared distance r2:
//
//   fFactor : scalar such that the force ON particle i (the one "at the origin"
//             of the separation vector) is  (dx, dy) * fFactor.
//   pe      : the pair potential energy.
//
// All potentials use a SHIFTED-FORCE cutoff: the force is modified so it reaches
// exactly zero at r = r_c (and the potential is shifted to zero there too),
// which removes the force discontinuity of a plain truncation and greatly
// improves energy conservation in the NVE ensemble.
//
// The shift is linear in r^2 (i.e. f -> f - f_c * r / r_c), which is a standard
// shifted-force form and — unlike the constant-force shift — needs no sqrt in
// the inner loop, keeping Lennard-Jones/WCA branch-free and fast.
// -----------------------------------------------------------------------------
#include <cmath>

namespace MathEngine
{

// Softening added to r^2 so that perfectly overlapping particles (r == 0) do not
// produce a division by zero. It makes the short-range repulsion finite but very
// large, which is exactly the behaviour we want for an atomistically "hard" core.
inline constexpr double kMDEpsSqr = 1.0e-10;

enum class PotentialType
{
    LennardJones = 0,
    WCA,           // Weeks-Chandler-Andersen: repulsive core of LJ
    Morse,
    SoftSphere,    // power-law repulsion   eps*(sigma/r)^n
    Yukawa,        // screened Coulomb      eps*(sigma/r)*exp(-kappa r)
    Coulomb2D      // soft-core logarithmic  eps*ln((r+sigma)/r)
};

struct PotentialParams
{
    double sigma       = 1.0;   // length scale: LJ/WCA sigma, Morse equilibrium r0
    double epsilon     = 1.0;   // energy scale: LJ/WCA epsilon, Morse depth D
    double cutoffCoeff = 2.5;   // LJ/Morse cutoff = sigma * cutoffCoeff (WCA is fixed)
    double morseAlpha  = 1.0;   // Morse width (1/length); > 0
    double powerN      = 9.0;   // SoftSphere exponent n (>= 1)
    double yukawaKappa = 1.0;   // Yukawa inverse screening length (>= 0)

    // Derived quantities (set by FinalizePotential).
    double distCutOff    = 2.5;
    double distCutOffSqr = 6.25;
    double fCutoff       = 0.0;   // dU/dr evaluated at r_c
    double uCutoff       = 0.0;   // U(r_c)
};

// Computes the derived cutoff / shift constants for the given potential type.
inline void FinalizePotential(PotentialType type, PotentialParams& p)
{
    if (p.morseAlpha <= 0.0) p.morseAlpha = 1.0;
    if (p.powerN < 1.0) p.powerN = 1.0;
    if (p.yukawaKappa < 0.0) p.yukawaKappa = 0.0;

    switch (type)
    {
        case PotentialType::WCA:
            p.distCutOff    = std::pow(2.0, 1.0 / 6.0) * p.sigma;   // LJ minimum
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double s6  = std::pow(p.sigma / p.distCutOff, 6.0);
                const double s12 = s6 * s6;
                p.uCutoff = 4.0 * p.epsilon * (s12 - s6);
                p.fCutoff = 24.0 * p.epsilon / p.distCutOff * (s6 - 2.0 * s12);
            }
            break;

        case PotentialType::Morse:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double e = std::exp(-p.morseAlpha * (p.distCutOff - p.sigma));
                p.uCutoff = p.epsilon * (e * e - 2.0 * e);
                p.fCutoff = 2.0 * p.epsilon * p.morseAlpha * (1.0 - e) * e;
            }
            break;

        case PotentialType::SoftSphere:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double n = p.powerN;
                const double u = std::pow(p.sigma / p.distCutOff, n);   // (sigma/r_c)^n
                p.uCutoff = p.epsilon * u;
                p.fCutoff = -n * p.epsilon * u / p.distCutOff;          // dU/dr at r_c
            }
            break;

        case PotentialType::Yukawa:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double e = std::exp(-p.yukawaKappa * p.distCutOff);
                p.uCutoff = p.epsilon * (p.sigma / p.distCutOff) * e;
                p.fCutoff = -p.epsilon * p.sigma * e * (p.yukawaKappa / p.distCutOff
                                                       + 1.0 / (p.distCutOff * p.distCutOff));
            }
            break;

        case PotentialType::Coulomb2D:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                p.uCutoff = p.epsilon * std::log((p.distCutOff + p.sigma) / p.distCutOff);
                p.fCutoff = -p.epsilon * p.sigma / (p.distCutOff * (p.distCutOff + p.sigma));
            }
            break;

        case PotentialType::LennardJones:
        default:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double s6  = std::pow(p.sigma / p.distCutOff, 6.0);
                const double s12 = s6 * s6;
                p.uCutoff = 4.0 * p.epsilon * (s12 - s6);
                p.fCutoff = 24.0 * p.epsilon / p.distCutOff * (s6 - 2.0 * s12);
            }
            break;
    }
}

// Computes the pair interaction (shifted force + potential). Returns early with
// fFactor = pe = 0 when the pair is beyond the cutoff.
inline void PairForceEnergy(double r2, PotentialType type, const PotentialParams& p,
                            double& fFactor, double& pe)
{
    fFactor = 0.0;
    pe      = 0.0;
    if (r2 >= p.distCutOffSqr) return;

    const double rsafe = r2 + kMDEpsSqr;

    double f0, pe0;
    switch (type)
    {
        case PotentialType::Morse:
        {
            const double r = std::sqrt(r2);
            const double e = std::exp(-p.morseAlpha * (r - p.sigma));
            f0  = 2.0 * p.epsilon * p.morseAlpha * (1.0 - e) * e / r;
            pe0 = p.epsilon * (e * e - 2.0 * e);
            break;
        }
        case PotentialType::SoftSphere:
        {
            const double r = std::sqrt(r2);
            const double u = std::pow(p.sigma / r, p.powerN);   // (sigma/r)^n
            f0  = -p.powerN * p.epsilon * u / r2;               // (1/r) dU/dr
            pe0 = p.epsilon * u;
            break;
        }
        case PotentialType::Yukawa:
        {
            const double r = std::sqrt(r2);
            const double e = std::exp(-p.yukawaKappa * r);
            f0  = -p.epsilon * p.sigma * e * (p.yukawaKappa / r2 + 1.0 / (r2 * r));
            pe0 = p.epsilon * (p.sigma / r) * e;
            break;
        }
        case PotentialType::Coulomb2D:
        {
            const double r = std::sqrt(r2);
            f0  = -p.epsilon * p.sigma / (r2 * (r + p.sigma));
            pe0 = p.epsilon * std::log((r + p.sigma) / r);
            break;
        }
        case PotentialType::LennardJones:
        case PotentialType::WCA:
        default:
        {
            const double s2r2 = p.sigma * p.sigma / rsafe;   // (sigma/r)^2
            const double s6r6 = s2r2 * s2r2 * s2r2;          // (sigma/r)^6
            const double s12  = s6r6 * s6r6;                 // (sigma/r)^12
            f0  = 24.0 * p.epsilon / rsafe * (s6r6 - 2.0 * s12);
            pe0 = 4.0 * p.epsilon * (s12 - s6r6);
            break;
        }
    }

    // Linear shifted-force: f -> f - (f_c / r_c) * r,  U -> U - U_c - (f_c/r_c)(r^2 - r_c^2)/2.
    const double slope = p.fCutoff / p.distCutOff;
    fFactor = f0 - slope;
    pe      = pe0 - p.uCutoff - 0.5 * slope * (r2 - p.distCutOffSqr);
}

} // namespace MathEngine

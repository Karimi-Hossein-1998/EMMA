#pragma once
// -----------------------------------------------------------------------------
// Pair potentials for the 3D molecular-dynamics model.
//
// Same radial convention as the 2D version: `PairForceEnergy3D` returns
//   fFactor : scalar such that the force ON particle i is (dx,dy,dz) * fFactor
//   pe      : pair potential energy
// All potentials use a linear-in-r^2 shifted-force cutoff. The only difference
// from 2D is Coulomb, which is 1/r in 3D (instead of logarithmic).
// -----------------------------------------------------------------------------
#include <cmath>

namespace MathEngine
{

inline constexpr double kMDEpsSqr3D = 1.0e-10;

enum class PotentialType3D
{
    LennardJones = 0,
    WCA,          // Weeks-Chandler-Andersen: repulsive core of LJ
    Morse,
    SoftSphere,   // power-law repulsion  eps*(sigma/r)^n
    Yukawa,       // screened Coulomb     eps*(sigma/r)*exp(-kappa r)
    Coulomb3D     // soft-core 1/r        eps*sigma/r
};

struct PotentialParams3D
{
    double sigma       = 1.0;
    double epsilon     = 1.0;
    double cutoffCoeff = 2.5;
    double morseAlpha  = 1.0;
    double powerN      = 9.0;
    double yukawaKappa = 1.0;

    double distCutOff    = 2.5;
    double distCutOffSqr = 6.25;
    double fCutoff       = 0.0;   // dU/dr at r_c
    double uCutoff       = 0.0;   // U(r_c)
};

inline void FinalizePotential3D(PotentialType3D type, PotentialParams3D& p)
{
    if (p.morseAlpha <= 0.0) p.morseAlpha = 1.0;
    if (p.powerN < 1.0) p.powerN = 1.0;
    if (p.yukawaKappa < 0.0) p.yukawaKappa = 0.0;

    switch (type)
    {
        case PotentialType3D::WCA:
            p.distCutOff    = std::pow(2.0, 1.0 / 6.0) * p.sigma;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double s6  = std::pow(p.sigma / p.distCutOff, 6.0);
                const double s12 = s6 * s6;
                p.uCutoff = 4.0 * p.epsilon * (s12 - s6);
                p.fCutoff = 24.0 * p.epsilon / p.distCutOff * (s6 - 2.0 * s12);
            }
            break;
        case PotentialType3D::Morse:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double e = std::exp(-p.morseAlpha * (p.distCutOff - p.sigma));
                p.uCutoff = p.epsilon * (e * e - 2.0 * e);
                p.fCutoff = 2.0 * p.epsilon * p.morseAlpha * (1.0 - e) * e;
            }
            break;
        case PotentialType3D::SoftSphere:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double n = p.powerN;
                const double u = std::pow(p.sigma / p.distCutOff, n);
                p.uCutoff = p.epsilon * u;
                p.fCutoff = -n * p.epsilon * u / p.distCutOff;
            }
            break;
        case PotentialType3D::Yukawa:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                const double e = std::exp(-p.yukawaKappa * p.distCutOff);
                p.uCutoff = p.epsilon * (p.sigma / p.distCutOff) * e;
                p.fCutoff = -p.epsilon * p.sigma * e * (p.yukawaKappa / p.distCutOff + 1.0 / (p.distCutOff * p.distCutOff));
            }
            break;
        case PotentialType3D::Coulomb3D:
            p.distCutOff    = p.sigma * p.cutoffCoeff;
            p.distCutOffSqr = p.distCutOff * p.distCutOff;
            {
                p.uCutoff = p.epsilon * p.sigma / p.distCutOff;
                p.fCutoff = -p.epsilon * p.sigma / (p.distCutOff * p.distCutOff);
            }
            break;
        case PotentialType3D::LennardJones:
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

inline void PairForceEnergy3D(double r2, PotentialType3D type, const PotentialParams3D& p,
                              double& fFactor, double& pe)
{
    fFactor = 0.0;
    pe      = 0.0;
    if (r2 >= p.distCutOffSqr) return;

    const double rsafe = r2 + kMDEpsSqr3D;

    double f0, pe0;
    switch (type)
    {
        case PotentialType3D::Morse:
        {
            const double r = std::sqrt(rsafe);
            const double e = std::exp(-p.morseAlpha * (r - p.sigma));
            f0  = 2.0 * p.epsilon * p.morseAlpha * (1.0 - e) * e / r;
            pe0 = p.epsilon * (e * e - 2.0 * e);
            break;
        }
        case PotentialType3D::SoftSphere:
        {
            const double r = std::sqrt(rsafe);
            const double u = std::pow(p.sigma / r, p.powerN);
            f0  = -p.powerN * p.epsilon * u / rsafe;
            pe0 = p.epsilon * u;
            break;
        }
        case PotentialType3D::Yukawa:
        {
            const double r = std::sqrt(rsafe);
            const double e = std::exp(-p.yukawaKappa * r);
            f0  = -p.epsilon * p.sigma * e * (p.yukawaKappa / rsafe + 1.0 / (rsafe * r));
            pe0 = p.epsilon * (p.sigma / r) * e;
            break;
        }
        case PotentialType3D::Coulomb3D:
        {
            const double r = std::sqrt(rsafe);
            f0  = -p.epsilon * p.sigma / (rsafe * r);
            pe0 = p.epsilon * p.sigma / r;
            break;
        }
        case PotentialType3D::LennardJones:
        case PotentialType3D::WCA:
        default:
        {
            const double s2r2 = p.sigma * p.sigma / rsafe;
            const double s6r6 = s2r2 * s2r2 * s2r2;
            const double s12  = s6r6 * s6r6;
            f0  = 24.0 * p.epsilon / rsafe * (s6r6 - 2.0 * s12);
            pe0 = 4.0 * p.epsilon * (s12 - s6r6);
            break;
        }
    }

    const double slope = p.fCutoff / p.distCutOff;
    fFactor = f0 - slope;
    pe      = pe0 - p.uCutoff - 0.5 * slope * (r2 - p.distCutOffSqr);
}

} // namespace MathEngine

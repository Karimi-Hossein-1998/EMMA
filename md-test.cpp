// -----------------------------------------------------------------------------
// Standalone check for the molecular-dynamics model (MathEngine port).
//
// Mirrors the Nature-of-Code validation harness, but runs inside the EMMA
// `MathEngine` namespace and dumps results through MathEngine::IO.
//
// Build (from the repo root):
//   g++ -std=c++23 -O2 -I. md-test.cpp -o md-test
// -----------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <vector>

#include "src/MM/typedefs/header.hpp"
#include "src/MM/utility/write.hpp"
#include "src/MM/models/molecular-dynamics.hpp"

using namespace MathEngine;

static void printPotentialSanity()
{
    std::printf("== 1. pair-potential sanity ==\n");
    const PotentialType types[] = { PotentialType::LennardJones, PotentialType::WCA, PotentialType::Morse };
    const char* names[] = { "LennardJones", "WCA", "Morse" };

    for (int t = 0; t < 3; ++t)
    {
        PotentialParams p;
        p.sigma = 1.0; p.epsilon = 1.0; p.cutoffCoeff = 2.5; p.morseAlpha = 1.0;
        FinalizePotential(types[t], p);

        double fFactor, pe;
        PairForceEnergy(0.81, types[t], p, fFactor, pe);
        const bool repulsive = fFactor < 0.0;

        double fFar, peFar;
        PairForceEnergy(9.0, types[t], p, fFar, peFar);

        std::printf("  %-14s r=0.9: fFactor=%+.4f (%s)   r=3.0: f=%.4f pe=%.4f (%s)\n",
                    names[t], fFactor, repulsive ? "repulsive OK" : "WRONG",
                    fFar, peFar, (fFar == 0.0 && peFar == 0.0) ? "cutoff OK" : "cutoff WRONG");
    }
    std::printf("\n");
}

int main()
{
    printPotentialSanity();

    // ---- 2. NVE energy conservation ----------------------------------------
    {
        std::printf("== 2. NVE energy conservation (velocity-Verlet, dt=0.001) ==\n");
        MDConfig cfg;
        cfg.numParticles = 64;
        cfg.width = cfg.height = 11.22;
        cfg.sigma = 1.0; cfg.epsilon = 1.0; cfg.cutoffCoeff = 2.5;
        cfg.temperature = 0.5; cfg.seed = 41;
        cfg.periodicBoundaryCondition = true; cfg.bounce = false;
        cfg.initialCondition = InitialConditionType::SquareLattice;
        cfg.potential = PotentialType::LennardJones;

        MolecularDynamics md(cfg);
        md.ComputeKineticEnergy(); md.ComputePotentialEnergy();
        const double E0 = md.kineticEnergy + md.potentialEnergy;
        double Emin = E0, Emax = E0;

        const double dt = 0.001;
        for (int s = 0; s < 5000; ++s)
        {
            md.Step(dt);
            if (s % 10 == 0)
            {
                md.ComputeKineticEnergy(); md.ComputePotentialEnergy();
                const double E = md.kineticEnergy + md.potentialEnergy;
                if (E < Emin) Emin = E;
                if (E > Emax) Emax = E;
            }
        }
        std::printf("  E0=%.6f  E in [%.6f, %.6f]  relative drift=%.3e  %s\n\n",
                    E0, Emin, Emax, (Emax - Emin) / std::abs(E0),
                    (Emax - Emin) / std::abs(E0) < 5e-3 ? "(OK)" : "(large)");
    }

    // ---- 3. crystallisation (psi6 rises) + collect time series -------------
    std::vector<double> tHist, THist, KEHist, PEHist, psi6Hist, psi4Hist, msdHist;
    {
        std::printf("== 3. Berendsen cooling -> psi6 ordering ==\n");
        MDConfig cfg;
        cfg.numParticles = 100;
        cfg.width = cfg.height = 11.22;
        cfg.sigma = 1.0; cfg.epsilon = 1.0; cfg.cutoffCoeff = 2.5;
        cfg.temperature = 1.5; cfg.seed = 7;
        cfg.periodicBoundaryCondition = true; cfg.bounce = false;
        cfg.initialCondition = InitialConditionType::SquareLattice;
        cfg.potential = PotentialType::LennardJones;

        MolecularDynamics md(cfg);
        const OrderParams start = ComputeBondOrientationalOrder(md, 1.4);

        const double dt = 0.001;
        for (int s = 0; s < 20000; ++s)
        {
            md.Step(dt);
            ApplyBerendsen(md, 0.10, dt, 2.0);
            if (s % 1000 == 0)
            {
                const Observables o = CollectObservables(md, s * dt, 1.4);
                tHist.push_back(o.time);    THist.push_back(o.temperature);
                KEHist.push_back(o.kineticEnergy); PEHist.push_back(o.potentialEnergy);
                psi6Hist.push_back(o.psi6); psi4Hist.push_back(o.psi4);
                msdHist.push_back(o.msd);
            }
        }
        const Observables o = CollectObservables(md, 20000 * dt, 1.4);
        std::printf("  psi6: start=%.4f -> end=%.4f   T: 1.5 -> %.4f  %s\n\n",
                    start.psi6, o.psi6, o.temperature,
                    o.psi6 > start.psi6 + 0.3 ? "(ordered, OK)" : "(check)");
    }

    // ---- 4. thermostat control ---------------------------------------------
    {
        std::printf("== 4. thermostat control (target T = 0.8) ==\n");
        const char* names[] = { "Rescale", "Berendsen", "Andersen", "Langevin" };
        for (int m = 0; m < 4; ++m)
        {
            MDConfig cfg;
            cfg.numParticles = 64;
            cfg.width = cfg.height = 11.22;
            cfg.sigma = 1.0; cfg.epsilon = 1.0; cfg.cutoffCoeff = 2.5;
            cfg.temperature = 2.0; cfg.seed = 41;
            cfg.periodicBoundaryCondition = true; cfg.bounce = false;
            cfg.initialCondition = InitialConditionType::SquareLattice;
            cfg.potential = PotentialType::LennardJones;

            MolecularDynamics md(cfg);
            size_t rngSeed = 123;
            const double dt = 0.001;
            double Tsum = 0.0; int Tcount = 0;
            for (int s = 0; s < 8000; ++s)
            {
                md.Step(dt);
                switch (m)
                {
                    case 0: ApplyVelocityRescale(md, 0.8); break;
                    case 1: ApplyBerendsen(md, 0.8, dt, 1.0); break;
                    case 2: ApplyAndersen(md, 0.8, dt, 5.0, rngSeed); break;
                    case 3: ApplyLangevin(md, 0.8, dt, 1.0, rngSeed); break;
                }
                if (s >= 4000) { md.ComputeKineticEnergy(); Tsum += md.currentTemperature; ++Tcount; }
            }
            std::printf("  %-10s <T> = %.4f (target 0.8)\n", names[m], Tsum / Tcount);
        }
        std::printf("\n");
    }

    // ---- 5. dump the cooling time series via MathEngine::IO ----------------
    {
        const size_t rows = tHist.size();
        const size_t cols = 7;
        dMatrix m(rows, cols, 0.0);
        for (size_t r = 0; r < rows; ++r)
        {
            m[r, 0] = tHist[r];   m[r, 1] = THist[r];  m[r, 2] = KEHist[r];
            m[r, 3] = PEHist[r];  m[r, 4] = psi6Hist[r];
            m[r, 5] = psi4Hist[r]; m[r, 6] = msdHist[r];
        }

        IO::WriteOptions wo;
        wo.path      = "EMMAOutput/md-test/observables.csv";
        wo.separator = ",";
        wo.header    = "time, temperature, kineticEnergy, potentialEnergy, psi6, psi4, msd";
        wo.comment   = {};
        wo.footer    = {};
        wo.colWidth  = 20;
        wo.precision = 15;
        wo.format    = IO::FPFormat::Scientific;
        wo.alignment = IO::Alignment::Center;
        wo.append    = false;
        wo.binary    = false;
        IO::WriteMatrix(m, wo);
        std::printf("Observables written to EMMAOutput/md-test/observables.csv\n");
    }

    std::printf("Done.\n");
    return 0;
}

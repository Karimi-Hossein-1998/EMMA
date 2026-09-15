#pragma once
#include "../../typedefs/header.hpp"

// SVD-based solver for Ax = b using Jacobi SVD (for small matrices, demonstration only).
// A: n x n matrix (dMatrix), b: vector of size n (dVec).
// Returns the solution vector x (dVec).
inline dVec svd_solve(dMatrix A, dVec b)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("dMatrix A must be square and match the size of b.");

    // Initialize U, S, V^T
    dMatrix U = A;
    dMatrix V(n, dVec(n, 0.0));
    dVec S(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        V[i][i] = 1.0;

    // Jacobi SVD (very basic, not optimized, for demonstration)
    for (size_t sweep = 0; sweep < 30; ++sweep)
    {
        for (size_t p = 0; p < n - 1; ++p)
        {
            for (size_t q = p + 1; q < n; ++q)
            {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (size_t k = 0; k < n; ++k)
                {
                    alpha += U[k][p] * U[k][p];
                    beta  += U[k][q] * U[k][q];
                    gamma += U[k][p] * U[k][q];
                }
                if (gamma == 0.0) continue;
                double zeta = (beta - alpha) / (2.0 * gamma);
                double t = ((zeta >= 0.0) ? 1.0 : -1.0) / (std::fabs(zeta) + std::sqrt(1.0 + zeta * zeta));
                double c = 1.0 / std::sqrt(1.0 + t * t);
                double s = c * t;
                for (size_t k = 0; k < n; ++k)
                {
                    double u_p = U[k][p];
                    double u_q = U[k][q];
                    U[k][p] = c * u_p - s * u_q;
                    U[k][q] = s * u_p + c * u_q;
                }
                for (size_t k = 0; k < n; ++k)
                {
                    double v_p = V[k][p];
                    double v_q = V[k][q];
                    V[k][p] = c * v_p - s * v_q;
                    V[k][q] = s * v_p + c * v_q;
                }
            }
        }
    }
    // Compute singular values
    for (size_t i = 0; i < n; ++i)
    {
        double norm = 0.0;
        for (size_t k = 0; k < n; ++k)
            norm += U[k][i] * U[k][i];
        S[i] = std::sqrt(norm);
        if (S[i] > 1e-12)
            for (size_t k = 0; k < n; ++k)
                U[k][i] /= S[i];
    }
    // Compute U^T b
    dVec Utb(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t k = 0; k < n; ++k)
            Utb[i] += U[k][i] * b[k];
    // Solve for y: S y = U^T b
    dVec y(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        y[i] = (S[i] > 1e-12) ? Utb[i] / S[i] : 0.0;
    // Compute x = V y
    dVec x(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t k = 0; k < n; ++k)
            x[i] += V[i][k] * y[k];
    return x;
} 
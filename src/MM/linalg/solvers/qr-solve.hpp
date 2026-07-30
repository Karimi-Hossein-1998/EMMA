#pragma once
#include "../../typedefs/header.hpp"

// QR solver for general A: solves A x = b using Gram-Schmidt QR decomposition.
// A: n x n matrix (dMatrix), b: vector of size n (dVec).
// Returns the solution vector x (dVec).
inline dVec qr_solve(dMatrix A, dVec b)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("dMatrix A must be square and match the size of b.");

    // Gram-Schmidt process to compute Q and R
    dMatrix Q(n, dVec(n, 0.0));
    dMatrix R(n, dVec(n, 0.0));
    for (size_t k = 0; k < n; ++k)
    {
        // Copy k-th column of A to q
        dVec q = A[k];
        for (size_t i = 0; i < k; ++i)
        {
            R[i][k] = 0.0;
            for (size_t j = 0; j < n; ++j)
                R[i][k] += Q[i][j] * A[k][j];
            for (size_t j = 0; j < n; ++j)
                q[j] -= R[i][k] * Q[i][j];
        }
        // Compute norm
        double norm = 0.0;
        for (size_t j = 0; j < n; ++j)
            norm += q[j] * q[j];
        norm = std::sqrt(norm);
        if (norm == 0.0)
            throw std::runtime_error("dMatrix is rank deficient.");
        for (size_t j = 0; j < n; ++j)
            Q[k][j] = q[j] / norm;
        R[k][k] = norm;
    }

    // Compute Q^T b
    dVec Qtb(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            Qtb[i] += Q[i][j] * b[j];

    // Backward substitution to solve Rx = Q^T b
    dVec x(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        x[i] = Qtb[i];
        for (size_t j = i + 1; j < n; ++j)
            x[i] -= R[i][j] * x[j];
        x[i] /= R[i][i];
    }
    return x;
} 
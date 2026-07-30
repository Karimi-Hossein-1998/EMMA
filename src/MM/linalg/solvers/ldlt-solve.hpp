#pragma once
#include "../../typedefs/header.hpp"

// LDL^T solver for symmetric A: solves A x = b using LDL^T decomposition.
// A: n x n matrix (dMatrix), b: vector of size n (dVec).
// Returns the solution vector x (dVec).
inline dVec ldlt_solve(dMatrix A, dVec b)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("dMatrix A must be square and match the size of b.");

    // LDL^T decomposition: A = L D L^T
    dMatrix L(n, dVec(n, 0.0));
    dVec D(n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            double sum = A[i][j];
            for (size_t k = 0; k < j; ++k)
                sum -= L[i][k] * D[k] * L[j][k];
            L[i][j] = sum / D[j];
        }
        double sum = A[i][i];
        for (size_t k = 0; k < i; ++k)
            sum -= L[i][k] * L[i][k] * D[k];
        D[i] = sum;
        L[i][i] = 1.0;
    }

    // Forward substitution: L y = b
    dVec y(n);
    for (size_t i = 0; i < n; ++i)
    {
        y[i] = b[i];
        for (size_t j = 0; j < i; ++j)
            y[i] -= L[i][j] * y[j];
    }

    // Diagonal solve: D z = y
    dVec z(n);
    for (size_t i = 0; i < n; ++i)
        z[i] = y[i] / D[i];

    // Backward substitution: L^T x = z
    dVec x(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        x[i] = z[i];
        for (size_t j = i + 1; j < n; ++j)
            x[i] -= L[j][i] * x[j];
    }
    return x;
} 
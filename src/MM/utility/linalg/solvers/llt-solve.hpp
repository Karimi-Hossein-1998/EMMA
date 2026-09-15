#pragma once
#include "../../typedefs/header.hpp"

// Cholesky (LL^T) solver for symmetric positive definite A: solves A x = b.
// A: n x n matrix (dMatrix), b: vector of size n (dVec).
// Returns the solution vector x (dVec).
inline dVec llt_solve(dMatrix A, dVec b)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("dMatrix A must be square and match the size of b.");

    // Cholesky decomposition: A = L L^T
    dMatrix L(n, dVec(n, 0.0));
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double sum = A[i][j];
            for (size_t k = 0; k < j; ++k)
                sum -= L[i][k] * L[j][k];
            if (i == j)
            {
                if (sum <= 0.0)
                    throw std::runtime_error("dMatrix is not positive definite.");
                L[i][j] = std::sqrt(sum);
            }
            else
            {
                L[i][j] = sum / L[j][j];
            }
        }
    }

    // Forward substitution: L y = b
    dVec y(n);
    for (size_t i = 0; i < n; ++i)
    {
        y[i] = b[i];
        for (size_t j = 0; j < i; ++j)
            y[i] -= L[i][j] * y[j];
        y[i] /= L[i][i];
    }

    // Backward substitution: L^T x = y
    dVec x(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        x[i] = y[i];
        for (size_t j = i + 1; j < n; ++j)
            x[i] -= L[j][i] * x[j];
        x[i] /= L[i][i];
    }
    return x;
} 
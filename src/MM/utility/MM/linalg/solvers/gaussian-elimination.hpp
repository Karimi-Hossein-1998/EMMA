#pragma once
#include "../../typedefs/header.hpp"

// Gaussian elimination with partial pivoting for solving A x = b.
// A: n x n matrix (dMatrix), b: vector of size n (dVec).
// Returns the solution vector x (dVec).
inline dVec gaussian_elimination(dMatrix A, dVec b)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("dMatrix A must be square and match the size of b.");

    // Forward elimination
    for (size_t k = 0; k < n; ++k)
    {
        // Partial pivoting
        size_t max_row = k;
        for (size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[max_row][k]))
                max_row = i;
        if (A[max_row][k] == 0.0)
            throw std::runtime_error("dMatrix is singular or nearly singular.");

        // Swap rows in A and b
        std::swap(A[k], A[max_row]);
        std::swap(b[k], b[max_row]);

        // Eliminate below
        for (size_t i = k + 1; i < n; ++i)
        {
            double factor = A[i][k] / A[k][k];
            for (size_t j = k; j < n; ++j)
                A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    dVec x(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        double sum = b[i];
        for (size_t j = i + 1; j < n; ++j)
            sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }
    return x;
} 

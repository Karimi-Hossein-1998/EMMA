#pragma once
#include "../../typedefs/header.hpp"

// LU factorization solver with partial pivoting for solving A x = b.
// A: n x n matrix (dMatrix), b: vector of size n (dVec).
// Returns the solution vector x (dVec).
inline dVec lu_solve(dMatrix A, dVec b)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("dMatrix A must be square and match the size of b.");

    // Permutation vector for partial pivoting
    std::vector<size_t> P(n);
    for (size_t i = 0; i < n; ++i) P[i] = i;

    // LU Decomposition with partial pivoting
    for (size_t k = 0; k < n; ++k)
    {
        // Find pivot
        size_t max_row = k;
        for (size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[max_row][k]))
                max_row = i;
        if (A[max_row][k] == 0.0)
            throw std::runtime_error("dMatrix is singular or nearly singular.");
        // Swap rows in A and permutation vector
        std::swap(A[k], A[max_row]);
        std::swap(P[k], P[max_row]);
        // Elimination
        for (size_t i = k + 1; i < n; ++i)
        {
            A[i][k] /= A[k][k];
            for (size_t j = k + 1; j < n; ++j)
                A[i][j] -= A[i][k] * A[k][j];
        }
    }

    // Apply permutation to b
    dVec bp(n);
    for (size_t i = 0; i < n; ++i)
        bp[i] = b[P[i]];

    // Forward substitution to solve Ly = Pb
    dVec y(n);
    for (size_t i = 0; i < n; ++i)
    {
        y[i] = bp[i];
        for (size_t j = 0; j < i; ++j)
            y[i] -= A[i][j] * y[j];
    }

    // Backward substitution to solve Ux = y
    dVec x(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        x[i] = y[i];
        for (size_t j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
} 
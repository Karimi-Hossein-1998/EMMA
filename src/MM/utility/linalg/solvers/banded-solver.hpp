#pragma once
#include "../../typedefs/header.hpp"

// Banded matrix solver for Ax = b using compact banded storage and Gaussian elimination.
// ab: banded matrix storage (2*kl+ku+1) x n, kl: lower bandwidth, ku: upper bandwidth, b: RHS (n)
// Returns the solution vector x (dVec).
// ab is stored so that ab[kl+ku+i-j][j] == A[i][j] for max(0, j-ku) <= i <= min(n-1, j+kl)
inline dVec banded_solve(dMatrix ab, int kl, int ku, dVec b)
{
    size_t n = b.size();
    if (ab.size() != static_cast<size_t>(kl+ku+1) || ab[0].size() != n || b.size() != n)
        throw std::invalid_argument("Invalid banded matrix or vector size.");
    // Forward elimination
    for (size_t k = 0; k < n; ++k)
    {
        // Find pivot (no pivoting for simplicity)
        double pivot = ab[kl+ku][k];
        if (pivot == 0.0)
            throw std::runtime_error("Zero pivot encountered in banded solver.");
        // Eliminate below
        size_t i_max = std::min(n-1, k+kl);
        for (size_t i = k+1; i <= i_max; ++i)
        {
            int band = kl+ku+i-k;
            double factor = ab[band][k] / pivot;
            ab[band][k] = 0.0;
            // Update the rest of the row
            for (size_t j = k+1; j <= std::min(n-1, k+ku); ++j)
                ab[band][j] -= factor * ab[kl+ku+j-k][k];
            b[i] -= factor * b[k];
        }
    }
    // Back substitution
    dVec x(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        double sum = b[i];
        for (size_t j = i+1; j <= std::min(n-1, static_cast<size_t>(i+ku)); ++j)
            sum -= ab[kl+ku+j-i][i] * x[j];
        x[i] = sum / ab[kl+ku][i];
    }
    return x;
} 
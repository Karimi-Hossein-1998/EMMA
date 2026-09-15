#pragma once
#include "../../typedefs/header.hpp"

// Tridiagonal solver (Thomas algorithm) for Ax = d, where A is tridiagonal.
// a: subdiagonal (n-1), b: diagonal (n), c: superdiagonal (n-1), d: right-hand side (n)
// Returns the solution vector x (dVec).
// For large sparse tridiagonal systems, consider using SparsedMatrix, but dVecs are most efficient for classic Thomas algorithm.
inline dVec tridiagonal_solve(const dVec& a, const dVec& b, const dVec& c, const dVec& d)
{
    size_t n = b.size();
    if (a.size() != n-1 || c.size() != n-1 || d.size() != n)
        throw std::invalid_argument("Invalid tridiagonal vector sizes.");
    dVec cp = c, dp = d, bp = b;
    // Forward sweep
    for (size_t i = 1; i < n; ++i)
    {
        double m = a[i-1] / bp[i-1];
        bp[i] -= m * cp[i-1];
        dp[i] -= m * dp[i-1];
    }
    // Back substitution
    dVec x(n);
    x[n-1] = dp[n-1] / bp[n-1];
    for (int i = static_cast<int>(n) - 2; i >= 0; --i)
        x[i] = (dp[i] - cp[i] * x[i+1]) / bp[i];
    return x;
} 
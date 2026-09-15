#pragma once
#include "../../typedefs/header.hpp"
#include "decomposers.hpp"
// Returns the trace of a square matrix (sum of diagonal elements)
inline double trace(const dMatrix& A)
{
    const size_t n = A.size();
    double tr = 0.0;
    for (size_t i = 0; i < n; ++i)
        tr += A[i][i];
    return tr;
}

// Optimized determinant using LU decomposition for n > 3
inline double determinant(const dMatrix& A)
{
    const size_t n = A.size();
    if (n == 1) return A[0][0];
    if (n == 2) return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (n == 3) 
    {
        return A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
             - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
             + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    }
    // For larger matrices, use LU decomposition
    auto [L, U, P] = lu_decomposition(A);
    double det = 1.0;
    int parity = 1;
    for (size_t i = 0; i < n; ++i)
    {
        det *= U[i][i];
        if (P[i] != i) parity *= -1;
    }
    return det * parity;
}

// Optimized minor: only use Laplace expansion for n <= 3
inline double minor(const dMatrix& A, size_t row, size_t col)
{
    const size_t n = A.size();
    dMatrix subA;
    for (size_t i = 0; i < n; ++i)
    {
        if (i == row) continue;
        dVec subrow;
        for (size_t j = 0; j < n; ++j)
        {
            if (j == col) continue;
            subrow.push_back(A[i][j]);
        }
        subA.push_back(subrow);
    }
    if (subA.size() == 1) return subA[0][0];
    if (subA.size() == 2) return subA[0][0]*subA[1][1] - subA[0][1]*subA[1][0];
    if (subA.size() == 3)
    {
        return subA[0][0]*(subA[1][1]*subA[2][2] - subA[1][2]*subA[2][1])
             - subA[0][1]*(subA[1][0]*subA[2][2] - subA[1][2]*subA[2][0])
             + subA[0][2]*(subA[1][0]*subA[2][1] - subA[1][1]*subA[2][0]);
    }
    // For larger matrices, use LU-based determinant
    return determinant(subA);
}

// Returns the transpose of a matrix
inline dMatrix transpose(const dMatrix& A)
{
    const size_t n = A.size();
    if (n == 0) return dMatrix();
    const size_t m = A[0].size();
    dMatrix At(m, dVec(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            At[j][i] = A[i][j];
    return At;
}

// Optimized cofactor_matrix and adjugate (no change needed, but now use optimized minor)
inline dMatrix cofactor_matrix(const dMatrix& A)
{
    const size_t n = A.size();
    dMatrix C(n, dVec(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            C[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * minor(A, i, j);
    return C;
}

// Returns the adjugate (classical adjoint) of A
inline dMatrix adjugate(const dMatrix& A)
{
    return transpose(cofactor_matrix(A));
}

// Optimized inverse: use LU decomposition for n > 3
inline dMatrix inverse(const dMatrix& A)
{
    const size_t n = A.size();
    if (n <= 3)
    {
        double det = determinant(A);
        if (std::abs(det) < 1e-12)
            throw std::runtime_error("dMatrix is singular or nearly singular.");
        dMatrix adj = adjugate(A);
        dMatrix inv(n, dVec(n));
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                inv[i][j] = adj[i][j] / det;
        return inv;
    }
    // For larger matrices, use LU decomposition
    auto [L, U, P] = lu_decomposition(A);
    dMatrix inv(n, dVec(n));
    // Forward and backward substitution for each column of identity
    for (size_t col = 0; col < n; ++col)
    {
        dVec e(n, 0.0); e[col] = 1.0;
        // Apply permutation
        dVec b(n);
        for (size_t i = 0; i < n; ++i) b[i] = e[P[i]];
        // Forward substitution (Ly = Pb)
        dVec y(n);
        for (size_t i = 0; i < n; ++i)
        {
            y[i] = b[i];
            for (size_t j = 0; j < i; ++j)
                y[i] -= L[i][j] * y[j];
        }
        // Backward substitution (Ux = y)
        dVec x(n);
        for (int i = static_cast<int>(n) - 1; i >= 0; --i)
        {
            x[i] = y[i];
            for (size_t j = i + 1; j < n; ++j)
                x[i] -= U[i][j] * x[j];
            x[i] /= U[i][i];
        }
        for (size_t i = 0; i < n; ++i)
            inv[i][col] = x[i];
    }
    return inv;
}

// Returns true if A is symmetric
inline bool is_symmetric(const dMatrix& A)
{
    const size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (A[i][j] != A[j][i])
                return false;
    return true;
}

// Returns true if A is the identity matrix
inline bool is_identity(const dMatrix& A)
{
    const size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if ((i == j && A[i][j] != 1.0) || (i != j && A[i][j] != 0.0))
                return false;
    return true;
}

// Returns true if A is diagonal
inline bool is_diagonal(const dMatrix& A)
{
    const size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (i != j && A[i][j] != 0.0)
                return false;
    return true;
}

// Returns true if A is upper triangular
inline bool is_upper_triangular(const dMatrix& A)
{
    const size_t n = A.size();
    for (size_t i = 1; i < n; ++i)
        for (size_t j = 0; j < i; ++j)
            if (A[i][j] != 0.0)
                return false;
    return true;
}

// Returns true if A is lower triangular
inline bool is_lower_triangular(const dMatrix& A)
{
    const size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            if (A[i][j] != 0.0)
                return false;
    return true;
}

// Returns the Frobenius norm of a matrix
inline double frobenius_norm(const dMatrix& A)
{
    double sum = 0.0;
    const size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < A[i].size(); ++j)
            sum += A[i][j] * A[i][j];
    return std::sqrt(sum);
}

// Returns the max absolute row sum norm (infinity norm)
inline double max_norm(const dMatrix& A)
{
    double max_sum = 0.0;
    const size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
    {
        double row_sum = 0.0;
        for (size_t j = 0; j < A[i].size(); ++j)
            row_sum += std::fabs(A[i][j]);
        if (row_sum > max_sum)
            max_sum = row_sum;
    }
    return max_sum;
}

// Returns the max absolute column sum norm (1-norm)
inline double one_norm(const dMatrix& A)
{
    if (A.empty()) return 0.0;
    const size_t n = A.size(), m = A[0].size();
    double max_sum = 0.0;
    for (size_t j = 0; j < m; ++j)
    {
        double col_sum = 0.0;
        for (size_t i = 0; i < n; ++i)
            col_sum += std::fabs(A[i][j]);
        if (col_sum > max_sum)
            max_sum = col_sum;
    }
    return max_sum;
}

// Returns the rank of a matrix (row echelon form, for small matrices)
inline size_t rank(const dMatrix& A)
{
    dMatrix M = A;
    const size_t n = M.size(), m = M[0].size();
    size_t r = 0;
    for (size_t col = 0, row = 0; col < m && row < n; ++col)
    {
        size_t sel = row;
        for (size_t i = row + 1; i < n; ++i)
            if (std::fabs(M[i][col]) > std::fabs(M[sel][col]))
                sel = i;
        if (std::fabs(M[sel][col]) < 1e-12)
            continue;
        std::swap(M[sel], M[row]);
        double inv = 1.0 / M[row][col];
        for (size_t j = col; j < m; ++j)
            M[row][j] *= inv;
        for (size_t i = 0; i < n; ++i)
        {
            if (i != row && std::fabs(M[i][col]) > 1e-12)
            {
                double factor = M[i][col];
                for (size_t j = col; j < m; ++j)
                    M[i][j] -= factor * M[row][j];
            }
        }
        ++row;
        ++r;
    }
    return r;
} 
#pragma once
#include "../../typedefs/header.hpp"
// #include "eigen.hpp"

// ------- REQUIRED FUNCTIONS FROM `matrix-properties.hpp` -------
// Returns the transpose_local of a matrix
inline dMatrix transpose_local(const dMatrix& A)
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

// ------- REQUIRED FUNCTIONS FROM `matrix-operations.hpp` -------
// Returns the product of a matrix and a vector
inline dVec matrix_multiply_local(const dMatrix& A, const dVec& x)
{
    const size_t n = A.size();
    if (n == 0) return dVec();
    const size_t m = A[0].size();
    if (x.size() != m)
        throw std::invalid_argument("dMatrix columns must match vector size for multiplication.");
    dVec y(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

// Returns the product of two matrices
inline dMatrix matrix_multiply_local(const dMatrix& A, const dMatrix& B)
{
    const size_t n = A.size();
    if (n == 0) return dMatrix();
    const size_t m = A[0].size();
    const size_t p = B[0].size();
    if (m != B.size())
        throw std::invalid_argument("A's columns must match B's rows for multiplication.");
    dMatrix C(n, dVec(p, 0.0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < p; ++j)
            for (size_t k = 0; k < m; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// ------- REQUIRED FUNCTIONS FROM `eigen.hpp` -------
// Reduces matrix to upper Hessenberg form using Householder reflections
inline dMatrix to_hessenberg_local(dMatrix A)
{
    const size_t n = A.size();
    for (size_t k = 0; k < n - 2; ++k)
    {
        double norm = 0.0;
        for (size_t i = k + 1; i < n; ++i)
            norm += A[i][k] * A[i][k];
        norm = std::sqrt(norm);

        if (norm > 0.0)
        {
            double s = (A[k+1][k] >= 0) ? norm : -norm;
            double u1 = A[k+1][k] + s;
            dVec v(n - k - 1);
            v[0] = u1;
            for (size_t i = k + 2; i < n; ++i)
                v[i-k-1] = A[i][k];
            
            double beta = u1 * s;
            if (beta != 0.0)
            {
                // Apply reflection to A
                for (size_t j = k; j < n; ++j)
                {
                    double sum = 0.0;
                    for (size_t i = k + 1; i < n; ++i)
                        sum += A[i][j] * v[i-k-1];
                    sum /= beta;
                    for (size_t i = k + 1; i < n; ++i)
                        A[i][j] -= sum * v[i-k-1];
                }
                for (size_t i = 0; i < n; ++i)
                {
                    double sum = 0.0;
                    for (size_t j = k + 1; j < n; ++j)
                        sum += A[i][j] * v[j-k-1];
                    sum /= beta;
                    for (size_t j = k + 1; j < n; ++j)
                        A[i][j] -= sum * v[j-k-1];
                }
                A[k+1][k] = -s;
                for (size_t i = k + 2; i < n; ++i)
                    A[i][k] = 0.0;
            }
        }
    }
    return A;
}

// Returns all eigenvalues and eigenvectors using QR iteration with shifts
// Returns pair of (eigenvalues, eigenvectors matrix)
inline std::pair<dVec, dMatrix> eigen_qr_local(dMatrix A, size_t max_iter = 100, double tol = 1e-10)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("dMatrix must be square for eigenvalue computation.");

    // Convert to Hessenberg form
    A = to_hessenberg_local(A);
    
    // Initialize eigenvector matrix as identity
    dMatrix V(n, dVec(n, 0.0));
    for (size_t i = 0; i < n; ++i)
        V[i][i] = 1.0;

    size_t m = n;
    while (m > 1)
    {
        for (size_t iter = 0; iter < max_iter; ++iter)
        {
            // Find Wilkinson shift
            double d = (A[m-2][m-2] - A[m-1][m-1]) / 2.0;
            double mu = A[m-1][m-1] + d - std::copysign(std::sqrt(d*d + A[m-1][m-2]*A[m-2][m-1]), d);
            
            // QR decomposition with implicit shift
            double x = A[0][0] - mu;
            double z = A[1][0];
            
            for (size_t k = 0; k < m - 1; ++k)
            {
                double c = 1.0, s = 0.0;
                if (z != 0.0)
                {
                    double r = std::sqrt(x*x + z*z);
                    c = x/r;
                    s = -z/r;
                }
                
                // Apply Givens rotation
                for (size_t j = k; j < m; ++j)
                {
                    double temp = c * A[k][j] - s * A[k+1][j];
                    A[k+1][j] = s * A[k][j] + c * A[k+1][j];
                    A[k][j] = temp;
                }
                for (size_t i = 0; i <= k + 2 && i < m; ++i)
                {
                    double temp = c * A[i][k] - s * A[i][k+1];
                    A[i][k+1] = s * A[i][k] + c * A[i][k+1];
                    A[i][k] = temp;
                }
                for (size_t i = 0; i < n; ++i)
                {
                    double temp = c * V[i][k] - s * V[i][k+1];
                    V[i][k+1] = s * V[i][k] + c * V[i][k+1];
                    V[i][k] = temp;
                }
                
                if (k < m - 2)
                {
                    x = A[k+1][k];
                    z = A[k+2][k];
                }
            }
            
            // Check for convergence
            if (std::abs(A[m-1][m-2]) < tol * (std::abs(A[m-2][m-2]) + std::abs(A[m-1][m-1])))
            {
                m--;
                break;
            }
            if (iter == max_iter - 1)
                throw std::runtime_error("QR iteration did not converge.");
        }
    }
    
    // Extract eigenvalues and normalize eigenvectors
    dVec eigenvalues(n);
    for (size_t i = 0; i < n; ++i)
    {
        eigenvalues[i] = A[i][i];
        double norm = 0.0;
        for (size_t j = 0; j < n; ++j)
            norm += V[j][i] * V[j][i];
        norm = std::sqrt(norm);
        for (size_t j = 0; j < n; ++j)
            V[j][i] /= norm;
    }
    
    return {eigenvalues, V};
} 
// ------- END OF REQUIRED FUNCTIONS -------


// Contains functions for performing various matrix decompositions.
// Returns LU decomposition with partial pivoting (PA = LU)
// Returns a tuple of {L, U, P} where P is a permutation vector.
inline std::tuple<dMatrix, dMatrix, wVec> lu_decomposition(dMatrix A)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("dMatrix must be square for LU decomposition.");

    wVec P(n);
    for (size_t i = 0; i < n; ++i) P[i] = i;

    for (size_t k = 0; k < n; ++k)
    {
        size_t max_row = k;
        for (size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[max_row][k]))
                max_row = i;
        if (std::abs(A[max_row][k]) < 1e-12)
            throw std::runtime_error("dMatrix is singular or nearly singular.");

        std::swap(A[k], A[max_row]);
        std::swap(P[k], P[max_row]);

        for (size_t i = k + 1; i < n; ++i)
        {
            A[i][k] /= A[k][k];
            for (size_t j = k + 1; j < n; ++j)
                A[i][j] -= A[i][k] * A[k][j];
        }
    }

    dMatrix L(n, dVec(n, 0.0));
    dMatrix U(n, dVec(n, 0.0));
    for (size_t i = 0; i < n; ++i)
    {
        L[i][i] = 1.0;
        for (size_t j = 0; j < n; ++j)
        {
            if (i > j)
                L[i][j] = A[i][j];
            else
                U[i][j] = A[i][j];
        }
    }
    return {L, U, P};
}

// Returns QR decomposition using modified Gram-Schmidt for an m x n matrix A.
// Returns a pair of {Q, R}.
inline std::pair<dMatrix, dMatrix> qr_decomposition(const dMatrix& A)
{
    const size_t m = A.size();
    if (m == 0) return {{}, {}};
    const size_t n = A[0].size();

    dMatrix Q(m, dVec(n, 0.0));
    dMatrix R(n, dVec(n, 0.0));
    dMatrix V = A; // Make a mutable copy of A

    for (size_t j = 0; j < n; ++j)
    {
        double norm_v = 0.0;
        for (size_t i = 0; i < m; ++i)
            norm_v += V[i][j] * V[i][j];
        R[j][j] = std::sqrt(norm_v);

        if (R[j][j] < 1e-12)
            throw std::runtime_error("dMatrix is rank deficient during QR decomposition.");

        for (size_t i = 0; i < m; ++i)
            Q[i][j] = V[i][j] / R[j][j];

        for (size_t k = j + 1; k < n; ++k)
        {
            R[j][k] = 0.0;
            for (size_t i = 0; i < m; ++i)
                R[j][k] += Q[i][j] * V[i][k];

            for (size_t i = 0; i < m; ++i)
                V[i][k] -= R[j][k] * Q[i][j];
        }
    }
    return {Q, R};
}

// Returns Cholesky decomposition L for a symmetric positive definite matrix A (A = LL^T).
inline dMatrix cholesky_decomposition(const dMatrix& A)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("dMatrix must be square for Cholesky decomposition.");

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
                if (std::abs(L[j][j]) < 1e-12)
                    throw std::runtime_error("Zero pivot in Cholesky decomposition.");
                L[i][j] = sum / L[j][j];
            }
        }
    }
    return L;
}

// Returns LDLT decomposition for a symmetric matrix A (A = LDL^T).
// Returns a pair of {L, D} where D is a diagonal vector.
inline std::pair<dMatrix, dVec> ldlt_decomposition(const dMatrix& A)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("dMatrix must be square for LDLT decomposition.");

    dMatrix L(n, dVec(n, 0.0));
    dVec D(n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            double sum = A[i][j];
            for (size_t k = 0; k < j; ++k)
                sum -= L[i][k] * D[k] * L[j][k];
            if (std::abs(D[j]) < 1e-12)
                throw std::runtime_error("Zero pivot in LDLT decomposition.");
            L[i][j] = sum / D[j];
        }
        double sum = A[i][i];
        for (size_t k = 0; k < i; ++k)
            sum -= L[i][k] * L[i][k] * D[k];
        D[i] = sum;
        L[i][i] = 1.0;
    }
    return {L, D};
}

// Returns Singular Value Decomposition (U, S, V) for a square matrix A (A = USV^T).
// Returns a tuple of {U, S, V}, where S is a vector of singular values.
inline std::tuple<dMatrix, dVec, dMatrix> svd(const dMatrix& A)
{
    if (A.empty() || A.size() != A[0].size())
        throw std::invalid_argument("SVD currently supports square matrices only.");

    dMatrix AtA = matrix_multiply_local(transpose_local(A), A);
    auto eigen_pair = eigen_qr_local(AtA);
    dVec s_values = eigen_pair.first;
    dMatrix V = eigen_pair.second;

    for(double& val : s_values)
        val = std::sqrt(std::max(0.0, val));
    
    dMatrix S_inv(s_values.size(), dVec(s_values.size(), 0.0));
    for(size_t i = 0; i < s_values.size(); ++i)
        if (s_values[i] > 1e-12)
            S_inv[i][i] = 1.0 / s_values[i];
    
    dMatrix U = matrix_multiply_local(A, matrix_multiply_local(V, S_inv));
    
    return {U, s_values, V};
}

// Returns the Eigen decomposition (eigenvalues, eigenvectors) for a square matrix A.
inline std::pair<dVec, dMatrix> eigen_decomposition(const dMatrix& A)
{
    if (A.empty() || A.size() != A[0].size())
        throw std::invalid_argument("dMatrix must be square for Eigen decomposition.");
    
    return eigen_qr_local(A);
} 
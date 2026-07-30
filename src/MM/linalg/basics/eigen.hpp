#pragma once
#include "../../typedefs/header.hpp"

// Returns a pair of (eigenvalue, eigenvector) using power method
// For finding the dominant eigenvalue and its corresponding eigenvector
inline std::pair<double, dVec> power_method(const dMatrix& A, size_t max_iter = 100, double tol = 1e-10)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("dMatrix must be square for eigenvalue computation.");

    // Start with random vector
    dVec v(n, 1.0);
    double lambda_old = 0.0;
    
    for (size_t iter = 0; iter < max_iter; ++iter)
    {
        // Normalize vector
        double norm = 0.0;
        for (size_t i = 0; i < n; ++i)
            norm += v[i] * v[i];
        norm = std::sqrt(norm);
        for (size_t i = 0; i < n; ++i)
            v[i] /= norm;

        // Power iteration
        dVec Av(n, 0.0);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                Av[i] += A[i][j] * v[j];

        // Rayleigh quotient for eigenvalue
        double lambda = 0.0;
        for (size_t i = 0; i < n; ++i)
            lambda += v[i] * Av[i];

        // Check convergence
        if (std::abs(lambda - lambda_old) < tol)
            return {lambda, v};

        lambda_old = lambda;
        v = Av;
    }
    throw std::runtime_error("Power method did not converge.");
}

// Reduces matrix to upper Hessenberg form using Householder reflections
inline dMatrix to_hessenberg(dMatrix A)
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
inline std::pair<dVec, dMatrix> eigen_qr(dMatrix A, size_t max_iter = 100, double tol = 1e-10)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("dMatrix must be square for eigenvalue computation.");

    // Convert to Hessenberg form
    A = to_hessenberg(A);
    
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
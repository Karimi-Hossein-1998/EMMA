#pragma once
#include "../../typedefs/header.hpp"
#include "matrix-properties.hpp"

// Returns the dot product of two vectors
inline double dot_product(const dVec& a, const dVec& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Vectors must have the same size for dot product.");
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
    return result;
}

// Returns the cross product of two 3D vectors
inline triple cross_product(const triple& a, const triple& b)
{
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}

// Returns the product of a matrix and a vector
inline dVec matrix_multiply(const dMatrix& A, const dVec& x)
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
inline dMatrix matrix_multiply(const dMatrix& A, const dMatrix& B)
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

// Returns the element-wise product of two matrices
inline dMatrix elementwise_multiply(const dMatrix& A, const dMatrix& B)
{
    const size_t n = A.size();
    if (n == 0) return dMatrix();
    const size_t m = A[0].size();
    if (n != B.size() || m != B[0].size())
        throw std::invalid_argument("Matrices must have the same dimensions for element-wise multiplication.");
    dMatrix C(n, dVec(m));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            C[i][j] = A[i][j] * B[i][j];
    return C;
}

// Returns the element-wise product of two vectors
inline dVec elementwise_multiply(const dVec& a, const dVec& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Vectors must have the same size for element-wise multiplication.");
    dVec c(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        c[i] = a[i] * b[i];
    return c;
}

// Returns the sum of two matrices
inline dMatrix add(const dMatrix& A, const dMatrix& B)
{
    const size_t n = A.size();
    if (n == 0) return dMatrix();
    const size_t m = A[0].size();
    if (n != B.size() || m != B[0].size())
        throw std::invalid_argument("Matrices must have the same dimensions for addition.");
    dMatrix C(n, dVec(m));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

// Returns the sum of two vectors
inline dVec add(const dVec& a, const dVec& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Vectors must have the same size for addition.");
    dVec c(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        c[i] = a[i] + b[i];
    return c;
}

// Returns the difference of two matrices
inline dMatrix subtract(const dMatrix& A, const dMatrix& B)
{
    const size_t n = A.size();
    if (n == 0) return dMatrix();
    const size_t m = A[0].size();
    if (n != B.size() || m != B[0].size())
        throw std::invalid_argument("Matrices must have the same dimensions for subtraction.");
    dMatrix C(n, dVec(m));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

// Returns the difference of two vectors
inline dVec subtract(const dVec& a, const dVec& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("Vectors must have the same size for subtraction.");
    dVec c(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        c[i] = a[i] - b[i];
    return c;
}

// Returns the scalar product of a matrix
inline dMatrix scalar_multiply(double s, const dMatrix& A)
{
    const size_t n = A.size();
    if (n == 0) return dMatrix();
    const size_t m = A[0].size();
    dMatrix C(n, dVec(m));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            C[i][j] = s * A[i][j];
    return C;
}

// Returns the scalar product of a vector
inline dVec scalar_multiply(double s, const dVec& v)
{
    dVec result(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        result[i] = s * v[i];
    return result;
}

// Returns the product of a transposed vector and a matrix (x^T A)
inline dVec vector_transpose_multiply(const dVec& x, const dMatrix& A)
{
    const size_t n = A.size();
    if (n == 0) return dVec();
    const size_t m = A[0].size();
    if (x.size() != n)
        throw std::invalid_argument("Vector size must match matrix rows for x^T A multiplication.");
    dVec result(m, 0.0);
    for (size_t j = 0; j < m; ++j)
        for (size_t i = 0; i < n; ++i)
            result[j] += x[i] * A[i][j];
    return result;
}

// Returns the outer product of two vectors (x ⊗ y or |x⟩⟨y|)
inline dMatrix outer_product(const dVec& x, const dVec& y)
{
    const size_t n = x.size();
    const size_t m = y.size();
    dMatrix result(n, dVec(m));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            result[i][j] = x[i] * y[j];
    return result;
}

// Returns the transpose of a vector as a 1xn matrix
inline dMatrix vector_transpose(const dVec& x)
{
    return dMatrix(1, x);
}

// Returns the outer product (Kronecker product) of two matrices
inline dMatrix matrix_outer_product(const dMatrix& A, const dMatrix& B)
{
    const size_t nA = A.size(), mA = A[0].size();
    const size_t nB = B.size(), mB = B[0].size();
    dMatrix result(nA * nB, dVec(mA * mB));
    for (size_t i = 0; i < nA; ++i)
        for (size_t j = 0; j < mA; ++j)
            for (size_t k = 0; k < nB; ++k)
                for (size_t l = 0; l < mB; ++l)
                    result[i * nB + k][j * mB + l] = A[i][j] * B[k][l];
    return result;
}

// Returns the 2D cross product (scalar result)
inline double cross_product_2d(const dVec& a, const dVec& b)
{
    if (a.size() != 2 || b.size() != 2)
        throw std::invalid_argument("2D cross product requires 2D vectors.");
    return a[0] * b[1] - a[1] * b[0];
}

// Returns a vector perpendicular to two input vectors in n-dimensional space
// Note: For n>3, this is not unique and returns one possible perpendicular vector
inline dVec perpendicular_vector(const dVec& a, const dVec& b)
{
    if (a.size() != b.size() || a.size() < 2)
        throw std::invalid_argument("Vectors must have same size >= 2 for perpendicular vector.");
    
    const size_t n = a.size();
    if (n == 2) return dVec{0.0}; // 2D cross product is scalar
    if (n == 3)
    {
        auto results = cross_product(triple{a[0], a[1], a[2]}, triple{b[0], b[1], b[2]});
        return dVec{results[0], results[1],results[2]};
    } 
    
    // For n>3, find a perpendicular vector using wedge product approach
    // This returns one possible perpendicular vector (not unique)
    dVec result(n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            if (i != j)
            {
                double cofactor = 0.0;
                for (size_t k = 0; k < n; ++k)
                {
                    if (k != i && k != j)
                    {
                        // Simple alternating sign pattern
                        int sign = ((i + j + k) % 2 == 0) ? 1 : -1;
                        cofactor += sign * a[k] * b[k];
                    }
                }
                result[i] += cofactor;
            }
        }
    }
    return result;
}

// Returns the generalized cross product for n-dimensional vectors (n-1 vectors needed)
// This is the proper generalization: for n-dimensional space, need n-1 vectors
inline dVec generalized_cross_product(const dMatrix& vectors)
{
    const size_t n = vectors.size() + 1;
    if (vectors.empty() || vectors[0].size() != n)
        throw std::invalid_argument("For n-dimensional cross product, need n-1 vectors of size n.");
    
    // Create matrix with unit vectors as first row and input vectors as remaining rows
    dMatrix M(n, dVec(n));
    for (size_t j = 0; j < n; ++j)
        M[0][j] = 1.0;  // First row: unit vector
    for (size_t i = 0; i < n - 1; ++i)
            M[i + 1] = vectors[i];
    
    // Compute cross product using cofactors of first row
    dVec result(n);
    for (size_t j = 0; j < n; ++j)
    {
        dMatrix subM;
        for (size_t i = 1; i < n; ++i)
        {
            dVec row;
            for (size_t k = 0; k < n; ++k)
                if (k != j) row.push_back(M[i][k]);
            subM.push_back(row);
        }
        result[j] = ((j % 2 == 0) ? 1 : -1) * determinant(subM);
    }
    return result;
}
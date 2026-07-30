#pragma once
#include "../typedefs/header.hpp"

// Lagrange interpolator: returns a function that interpolates at any x
inline func lagrange_interpolator(const dVec& x_points, const dVec& y_points)
{
    return [x_points, y_points](double x) -> double {
        double result = 0.0;
        size_t n = x_points.size();
        for (size_t i = 0; i < n; ++i)
        {
            double term = y_points[i];
            for (size_t j = 0; j < n; ++j)
            {
                if (i != j)
                    term *= (x - x_points[j]) / (x_points[i] - x_points[j]);
            }
            result += term;
        }
        return result;
    };
}

// Vector-valued Lagrange interpolator: returns a function that interpolates at any x
inline Func lagrange_interpolator(const dVec& x_points, const dMatrix& y_points)
{
    return [x_points, y_points](double x) -> dVec {
        size_t n = x_points.size();
        size_t N = y_points.size();
        dVec result(N, 0.0);
        for (size_t k = 0; k < N; ++k)
        {
            // Interpolate the k-th component
            double value = 0.0;
            for (size_t i = 0; i < n; ++i)
            {
                double term = y_points[k][i];
                for (size_t j = 0; j < n; ++j)
                {
                    if (i != j)
                        term *= (x - x_points[j]) / (x_points[i] - x_points[j]);
                }
                value += term;
            }
            result[k] = value;
        }
        return result;
    };
}

// Barycentric Lagrange interpolator: returns a function that interpolates at any x
inline func barycentric_lagrange_interpolator(const dVec& x_points, const dVec& y_points)
{
    size_t n = x_points.size();
    dVec w(n, 1.0);
    // Compute barycentric weights
    for (size_t j = 0; j < n; ++j)
    {
        for (size_t k = 0; k < n; ++k)
        {
            if (j != k)
                w[j] /= (x_points[j] - x_points[k]);
        }
    }

    return [x_points, y_points, w, n](double x) -> double
    {
        double num = 0.0, denom = 0.0;
        for (size_t j = 0; j < n; ++j)
        {
            if (x == x_points[j])
                return y_points[j];
            double temp = w[j] / (x - x_points[j]);
            num += temp * y_points[j];
            denom += temp;
        }
        return num / denom;
    };
}

// Vector-valued Barycentric Lagrange interpolator: returns a function that interpolates at any x
inline Func barycentric_lagrange_interpolator(const dVec& x_points, const dMatrix& y_points)
{
    size_t n = x_points.size();
    size_t N = y_points.size();
    dVec w(n, 1.0);
    // Compute barycentric weights
    for (size_t j = 0; j < n; ++j)
    {
        for (size_t k = 0; k < n; ++k)
        {
            if (j != k)
                w[j] /= (x_points[j] - x_points[k]);
        }
    }
    return [x_points, y_points, w, n, N](double x) -> dVec {
        dVec result(N, 0.0);
        for (size_t k = 0; k < N; ++k)
        {
            double num = 0.0, denom = 0.0;
            bool found = false;
            for (size_t j = 0; j < n; ++j)
            {
                if (x == x_points[j])
                {
                    result[k] = y_points[k][j];
                    found = true;
                    break;
                }
                double temp = w[j] / (x - x_points[j]);
                num += temp * y_points[k][j];
                denom += temp;
            }
            if (!found)
                result[k] = num / denom;
        }
        return result;
    };
}


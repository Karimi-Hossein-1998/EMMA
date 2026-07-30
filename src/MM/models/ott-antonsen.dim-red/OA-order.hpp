#pragma once
#include "../../typedefs/header.hpp"

inline std::tuple<dMatrix, dMatrix, dMatrix> OA_order(
    const SolverResults& results,
    size_t               num_modules
)
{
    dMatrix sol        = results.solution;
    size_t  length     = sol.size();
    size_t  width      = sol[0].size();
    dMatrix order      = dMatrix(length, dVec(3*num_modules+4, 0.0));
    dMatrix order_sqr  = dMatrix(length, dVec(num_modules+2, 0.0));
    dMatrix magnitudes = dMatrix(length,dVec(num_modules+2,0.0));

    double dummy_order_sqr = 0.0, sum_im = 0.0, sum_re = 0.0;
    for ( size_t i = 0; i < length; ++i )
    {
        sum_re = 0.0;
        sum_im = 0.0;
        for ( size_t j = 0; j < num_modules; ++j)
        {
            dummy_order_sqr = sol[i][2*j]*sol[i][2*j] + sol[i][2*j+1]*sol[i][2*j+1];

            order[i][3*j]   = sol[i][2*j];
            order[i][3*j+1] = sol[i][2*j+1];
            order[i][3*j+2] = std::sqrt(dummy_order_sqr);

            magnitudes[i][j] = order[i][3*j+2];

            order_sqr[i][j] = dummy_order_sqr;

            sum_re         += sol[i][2*j];
            sum_im         += sol[i][2*j+1];
        }
        sum_re                     /= num_modules;
        sum_im                     /= num_modules;
        dummy_order_sqr             = sum_re*sum_re + sum_im*sum_im;

        order[i][3*num_modules]     = sum_re;
        order[i][3*num_modules+1]   = sum_im;
        order[i][3*num_modules+2]   = std::sqrt(dummy_order_sqr);
        order[i][3*num_modules+3]   = results.time_points[i];

        magnitudes[i][num_modules]   = order[i][3*num_modules+2];
        magnitudes[i][num_modules+1] = order[i][3*num_modules+3];
 
        order_sqr[i][num_modules]   = dummy_order_sqr;
        order_sqr[i][num_modules+1] = results.time_points[i];
    }
    return {magnitudes, order_sqr, order};
}

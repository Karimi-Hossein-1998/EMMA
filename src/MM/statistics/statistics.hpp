#pragma once
#include "../typedefs/header.hpp"
#include "../linalg/basics/matrix-properties.hpp"

inline double mean_value(const dVec& data)
{
    size_t N = data.size();
    double mean_val = 0.0;
    for ( size_t i = 0; i < N; ++i )
    {
        mean_val += data[i];
    }
    mean_val /= N;
    return mean_val;
}

inline dVec mean_value0(const dMatrix& data)
{
    size_t length = data.size();
    dVec mean_vals(length, 0.0);
    for ( size_t i = 0; i < length; ++i)
    {
        mean_vals[i] = mean_value(data[i]);
    }
    return mean_vals;
}

inline dVec mean_value1(const dMatrix& data)
{
    dMatrix data_transpose = transpose(data);
    size_t length = data_transpose.size();
    dVec mean_vals(length, 0.0);
    for ( size_t i = 0; i < length; ++i)
    {
        mean_vals[i] = mean_value(data_transpose[i]);
    }
    return mean_vals;
}

// mode = 0 a vector with mean values calculated over rows
// mode = 1 a vector with mean values calculated over columnss
inline dVec mean_value(const dMatrix& data, size_t mode)
{
    if ( mode == 0 )
    {
        return mean_value0(data);
    }
    else if ( mode == 1 )
    {
        return mean_value1(data);
    }
    else
    {
        throw std::invalid_argument("The Value for mode ("+std::to_string(mode)+") is invalid!\n");
    }
}

inline double variance(const dVec& data)
{
    double mean_local = 0.0;
    mean_local = mean_value(data);

    size_t N = data.size();
    dVec sqr_data = dVec(N, 0.0);
    for ( size_t i = 0; i < N; ++i )
    {
        sqr_data[i] = data[i]*data[i];
    }
    double sqr_mean = mean_value(sqr_data);

    double sigma2 = sqr_mean - mean_local*mean_local;
    return sigma2;
}

inline dVec variance0(const dMatrix& data)
{
    size_t length   = data.size();
    dVec   var_vals = dVec(length, 0.0);
    for ( size_t i = 0; i < length; ++i )
    {
        var_vals[i]   = variance(data[i]);
    }
    return var_vals;
}

inline dVec variance1(const dMatrix& data)
{
    dMatrix data_transpose = transpose(data);
    size_t  length         = data.size();
    dVec    var_vals       = dVec(length, 0.0);
    for ( size_t i = 0; i < length; ++i )
    {
        var_vals[i]   = variance(data_transpose[i]);
    }
    return var_vals;
}

inline dVec variance(const dMatrix& data, size_t mode)
{
    if ( mode == 0 )
    {
        return variance0(data);
    }
    else if ( mode == 1 )
    {
        return variance1(data);
    }
    else
    {
        throw std::invalid_argument("The Value for mode ("+std::to_string(mode)+") is invalid!\n");
    }
}

inline double stddev(const dVec& data)
{
    double sigma2_local = 0.0;
    sigma2_local = variance(data);

    double sigma = std::sqrt(sigma2_local);
    return sigma;
}

inline dVec stddev0(const dMatrix& data)
{
    size_t length      = data.size();
    dVec   stddev_vals = dVec(length, 0.0);
    for ( size_t i = 0; i < length; ++i )
    {
        stddev_vals[i] = stddev(data[i]);
    }
    return stddev_vals;
}

inline dVec stddev1(const dMatrix& data)
{
    dMatrix data_transpose = transpose(data);
    size_t  length         = data.size();
    dVec    stddev_vals    = dVec(length, 0.0);
    for ( size_t i = 0; i < length; ++i )
    {
        stddev_vals[i] = stddev(data_transpose[i]);
    }
    return stddev_vals;
}

inline dVec stddev(const dMatrix& data, size_t mode)
{
    if ( mode == 0 )
    {
        return stddev0(data);
    }
    else if ( mode == 1 )
    {
        return stddev1(data);
    }
    else
    {
        throw std::invalid_argument("The Value for mode ("+std::to_string(mode)+") is invalid!\n");
    }
}

inline double coefficient_of_variability(const dVec& data)
{
    double local_mean = 0.0;
    local_mean = mean_value(data);
    double local_sigma = 0.0;
    local_sigma = stddev(data);
    double cv = local_sigma/local_mean;
    return cv;
}

inline dVec coefficient_of_variability0(const dMatrix& data)
{
    size_t length = data.size();
    dVec   cv     = dVec(length, 0.0);
    for ( size_t i = 0; i < length; ++i )
    {
        cv[i] = coefficient_of_variability(data[i]);
    }
    return cv;
}

inline dVec coefficient_of_variability1(const dMatrix& data)
{
    dMatrix data_transpose = transpose(data);
    size_t  length         = data.size();
    dVec    cv             = dVec(length, 0.0);
    for ( size_t i = 0; i < length; ++i )
    {
        cv[i] = coefficient_of_variability(data_transpose[i]);
    }
    return cv;
}

inline dVec coefficient_of_variability(const dMatrix& data, size_t mode)
{
    if ( mode == 0 )
    {
        return coefficient_of_variability0(data);
    }
    else if ( mode == 1 )
    {
        return coefficient_of_variability1(data);
    }
    else
    {
        throw std::invalid_argument("The Value for mode ("+std::to_string(mode)+") is invalid!\n");
    }
}

inline dVec statistics(const dVec& data)
{
    size_t length   = data.size();
    double sum      = 0.0;
    double sum2     = 0.0;
    double stddev2  = 0.0;
    double stddev   = 0.0;
    double cv       = 0.0;
    for ( size_t i = 0; i < length; ++i )
    {
        sum  += data[i];
        sum2 += data[i]*data[i];
    }
    sum  /= length;
    sum2 /= length;

    stddev2    = sum2 - sum*sum;
    stddev     = std::sqrt(stddev2);
    cv         = stddev / sum;

    return {sum, stddev, cv};
}

inline dMatrix statistics0(const dMatrix& data)
{
    size_t  length = data.size();
    dMatrix stats  = dMatrix(length, {0.0,0.0,0.0});
    for ( size_t i = 0; i < length; ++i )
    {
        stats[i]      = statistics(data[i]);
    }

    return stats;
}

inline dMatrix statistics1(const dMatrix& data)
{
    dMatrix data_transpose = transpose(data);
    size_t  length         = data_transpose.size();
    dMatrix stats          = dMatrix(length, {0.0,0.0,0.0});
    for ( size_t i = 0; i < length; ++i )
    {
        stats[i]      = statistics(data_transpose[i]);
    }

    return stats;
}

inline dMatrix statistics(const dMatrix& data, size_t mode)
{
    if ( mode == 0 )
    {
        return statistics0(data);
    }
    else if ( mode == 1 )
    {
        return statistics1(data);
    }
    else
    {
        throw std::invalid_argument("The Value for mode ("+std::to_string(mode)+") is invalid!\n");
    }
}

#pragma once
#include "../typedefs/header.hpp"
#include "../linalg/basics/matrix-properties.hpp"

inline dMatrix statistics_last_n_percent(const dVec& data, const dVec& percentages)
{
    size_t  N_p     = percentages.size();
    size_t  length  = data.size();
    dMatrix stats   = dMatrix(N_p, dVec(3, 0.0));
    wVec    entries = wVec(N_p, 0);

    double sum     = 0.0;
    double sum2    = 0.0;
    double stddev2 = 0.0;
    double stddev  = 0.0;
    double cv      = 0.0;
    for ( size_t i = 0; i < N_p; ++i)
    {
        entries[i] = static_cast<size_t>(percentages[i]*length);
        sum        = 0.0;
        sum2       = 0.0;
        stddev2    = 0.0;
        stddev     = 0.0;
        cv         = 0.0;
        size_t N_i = length - entries[i];
        for ( size_t k = N_i; k < length; ++k )
        {
            sum  += data[k];
            sum2 += data[k]*data[k];
        }
        sum  /= entries[i];
        sum2 /= entries[i];

        stddev2  = sum2 - sum*sum;
        stddev   = std::sqrt(stddev2);
        cv       = stddev / sum;
        stats[i] = {sum, stddev, cv};
    }
    return stats;
}

inline Vec<dMatrix> statistics_last_n_percent0(const dMatrix& data, const dVec& percentages)
{
    size_t length = data.size();
    Vec<dMatrix> stats = Vec(length, dMatrix());

    for ( size_t i = 0; i < length; ++i )
    {
        stats[i] = statistics_last_n_percent(data[i],percentages);
    }
    return stats;
}

inline Vec<dMatrix> statistics_last_n_percent1(const dMatrix& data, const dVec& percentages)
{
    dMatrix data_transpose = transpose(data);
    size_t  length         = data_transpose.size();
    Vec<dMatrix> stats = Vec(length, dMatrix());

    for ( size_t i = 0; i < length; ++i )
    {
        stats[i] = statistics_last_n_percent(data_transpose[i],percentages);
    }
    return stats;
}

inline Vec<dMatrix> statistics_last_n_percent(const dMatrix& data, const dVec& percentages, size_t mode)
{
    if ( mode == 0 )
    {
        return statistics_last_n_percent0(data,percentages);
    }
    else if ( mode == 1 )
    {
        return statistics_last_n_percent1(data,percentages);
    }
    else
    {
        throw std::invalid_argument("The Value for mode ("+std::to_string(mode)+") is invalid!\n");
    }
}

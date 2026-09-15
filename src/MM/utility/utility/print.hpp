#pragma once
#include "../typedefs/header.hpp"

// -----------------------------------------------------------------------------
// Utility functions for printing matrices/vectors
//
// Note: These are useful for debugging and data export. For large matrices/vectors,
// consider writing to file instead of printing to console to avoid overwhelming output.

// Print a vector to std::cout (space-separated)
template <typename t>
void print_vector(const Vec<t>& vector)
{
    size_t size = vector.empty() ? 0 : vector.size();
    std::cout << "Vector (" << size << "):" << std::endl;
    std::cout << "[ ";
    for (size_t v = 0; v < vector.size(); ++v)
    {
        std::cout << vector[v];
        if ( v + 1 < vector.size()) std::cout << ", ";
    }
    std::cout << " ]";
    std::cout << std::endl;
}

template <typename t>
void print_matrix(const Vec<Vec<t>>& matrix)
{
    std::cout << "dMatrix (" << matrix.size() << "x" << (matrix.empty() ? 0 : matrix[0].size()) << "): " << std::endl;
    std::cout << "[ ";
    for (size_t r = 0; r < matrix.size(); ++r)
    {
        const auto& row = matrix[r];
        std::cout << "[ ";
        for (size_t i = 0; i < row.size(); ++i)
        {
            std::cout << row[i];
            if (i + 1 < row.size()) std::cout << ", ";
        }
        std::cout << " ]";
        if (r + 1 < matrix.size()) std::cout << ",";
        std::cout << std::endl;
    }
    std::cout << " ]" << std::endl;
}


// Call this function to print the progress bar
template <typename T>
void print_progress_bar_eta(size_t current, size_t total, T eta, size_t bar_width = 50)
{
    double progress = static_cast<double>(current) / total;
    size_t pos = static_cast<size_t>(bar_width * progress);
    std::cout << "\rProgress: [";
    for (size_t i = 0; i < bar_width; ++i)
    {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << ".";
    }
    std::cout << "] ";
    auto percentage = static_cast<size_t>(progress * 100.0);
    if ( percentage != 100)
        std::cout << std::setw(3) << percentage << "%";
    else
        std::cout << "Done!";
    std::cout << " : eta = " << eta;
    std::cout.flush();
}

void print_progress_bar(size_t current, size_t total, size_t bar_width = 50)
{
    double progress = static_cast<double>(current) / total;
    size_t pos = static_cast<size_t>(bar_width * progress);
    std::cout << "\rProgress: [";
    for (size_t i = 0; i < bar_width; ++i)
    {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << ".";
    }
    std::cout << "] ";
    auto percentage = static_cast<size_t>(progress * 100.0);
    if ( percentage != 100)
        std::cout << std::setw(3) << percentage << "%";
    else
        std::cout << "Done!";
    std::cout.flush();
}
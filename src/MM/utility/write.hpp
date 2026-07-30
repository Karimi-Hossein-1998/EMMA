#pragma once
#include "../typedefs/header.hpp"

// -----------------------------------------------------------------------------
// Utility functions for writing matrices/vectors
//
// Note: These are useful for debugging and data export. For large matrices/vectors,
// consider writing to file instead of printing to console to avoid overwhelming output.


// Write a matrix to a file with flexible formatting options
// Parameters:
//   matrix:      dMatrix to write
//   directory:   Output directory
//   filename:    Output filename (without extension)
//   filetype:    File extension (default: "txt")
//   separator:   Separator between values (default: " ")
//   comment:     Optional comment string (added after header)
//   precision:   String for precision (number of digits after decimal)
//   format:      String for format ("fixed" or "scientific")
//   alignment:   String for alignment ("left", "right", or "")
//   header:      Optional header line
//   footer:      Optional footer line
//   append:      If true, append to file instead of overwriting (default: false)
//   col_width:   Column width for alignment (default: 0 = no width)
//   skip_trailing_separator: If true, do not write separator after last value in row (default: false)
//   binary:      If true, write as raw binary (default: false)
//
// Notes:
//   - If separator is ',' or '\t', values containing separator or special chars are quoted (CSV/TSV escaping)
template <typename t>
void write_matrix_to_file(
    const Vec<Vec<t>>& matrix, 
    const std::string& directory, 
    const std::string& filename, 
    const std::string& filetype                = "txt", 
    const std::string& separator               = " ", 
    const std::string& comment                 = "", 
    const std::string& precision               = "16", 
    const std::string& format                  = "fixed", 
    const std::string& alignment               = "", 
    const std::string& header                  = "", 
    const std::string& footer                  = "",
    bool               append                  = false,
    size_t             col_width               = 0,
    bool               skip_trailing_separator = false,
    bool               binary                  = false
)
{
    std::ios_base::openmode mode = std::ios::out | (append ? std::ios::app : std::ios::trunc);
    if (binary) mode |= std::ios::binary;
    std::ofstream file(directory+"/"+filename+"."+filetype, mode);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + directory + "/" + filename + "." + filetype);
    }
    if (binary) 
    {
        // Write matrix dimensions first
        size_t rows = matrix.size();
        size_t cols = rows ? matrix[0].size() : 0;
        file.write(reinterpret_cast<const char*>(&rows), sizeof(size_t));
        file.write(reinterpret_cast<const char*>(&cols), sizeof(size_t));
        for (const auto& row : matrix) 
        {
            file.write(reinterpret_cast<const char*>(row.data()), sizeof(double) * row.size());
        }
        file.close();
        return;
    }
    // Set floating-point format
    if (format == "fixed")
        file << std::fixed;
    else if (format == "scientific")
        file << std::scientific;
    // Set precision
    file << std::setprecision(std::stoi(precision));
    // Alignment (left/right/none)
    if (alignment == "left")
        file << std::left;
    else if (alignment == "right")
        file << std::right;
    // Header and comment
    if (!header.empty())
    {
        if (comment.empty())
        {
            file << header << std::endl;
        }
        else
        {
            file << "# " << header << std::endl;
            file << "## " << comment << std::endl;
        }
    }
    auto needs_quoting = [&](const std::string& s) {
        return (separator == "," || separator == "\t") && (s.find(separator) != std::string::npos || s.find('"') != std::string::npos || s.find('\n') != std::string::npos);
    };
    for (const auto& row : matrix)
    {
        for (size_t i = 0; i < row.size(); ++i)
        {
            std::ostringstream oss;
            if (col_width > 0) oss << std::setw(col_width);
            oss << row[i];
            std::string sval = oss.str();
            if (needs_quoting(sval)) 
            {
                std::string quoted = '"' + sval + '"';
                sval = quoted;
            }
            file << sval;
            if (!(skip_trailing_separator && i == row.size() - 1))
            {
                file << separator;
            }
        }
        file << std::endl;
    }
    if (!footer.empty())
    {
        file << footer << std::endl;
    }
    file.close();
}

// Write a vector to a file with flexible formatting options
// Parameters:
// #  vector:      Vector to write
// #  directory:   Output directory
// #  filename:    Output filename (without extension)
// #  filetype:    File extension (default: "txt")
// #  separator:   Separator between values (default: " ")
// #  comment:     Optional comment string (added after header)
// #  precision:   String for precision (number of digits after decimal)
// #  format:      String for format ("fixed" or "scientific")
// #  alignment:   String for alignment ("left", "right", or "")
// #  header:      Optional header line
// #  footer:      Optional footer line
// #  append:      If true, append to file instead of overwriting (default: false)
// #  col_width:   Column width for alignment (default: 0 = no width)
// #  skip_trailing_separator: If true, do not write separator after last value (default: false)
// #  binary:      If true, write as raw binary (default: false)
//
// Notes:
//   - If separator is ',' or '\t', values containing separator or special chars are quoted (CSV/TSV escaping)
template <typename t>
void write_vector_to_file(
    const Vec<t>&      vector,
    const std::string& directory,
    const std::string& filename,
    const std::string& filetype                = "txt",
    const std::string& separator               = " ",
    const std::string& comment                 = "",
    const std::string& precision               = "16",
    const std::string& format                  = "fixed",
    const std::string& alignment               = "",
    const std::string& header                  = "",
    const std::string& footer                  = "",
    bool               append                  = false,
    size_t             col_width               = 0,
    bool               skip_trailing_separator = false,
    bool               binary                  = false
)
{
    std::ios_base::openmode mode = std::ios::out | (append ? std::ios::app : std::ios::trunc);
    if (binary) mode |= std::ios::binary;
    std::ofstream file(directory+"/"+filename+"."+filetype, mode);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + directory + "/" + filename + "." + filetype);
    }
    if (binary) {
        // Write vector size first
        size_t n = vector.size();
        file.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
        file.write(reinterpret_cast<const char*>(vector.data()), sizeof(double) * n);
        file.close();
        return;
    }
    // Set floating-point format
    if (format == "fixed")
        file << std::fixed;
    else if (format == "scientific")
        file << std::scientific;
    // Set precision
    file << std::setprecision(std::stoi(precision));
    // Alignment (left/right/none)
    if (alignment == "left")
        file << std::left;
    else if (alignment == "right")
        file << std::right;
    // Header and comment
    if (!header.empty())
    {
        if (comment.empty())
        {
            file << header << std::endl;
        }
        else
        {
            file << "# " << header << std::endl;
            file << "## " << comment << std::endl;
        }
    }
    auto needs_quoting = [&](const std::string& s) {
        return (separator == "," || separator == "\t") && (s.find(separator) != std::string::npos || s.find('"') != std::string::npos || s.find('\n') != std::string::npos);
    };
    for (size_t i = 0; i < vector.size(); ++i)
    {
        std::ostringstream oss;
        if (col_width > 0) oss << std::setw(col_width);
        oss << vector[i];
        std::string sval = oss.str();
        if (needs_quoting(sval)) {
            std::string quoted = '"' + sval + '"';
            sval = quoted;
        }
        file << sval;
        if (!(skip_trailing_separator && i == vector.size() - 1))
            file << separator;
    }
    if (!footer.empty())
    {
        file << footer << std::endl;
    }
    file.close();
}

struct Write_Read
{
    std::string filename;
    std::string filepath;
    std::string header;
    std::string footer;
    std::string separator;
};

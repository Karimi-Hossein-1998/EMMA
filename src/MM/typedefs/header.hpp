#pragma once
#include <cassert>
#include <vector>
#include <thread>
#include <array>
#include <functional>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <ostream>
#include <complex>

#ifdef PI
#undef PI
#endif
namespace MathEngine
{ // MathEngine namespace
// -----------------------------------------------------------------------------
// Constant expression for pi
constexpr double PI        = 3.14159265358979323846;
constexpr int    max_order = 10;
constexpr double myexp     = 2.7;

// -----------------------------------------------------------------------------
// Type aliases for clarity and convenience
template<typename T>                                                              // A template for std::vector
using Vec            = std::vector<T>;                                            // that is basically a renaming (Vec = std::vector)
template<typename U>
using Matrix         = std::vector<std::vector<U>>;
using dVec           = Vec<double>;                                               // 1D vector of doubles
using wVec           = Vec<size_t>;
using dcVec          = Vec<std::complex<double>>;                                 // 1D vector of complex doubles
using iVec           = Vec<int>;                                                  // 1D vector of integers
using uVec           = Vec<unsigned int>;                                         // 1D vector of unsigned integers
using bVec           = Vec<bool>;                                                 // 1D vector of booleans
using sVec           = Vec<std::string>;                                          // 1D vector of strings
using dpairVec       = Vec<std::pair<double, double>>;                            // 1D vector of (double, double) pairs
using spairVec       = Vec<std::pair<std::string, std::string>>;                  // 1D vector of (string, string) pairs
using triple         = std::array<double, 3>;
using dMatrix        = Matrix<double>;                                            // 2D matrix (vector of vectors) -- double
using dcMatrix       = Matrix<std::complex<double>>;                              // 2D matrix (vector of vectors) -- std::complex
using wMatrix        = Matrix<size_t>;
using iMatrix        = Matrix<int>;
using uMatrix        = Matrix<unsigned int>;
using bMatrix        = Matrix<bool>;
using sMatrix        = Matrix<std::string>;
using func           = std::function<double(double)>;                             // Function type taking double and returning double
using Func           = std::function<dVec(double)>;                               // Function type taking double and returning dVec
using MyFunc         = std::function<dVec(double, const dVec&)>;                  // ODE function type taking (double, dVec) and returning dVec
using DelayFunc      = std::function<dVec(double, const dVec&, const dVec&)>;     // DDE function type (single delay) taking (double, dVec, dVec) and returning dVec
using DelayFuncMulti = std::function<dVec(double, const dVec&, const dMatrix&)>;   // DDE function type (multiple delays) taking (double, dVec, dMatrix) and returning dVec

using interPolater   = std::function<func(const dVec&, const dVec&)>;
using InterPolater   = std::function<Func(const dVec&, const dMatrix&)>;

// Type definition for solvers
// -------------------------------------------------------------------------------------------------------------
struct SolverParameters 
{
    MyFunc derivative;                           // ODE function
    dVec   initial_conditions;
    double t0                    = 0.0;         // Start time
    double t1                    = 10.0;        // End time
    double dt                    = 1.0 / 16.0;  // Initial step size
    int    order                 = 4;           // Order for AB/ABM
    int    iterations            = 2;           // Iterations for ABM
    double min_dt                = 1.0 / 128.0; // Minimum step size (for variable step size control)
    double max_dt                = 0.5;         // Maximum step size (for variable step size control)
    double decrease_factor       = 0.9;
    double increase_factor       = 2.0;
    double local_tol             = 1.0e-7;      // Local error tolerance
    double absolute_tol          = 1.0e-10;     // Takes care of division by zero (NaN and/or inf values)
    double local_tol_error_ratio = 0.1;  
    size_t max_trial             = 25;
    bool   error_estimate        = false;      // Wether to use error estimation (or not)
    bool   variable_steps        = false;      // Wether to use variable step size (or not)
    bool   attempts_history      = false;      // Controls wether the errors and step sizes are stored for all attempts or not.
    bool   weighted_error        = true;       // Calculate error based on weighted error formula
    bool   norm_error            = false;      // Calculate error based on norm error formula (norm or the solution is used instead of maximum of the absolute value of the elements of the solution)
    // DDE-specific parameters
    DelayFunc delay_derivative            = nullptr; // DDE function f(t, y, y_delayed)
    DelayFuncMulti delay_derivative_multi = nullptr; // DDE function f(t, y, y_delayed_states) for multiple delays
    dVec      delay_times;                           // Vector of delay times τ_i
    Func      history_function;                      // Initial history function y(t) for t < t0
    bool      is_dde                      = false;   // Flag to indicate DDE vs ODE
    size_t    max_delay_order             = 10;      // Maximum order for delay terms
    double    interpolation_tol           = 1.0e-6;  // Tolerance for delay interpolation
    double    dt_scale                    = 1.0e-1;  // tau_min is compared to dt_scale*dt to determine a fine time step.
};

struct SolverResults 
{
    dMatrix      solution;       // The solution matrix.
    dVec        time_points;    // Time points where solution is computed.
    dVec        errors;         // Error estimates.
    dVec        step_sizes;     // Actual step sizes used.
    dMatrix      steps_history;  // History of step sizes used (contains failed steps attempted).
    dMatrix      errors_history; // History of errors (contains failed steps attempted).
    wVec failed_trials;  // Number of failed trials in each step.
};

// -----------------------------------------------------------------------------
// Sparse vector and matrix structures for efficient storage of sparse data
struct SparsedVec  
{
    std::vector<std::pair<size_t, double>> values; // (index, value) pairs
    SparsedVec(size_t N) : values(N) {}
    void add(size_t i, double value) { values.emplace_back(i, value); }
};

struct SparsedMatrix {
    std::vector<std::vector<std::pair<size_t, double>>> rows; // Each row: (col, value) pairs
    SparsedMatrix(size_t N) : rows(N) {}
    void add(size_t i, size_t j, double weight) { rows[i].emplace_back(j, weight); }
};

} // End MathEngine namespace

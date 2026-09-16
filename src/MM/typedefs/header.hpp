#pragma once
#include "NumCpp/NumCpp.hpp"
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
static constexpr const double PI        = 3.14159265358979323846;
static constexpr const std::uint8_t   max_order = 10;

// -----------------------------------------------------------------------------
// Type aliases for clarity and convenience
template<typename T>                                                              // A template for std::vector
using Vec       = std::vector<T>;                                            // that is basically a renaming (Vec = std::vector)
using dVec      = Vec<double>;                                               // 1D vector of doubles
using wVec      = Vec<size_t>;
using dcVec     = Vec<std::complex<double>>;                                 // 1D vector of complex doubles
using iVec      = Vec<int>;                                                  // 1D vector of integers
using uVec      = Vec<unsigned int>;                                         // 1D vector of unsigned integers
using bVec      = Vec<bool>;                                                 // 1D vector of booleans
using sVec      = Vec<std::string>;                                          // 1D vector of strings
using wpairVec  = Vec<std::pair<size_t,size_t>>;
using dpairVec  = Vec<std::pair<double, double>>;                            // 1D vector of (double, double) pairs
using spairVec  = Vec<std::pair<std::string, std::string>>;                  // 1D vector of (string, string) pairs

using dMatrix = Matrix<double>;


// Call Back Function to get data on the fly
struct OneStepSolverResult
{
    Vec<double> sol;
    Vec<double> stepHistory;
    Vec<double> errorHistory;
    double      timePoint;
    double      error;
    double      stepSize;
    size_t      failedTrial;
};
using func           = std::function<double(double)>;
using Func           = std::function<dVec(double)>;
using MyFunc         = std::function<void(double, const dVec&, dVec&)>;
using DelayFunc      = std::function<void(double, const dVec&, const dVec&, dVec&)>;
using MultiDelayFunc = std::function<void(double, const dVec&, const dMatrix&, dVec&)>;
using CallBackFunc   = std::function<void(const OneStepSolverResult&)>;
struct ODESolverParameters
{
    MyFunc       derivative;                           // ODE function
    CallBackFunc onStep = nullptr;
    Vec<double>  initialConditions;
    double       t0                 = static_cast<double>(0.0);         // Start time
    double       t1                 = static_cast<double>(10.0);        // End time
    double       dt                 = static_cast<double>(1.0 / 16.0);  // Initial step size
    double       minDt              = static_cast<double>(1.0 / 128.0); // Minimum step size (for variable step size control)
    double       maxDt              = static_cast<double>(0.5);         // Maximum step size (for variable step size control)
    double       decreaseFactor     = static_cast<double>(0.9);
    double       increaseFactor     = static_cast<double>(2.0);
    double       localTol           = static_cast<double>(1.0e-7);      // Local error tolerance
    double       absolute_tol       = static_cast<double>(1.0e-10);     // Takes care of division by zero (NaN and/or inf values)
    double       localTolErrorRatio = static_cast<double>(0.1);
    std::uint8_t order              = 4;           // Order for AB/ABM
    std::uint8_t iterations         = 2;           // Iterations for ABM
    std::uint8_t maxTrial           = 25;
    bool         errorEstimate      = false;      // Wether to use error estimation (or not)
    bool         variableSteps      = false;      // Wether to use variable step size (or not)
    bool         attemptsHistory    = false;      // Controls wether the errors and step sizes are stored for all attempts or not.
    bool         weightedError      = true;       // Calculate error based on weighted error formula
    bool         normError          = false;      // Calculate error based on norm error formula (norm or the solution is used instead of maximum of the absolute value of the elements of the solut[...]
};

struct SolverResults 
{
    Matrix<double> solution;      // The solution matrix.
    Matrix<double> stepsHistory;  // History of step sizes used (contains failed steps attempted).
    Matrix<double> errorsHistory; // History of errors (contains failed steps attempted).
    Vec<double> timePoints;       // Time points where solution is computed.
    Vec<double> errors;           // Error estimates.
    Vec<double> stepSizes;        // Actual step sizes used.
    wVec failedTrials;            // doubleber of failed trials in each step.
};

using SolverFunc = std::function<SolverResults(const ODESolverParameters&)>;

struct SparsedMatrix
{
    Vec<Vec<std::pair<std::uint64_t, double>>> rows; // Each row: (col, value) pairs
    SparsedMatrix(std::uint64_t N=0) : rows(N) {}
    void add(std::uint64_t i, std::uint64_t j, double weight) { rows[i].emplace_back(j, weight); }
};

} // End MathEngine namespace

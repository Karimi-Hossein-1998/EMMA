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
static constexpr const std::float64_t PI        = 3.14159265358979323846;
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

// Vec<T> = std::vector<T>
// Matrix = a container for contiguous memory designed for cache safety.
// FPNumber = concept for Floating Point Number(s) (A concept for arithmetic and complex types);
template <typename F, typename Num>
concept func = FPNumber<Num> && std::invocable<F, Num> && std::same_as<std::invoke_result_t<F, Num>, Num>;

template <typename F, typename Num>
concept Func = FPNumber<Num> && std::invocable<F, Num> && std::same_as<std::invoke_result_t<F, Num>, Vec<Num>>;

template <typename F, typename Num>
concept MyFunc = FPNumber<Num> && std::invocable<F, Num, const Vec<Num>&, Vec<Num>&> && std::same_as<std::invoke_result_t<F, Num, const Vec<Num>&, Vec<Num>&>, void>;

template <typename F, typename Num>
concept DelayFunc = FPNumber<Num> && std::invocable<F, Num, const Vec<Num>&, const Vec<Num>&, Vec<Num>&> &&  std::same_as<std::invoke_result_t<F, Num, const Vec<Num>&, const Vec<Num>&, Vec<Num>&>, void>;

template <typename F, typename Num>
concept MultiDelayFunc = FPNumber<Num> && std::invocable<F, Num, const Vec<Num>&, const Matrix<Num>&, Vec<Num>&> && std::same_as<std::invoke_result_t<F, Num, const Vec<Num>&, const Matrix<Num>&, Vec<Num>&>, void>;

// using interPolater = std::function<func(const dVec&, const dVec&)>;
// using InterPolater = std::function<Func(const dVec&, const dMatrix&)>;

// Call Back Function to get data on the fly
template <FPNumber Num>
struct OneStepSolverResult
{
    Vec<Num> sol;
    Vec<Num> stepHistory;
    Vec<Num> errorHistory;
    Num timePoint;
    Num error;
    Num stepSize;
    size_t failedTrial;
};

template <typename F, typename Num>
concept CallBackFunc = FPNumber<Num> && std::invocable<F, const OneStepSolverResult<Num>&> && std::same_as<std::invoke_result_t<F, const OneStepSolverResult<Num>&>, void>;

struct DefaultCallBackFunc
{
    template <typename ResulT>
    void operator()(const ResulT&) const noexcept {}
};

// Type definition for solvers
// -------------------------------------------------------------------------------------------------------------
template <FPNumber Num, MyFunc<Num> MF, CallBackFunc<Num> CBF=DefaultCallBackFunc>
struct ODESolverParameters
{
    MF  derivative;                           // ODE function
    CBF onStep;
    Vec<Num> initialConditions;
    Num t0                 = static_cast<Num>(0.0);         // Start time
    Num t1                 = static_cast<Num>(10.0);        // End time
    Num dt                 = static_cast<Num>(1.0 / 16.0);  // Initial step size
    Num minDt              = static_cast<Num>(1.0 / 128.0); // Minimum step size (for variable step size control)
    Num maxDt              = static_cast<Num>(0.5);         // Maximum step size (for variable step size control)
    Num decreaseFactor     = static_cast<Num>(0.9);
    Num increaseFactor     = static_cast<Num>(2.0);
    Num localTol           = static_cast<Num>(1.0e-7);      // Local error tolerance
    Num absolute_tol       = static_cast<Num>(1.0e-10);     // Takes care of division by zero (NaN and/or inf values)
    Num localTolErrorRatio = static_cast<Num>(0.1);
    std::uint8_t order      = 4;           // Order for AB/ABM
    std::uint8_t iterations = 2;           // Iterations for ABM
    std::uint8_t maxTrial   = 25;
    bool errorEstimate   = false;      // Wether to use error estimation (or not)
    bool variableSteps   = false;      // Wether to use variable step size (or not)
    bool attemptsHistory = false;      // Controls wether the errors and step sizes are stored for all attempts or not.
    bool weightedError   = true;       // Calculate error based on weighted error formula
    bool normError       = false;      // Calculate error based on norm error formula (norm or the solution is used instead of maximum of the absolute value of the elements of the solution)
};

template <FPNumber Num>
struct SolverResults 
{
    Matrix<Num> solution;      // The solution matrix.
    Matrix<Num> stepsHistory;  // History of step sizes used (contains failed steps attempted).
    Matrix<Num> errorsHistory; // History of errors (contains failed steps attempted).
    Vec<Num> timePoints; // Time points where solution is computed.
    Vec<Num> errors;     // Error estimates.
    Vec<Num> stepSizes;  // Actual step sizes used.
    wVec failedTrials;   // Number of failed trials in each step.
};

// using SolverFunc = std::function<MathEngine::SolverResults(const MathEngine::SolverParameters&)>;
// using SolverFunc = std::function<SolverResults(const SolverParameters&)>;

template <FPNumber Num>
struct SparsedMatrix
{
    Vec<Vec<std::pair<std::uint64_t, Num>>> rows; // Each row: (col, value) pairs
    SparsedMatrix(std::uint64_t N) : rows(N) {}
    void add(std::uint64_t i, std::uint64_t j, Num weight) { rows[i].emplace_back(j, weight); }
};

} // End MathEngine namespace

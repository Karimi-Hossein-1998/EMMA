# Types & Conventions (`typedefs/`)

`src/MM/typedefs/header.hpp` defines the fundamental containers and callback
types used across the toolkit. The `Matrix<T>` container is the NumCpp row-major
contiguous matrix (see `src/MM/typedefs/NumCpp/`).

## Vector aliases

```cpp
template <typename T> using Vec = std::vector<T>;
using dVec  = Vec<double>;
using wVec  = Vec<size_t>;
using dcVec = Vec<std::complex<double>>;
using iVec  = Vec<int>;
using uVec  = Vec<unsigned int>;
using bVec  = Vec<bool>;
using sVec  = Vec<std::string>;
using dMatrix = Matrix<double>;
```

## Matrix

- `Matrix<T>(rows, cols, value)` constructs a row-major, contiguous matrix.
- `m[i, j]` accesses an element; `m[i]` returns a row span; `m.ptr()` the flat
  buffer; `m.Rows()`, `m.Cols()`, `m.empty()`.
- `AppendRows` / `AppendCols` / `resize` / `SetRow` / `reshape` mutate it.

## Solver callbacks

```cpp
using func     = std::function<double(double)>;                  // scalar 1D
using Func     = std::function<dVec(double)>;                    // vector 1D
using MyFunc   = std::function<void(double, const dVec&, dVec&)>; // (t, y) -> dy/dt
using CallBackFunc = std::function<void(const OneStepSolverResult&)>;

struct OneStepSolverResult { dVec sol, stepHistory, errorHistory; double timePoint, error, stepSize; size_t failedTrial; };

struct ODESolverParameters {
    MyFunc derivative; CallBackFunc onStep; dVec initialConditions;
    double t0, t1, dt, minDt, maxDt, decreaseFactor, increaseFactor;
    double localTol, absolute_tol, localTolErrorRatio;
    uint8_t order, iterations, maxTrial;
    bool errorEstimate, variableSteps, attemptsHistory, weightedError, normError;
};

struct SolverResults { dMatrix solution, stepsHistory, errorsHistory; dVec timePoints, errors, stepSizes; wVec failedTrials; };
using SolverFunc = std::function<SolverResults(const ODESolverParameters&)>;
```

## Sparse matrix

```cpp
struct SparsedMatrix { Vec<Vec<std::pair<uint64_t, double>>> rows; ... };
```

used by the sparse Kuramoto model for memory-efficient adjacency storage.

## Conventions

- `PI` (`3.14159…`) is defined as `MathEngine::PI`.
- All model/solver code lives in `namespace MathEngine`; serialisation helpers in
  `namespace MathEngine::IO`.

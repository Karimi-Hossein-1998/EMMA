# Serialisation (`utility/write.hpp`)

`src/MM/utility/write.hpp` (namespace `MathEngine::IO`) writes matrices and
vectors to disk as human-readable text or raw binary.

## API

```cpp
namespace MathEngine::IO {
enum class FPFormat   { Fixed, Scientific, Default };
enum class Alignment  { Left, Center, Right, None };

struct WriteOptions {
    std::filesystem::path path;
    std::string_view separator;   // e.g. ","
    std::string_view header, comment, footer;
    size_t colWidth; int precision;
    FPFormat format; Alignment alignment;
    bool append; bool binary;
};

template <Number T> void WriteMatrix(const Matrix<T>& mat, const WriteOptions&);
template <typename T> void WriteVector(std::span<const T> vec, const WriteOptions&);
}
```

## Behaviour

- **Text mode** — emits an optional `# header` line, then one value per cell with
  the configured column width, precision, float format (fixed/scientific) and
  alignment (using `std::format`).
- **Binary mode** — writes row/col counts followed by the raw contiguous buffer
  (fast, exact round-trip; matrix-major ordering).
- Parent directories are created automatically.

## Use cases

The GUI's Save tab serialises every model's artefacts through these functions:
solutions, time points, order parameters, coupling matrices, and — for molecular
dynamics — the final particle state and the observables time series
(temperature, energies, pressure, $\psi_4$, $\psi_6$, MSD).

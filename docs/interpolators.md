# Interpolators

`src/MM/interpolators/` provides Lagrange and Newton interpolation that return a
callable `func` (scalar) or `Func` (vector-valued) interpolating a set of points
$(x_i, y_i)$.

## Lagrange interpolation

The interpolating polynomial of degree $n-1$ through $n$ points is

$$
P(x) = \sum_{i=1}^{n} y_i\, L_i(x),
\qquad
L_i(x) = \prod_{j\neq i} \frac{x - x_j}{x_i - x_j}.
$$

`lagrange_interpolator(x_points, y_points)` builds the $L_i$ on the fly.
`barycentric_lagrange_interpolator` precomputes the barycentric weights

$$
w_i = \prod_{j\neq i}(x_i - x_j)^{-1}
$$

and evaluates via the numerically stable second barycentric form

$$
P(x) = \frac{\sum_i \frac{w_i}{x-x_i} y_i}{\sum_i \frac{w_i}{x-x_i}},
$$

which is $O(n)$ per evaluation (after the $O(n^2)$ weight precomputation) and
avoids the repeated products of the naive form.

## Newton divided differences

Newton's form is

$$
P(x) = f[x_1] + f[x_1,x_2](x-x_1) + f[x_1,x_2,x_3](x-x_1)(x-x_2) + \dots,
$$

where the divided differences are defined recursively

$$
f[x_i] = y_i,
\qquad
f[x_i,\dots,x_{i+k}] = \frac{f[x_{i+1},\dots,x_{i+k}] - f[x_i,\dots,x_{i+k-1}]}{x_{i+k} - x_i}.
$$

`newton_interpolator` computes the table and evaluates with Horner's rule in
$O(n)$. It is convenient for adding points incrementally (a new point only adds
one divided difference).

## Usage

```cpp
MathEngine::func  f = MathEngine::lagrange_interpolator(xs, ys);
MathEngine::Func  g = MathEngine::newton_interpolator(xs, yMatrix);  // vector-valued
double y = f(x);
```

## Notes

- Both interpolators reproduce the data exactly and are exact for polynomials of
  degree $< n$.
- For many points, high-degree polynomials can exhibit Runge's phenomenon;
  prefer piecewise interpolation or the barycentric form for stability.

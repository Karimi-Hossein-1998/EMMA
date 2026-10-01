# Runge–Kutta Methods (RK1–RK4 and variants)

Explicit one-step methods for the initial-value problem

$$
\frac{\mathrm{d}y}{\mathrm{d}t} = f(t,y), \qquad y(t_0) = y_0.
$$

The MM toolkit provides RK1 (Euler), RK2 (midpoint), RK3, the classical RK4, and
three RK4 variants (3/8 rule, Gill, Ralston). All follow the generic $s$-stage
form

$$
y_{n+1} = y_n + \mathrm{d}t\sum_{i=1}^{s} b_i k_i,
\qquad
k_i = f\!\left(t_n + c_i\,\mathrm{d}t,\; y_n + \mathrm{d}t\sum_{j=1}^{i-1} a_{ij}k_j\right).
$$

## Butcher tableaus

**RK1 (Euler), $s=1$:**

$$
\begin{array}{c|c} 0 & \\ \hline & 1 \end{array}
$$

**RK2 (midpoint), $s=2$:**

$$
\begin{array}{c|cc}
0 & & \\
1/2 & 1/2 & \\
\hline & 0 & 1
\end{array}
$$

**RK3 (Kutta's third-order), $s=3$:**

$$
\begin{array}{c|ccc}
0 & & & \\
1/2 & 1/2 & & \\
1 & -1 & 2 & \\
\hline & 1/6 & 2/3 & 1/6
\end{array}
$$

**RK4 (classical), $s=4$:**

$$
\begin{array}{c|cccc}
0 & & & & \\
1/2 & 1/2 & & & \\
1/2 & 0 & 1/2 & & \\
1 & 0 & 0 & 1 & \\
\hline & 1/6 & 1/3 & 1/3 & 1/6
\end{array}
$$

**RK4 3/8 rule:**

$$
\begin{array}{c|cccc}
0 & & & & \\
1/3 & 1/3 & & & \\
2/3 & -1/3 & 1 & & \\
1 & 1 & -1 & 1 & \\
\hline & 1/8 & 3/8 & 3/8 & 1/8
\end{array}
$$

**RK4 Gill:**

$$
\begin{array}{c|cccc}
0 & & & & \\
1/2 & 1/2 & & & \\
1/2 & (\sqrt2-1)/2 & (2-\sqrt2)/2 & & \\
1 & 0 & -\sqrt2/2 & (2+\sqrt2)/2 & \\
\hline & 1/6 & (2-\sqrt2)/6 & (2+\sqrt2)/6 & 1/6
\end{array}
$$

**RK4 Ralston:**

$$
\begin{array}{c|cccc}
0 & & & & \\
0.4 & 0.4 & & & \\
0.45573725 & 0.29697761 & 0.15875964 & & \\
1 & 0.21810040 & -3.05096516 & 3.83286476 & \\
\hline & 0.17476028 & -0.55148066 & 1.20553560 & 0.17118478
\end{array}
$$

## Order conditions

A Runge–Kutta method has order $p$ if its coefficients satisfy the *Butcher order
conditions*, obtained by matching the Taylor expansion of the numerical solution
to that of the exact solution. With the row sums $c_i=\sum_j a_{ij}$, the
conditions up to order 4 are:

| Order | Conditions |
|---|---|
| 1 | $\sum_i b_i = 1$ |
| 2 | $\sum_i b_i c_i = \tfrac12$ |
| 3 | $\sum_i b_i c_i^2 = \tfrac13$, $\quad \sum_{i,j} b_i a_{ij} c_j = \tfrac16$ |
| 4 | $\sum_i b_i c_i^3 = \tfrac14$, $\quad \sum_{i,j} b_i c_i a_{ij} c_j = \tfrac18$, $\quad \sum_{i,j} b_i a_{ij} c_j^2 = \tfrac1{12}$, $\quad \sum_{i,j,k} b_i a_{ij} a_{jk} c_k = \tfrac1{24}$ |

(The conditions form a tree-ordered set; the $4^{\text{th}}$-order one
$\sum b_i a_{ij} a_{jk} c_k = \tfrac1{24}$ is the "binary tree" condition.) The
classical RK4 tableau can be checked against these: $c=(0,\tfrac12,\tfrac12,1)$,
$b=(\tfrac16,\tfrac13,\tfrac13,\tfrac16)$ satisfy them, which is why RK4 is
fourth-order.

## Order, error and stability

A method of order $p$ has local truncation error $\mathcal{O}(\mathrm{d}t^{p+1})$
and global error $\mathcal{O}(\mathrm{d}t^p)$. Applied to the linear test equation
$y'=\lambda y$ with $z=\lambda\,\mathrm{d}t$, each method has an amplification
factor equal to the degree-$p$ Taylor truncation of $e^{z}$:

$$
R(z) = 1 + z + \frac{z^2}{2!} + \dots + \frac{z^p}{p!}.
$$

Stability requires $|R(z)|<1$. All these methods are explicit and therefore only
conditionally stable; for stiff problems use smaller $\mathrm{d}t$ or an implicit
method.

- RK1: $p=1$ (Euler).
- RK2 midpoint: $p=2$.
- RK3: $p=3$.
- RK4 family: $p=4$. The variants (3/8, Gill, Ralston) share the classical RK4
  order but differ in stage structure — Gill minimises round-off error, Ralston
  minimises the leading error bound, and 3/8 uses symmetric stages.

## Usage

All solvers implement the `SolverFunc` interface

```cpp
using SolverFunc = std::function<SolverResults(const ODESolverParameters&)>;
```

and are selected in the GUI's Solver tab (`RK1`…`RK4`, `RK4_38`, `RK4_Gill`,
`RK4_Ralston`). Implementations live in `src/MM/solvers/ODE/rk/explicit/`.

## References

- Butcher, *Numerical Methods for Ordinary Differential Equations* (2016).
- Hairer, Nørsett & Wanner, *Solving Ordinary Differential Equations I* (1993).

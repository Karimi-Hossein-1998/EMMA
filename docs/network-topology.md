# Network Topologies

`src/MM/network/topology.hpp` generates adjacency (coupling) matrices for the
Kuramoto model. The GUI's Topology tab exposes these through `NetworkTopology`.

| Enum | Topology | Weights |
|---|---|---|
| `Uniform` / `UniformSymmetric` | dense random | $[w_{\min}, w_{\max}]$ |
| `ErdosRenyi` (+ `Uniform`, `Symmetric`, `SymmetricUniform`) | $G(N,p)$ | $[w_{\min}, w_{\max}]$ |
| `SmallWorld` / `SmallWorldDirected` | Watts–Strogatz | ring + rewiring |
| `Modular` | block structure | dense intra / sparse inter |
| `Hierarchical` | nested modules | level-dependent probabilities |

## Mathematical descriptions

**Erdős–Rényi $G(N,p)$** — each directed edge exists independently with
probability $p$; weights are uniform in $[w_{\min}, w_{\max}]$. The `Uniform`
variant fixes the edge count to exactly $N(N-1)p$.

**Watts–Strogatz small-world** — start from a ring lattice and rewire each edge
with probability $p$. Produces high clustering with short average path length
(small-world behaviour).

**Modular** — $N$ nodes are split into modules; intra-module edges are dense,
inter-module edges sparse, with distinct probabilities $p_{\mathrm{in}}$,
$p_{\mathrm{out}}$ and weight ranges.

**Hierarchical** — modules are nested recursively over $L$ levels with
level-dependent connection probabilities and a decay ratio.

## Utilities

- `dense_to_sparse(adj)` / `dense_to_sparse_conditional(...)` — convert a dense
  matrix to a `SparsedMatrix` when the density is low enough.
- `density(adj, threshold)` — fraction of non-zero off-diagonal entries.
- `multilayered(layers)` / `effective_multiplex(...)` — multiplex aggregation.

## Use cases

- Erdős–Rényi and small-world networks reproduce the standard synchronisation
  and chimera-state experiments.
- Modular/hierarchical topologies model community structure; the GUI's
  per-module order-parameter plot reveals partial synchrony between communities.

## References

- Erdős & Rényi (1959), *Publ. Math. Debrecen* **6**, 290.
- Watts & Strogatz (1998), *Nature* **393**, 440.

# MeanderChaos

Code and data for: Noh, B. & Wani, O. (2026). Cutoffs as a sufficient condition for chaos in kinematic river channel evolution. *Communications Earth & Environment*.

## The variable-dimensionality problem

The Howard-Knutson kinematic model represents a river as a Lagrangian centerline of $n$ nodes that migrate normal to the local tangent at a rate determined by curvature and an upstream-weighted convolution. After each migration step, the centerline is resampled to uniform spacing:

```python
unew = np.linspace(0, 1, 1 + int(round(s[-1] / deltas)))
```

The `round()` means an infinitesimal perturbation in arc length can change $n$ by one. In Lagrangian coordinates, this shifts every node index downstream, injecting a spurious $O(1)$ divergence that has nothing to do with geometry. Any chaos diagnostic built on node-to-node comparison is contaminated by this artifact.

## Eulerian state representation

We sidestep variable dimensionality by projecting the Lagrangian centerline onto a fixed Eulerian binary grid. The state becomes a matrix $S_{k\ell}(t) \in \{0,1\}$ of constant dimension, independent of how many nodes parameterize the curve:

```python
def rasterize_channel(ch):
    g = np.zeros((rows, cols), dtype=bool)
    xs = np.linspace(ch.x[:-1], ch.x[1:], 10)
    ys = np.linspace(ch.y[:-1], ch.y[1:], 10)
    col_idx = ((xs - xmin) / cell_size).astype(int)
    row_idx = ((ys - ymin) / cell_size).astype(int)
    mask = (0 <= col_idx) & (col_idx < cols) & (0 <= row_idx) & (row_idx < rows)
    g[row_idx[mask], col_idx[mask]] = True
    return g
```

A curve with 1000 nodes and the same curve with 1001 nodes occupy identical grid cells. The Hamming distance $d_H(t) = \|S^*(t) - S(t)\|_1$ between reference and perturbed occupancy fields measures purely geometric divergence.

<p align="center">
  <img src="figures/fig1.png" width="95%">
</p>

Fig. 1. (a) Lagrangian ensemble of 100 realizations from near-identical initial conditions. (b-d) The same state projected onto Eulerian grids at 10 m, 50 m, and 100 m resolution.

## The cutoff algorithm is deterministic

The reconnection logic in `meanderpy` contains no stochastic elements:

1. Pairwise Euclidean distances are computed exactly via `scipy.spatial.distance.cdist`. The first pair satisfying $\|x_i - x_j\| < d_c$ in row-major order is selected deterministically.

2. Topological surgery is an exact array splice: `x = np.hstack((x[:i+1], x[j:]))`. No floating-point noise is introduced.

3. Resampling uses `scipy.interpolate.splprep` with `s=0`, forcing the cubic B-spline through every remaining node without smoothing.

Given identical inputs, the algorithm produces identical outputs. The divergence is not a numerical artifact.

## Mechanism: event-driven chaos

The system is a hybrid continuous-discrete dynamical system. The continuous phase (curvature-driven migration) is smooth and non-chaotic. Chaos enters through the discrete phase (cutoffs) via the following cascade:

1. Two trajectories separated by $10^{-5}$ m approach a cutoff threshold. Because they differ microscopically, one crosses $\|x_i - x_j\| < d_c$ at timestep $N$; the other misses by a fraction of a millimeter and waits until $N+1$.

2. During that one-step delay, the trajectories evolve under different topologies. The first has dropped its oxbow loop; its downstream nodes integrate over a short, straight neck. The second still retains the loop; its downstream nodes integrate over high-curvature geometry via the upstream convolution $\sum R_0[i::-1] \cdot G$.

3. For that single timestep, the downstream migration rates differ by $O(1)$, injecting a meter-scale geometric separation.

4. Once separated by meters, the next cutoff threshold is straddled at an even wider time offset, and the process cascades.

This is the same threshold-reset mechanism that produces chaos in impact oscillators and other hybrid systems. The continuous dynamics provide the stretching; the discrete events provide the folding.

## Counterfactual experiment

The test is a single binary switch: cutoffs on or off. With cutoffs disabled, the two trajectories remain coincident indefinitely ($d_H = 0$). With cutoffs enabled, $\ln d_H$ grows linearly, yielding a positive finite-time Lyapunov exponent:

$$\lambda_{\mathrm{FT}} = \frac{1}{t_2 - t_1} \ln\!\left(\frac{d_H(t_2)}{d_H(t_1)}\right)$$

```python
def get_log_diff(ch1, ch2, cell_size):
    g1, g2 = raster(ch1), raster(ch2)
    diff = np.count_nonzero(g1 != g2)
    return np.log(diff) if diff > 0 else np.nan

log_norms = [get_log_diff(chb_ref.channels[t], chb_pert.channels[t], 50.0)
             for t in range(NIT)]
```

<p align="center">
  <img src="figures/hamming.png" width="80%">
</p>

The growth rate converges with grid refinement, is insensitive to perturbation magnitude over ten orders of magnitude ($10^{0}$ to $10^{-10}$ m), and persists across Kinoshita planforms with $\theta_0 \in \{0.5, 1.0, 1.5, 2.0\}$. The Lyapunov exponent scales with migration rate $k_\ell$ but is invariant to the cutoff threshold $d_c$, while $d_c$ controls the frequency of topological resets. The predictability horizon, defined as the number of cutoffs per Lyapunov time, saturates at approximately 10 events in the neck-cutoff regime.

## Repository structure

```
MeanderChaos/
├── MeanderChaos_Tutorial.ipynb           # Executable version of the above
├── scripts/
│   ├── MeanderChaos_Benettin.py          # Benettin algorithm for Lyapunov exponents
│   ├── Gridsize_and_Perturbation_Test.py # Resolution and perturbation sensitivity
│   └── Recurrence_Plot.py               # Recurrence quantification analysis
└── figures/
    └── codes/                            # Paper figure reproduction scripts
```

## Interactive demo

[braydennoh.github.io/chaotic-rivers](https://braydennoh.github.io/chaotic-rivers.html)

## Citation

```bibtex
@article{noh2026cutoffs,
  title     = {Cutoffs as a sufficient condition for chaos in kinematic
               river channel evolution},
  author    = {Noh, Brayden and Wani, Omar},
  journal   = {Communications Earth \& Environment},
  year      = {2026},
  publisher = {Nature Publishing Group}
}
```

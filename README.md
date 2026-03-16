# MeanderChaos

**Cutoffs as a sufficient condition for chaos in kinematic river channel evolution**

Noh, B. & Wani, O. (2026). *Communications Earth & Environment*.

## Overview

Rivers shape their floodplains through meander growth and cutoffs, which reorganize channel geometry. We test whether cutoffs alone are sufficient to generate deterministic chaos using a kinematic meander model formulated at fixed spatial resolution.

**Key finding:** Trajectories with cutoffs exhibit sustained exponential divergence, whereas those without cutoffs do not. The inferred Lyapunov exponent converges with grid resolution, is insensitive to perturbation magnitude, and is consistent across diverse initial planforms.

## Usage

This repository uses [`meanderpy`](https://github.com/zsylvester/meanderpy) to simulate river planform evolution from infinitesimally perturbed initial conditions. Each evolving centerline is rasterized onto a fixed Eulerian grid, enabling the computation of Hamming distance between binary occupancy fields as a measure of divergence.

### Requirements

```bash
pip install numpy matplotlib scipy cmocean meanderpy
```

### Quick Start

```bash
jupyter notebook MeanderChaos_Tutorial.ipynb
```

---

### 1. Simulation Setup

We run two deterministic simulations from near-identical initial conditions — a reference and a single-node perturbation ($\delta = 10^{-5}$ m):

```python
SECONDS_PER_YEAR = 365.25 * 24 * 3600
NIT = 1001          # Iterations
W = 100.0           # Channel width (m)
MAG = 1e-5          # Perturbation magnitude (m)
CRDIST = 2 * W      # Cutoff threshold

KL_M_PER_YR = 100.0
DT_YEARS = 0.1

kl = KL_M_PER_YR / SECONDS_PER_YEAR
dt = DT_YEARS * SECONDS_PER_YEAR

chb_ref  = run_sim(0.0, CRDIST)   # Reference
chb_pert = run_sim(MAG, CRDIST)   # Perturbed by 1e-5 m
```

---

### 2. Lagrangian → Eulerian Transformation

Two completely different planforms can have the same sinuosity — so tracking global parameters cannot detect chaos. Instead, we rasterize each Lagrangian centerline onto a fixed Eulerian binary grid, giving a fixed-dimensional state for comparison:

```python
def rasterize_channel(ch):
    g = np.zeros((rows, cols), dtype=bool)

    # Densify polyline segments so thin lines don't skip cells
    xs = np.linspace(ch.x[:-1], ch.x[1:], 10)
    ys = np.linspace(ch.y[:-1], ch.y[1:], 10)

    col_idx = ((xs - xmin) / cell_size).astype(int)
    row_idx = ((ys - ymin) / cell_size).astype(int)

    mask = (0 <= col_idx) & (col_idx < cols) & (0 <= row_idx) & (row_idx < rows)
    g[row_idx[mask], col_idx[mask]] = True
    return g

G1 = rasterize_channel(chb_ref.channels[t_idx])   # Reference occupancy
G2 = rasterize_channel(chb_pert.channels[t_idx])   # Perturbed occupancy
```

<p align="center">
  <img src="figures/fig1.png" width="95%" alt="Lagrangian and Eulerian representations">
</p>

> **Fig. 1.** Lagrangian ensemble of 100 realizations (a), and Eulerian representation on grids with resolutions of 10 m (b), 50 m (c), and 100 m (d).

---

### 3. Hamming Distance

Divergence is measured via the Hamming distance between occupancy fields:

$$d_H(t) = \sum_{i,j} \left| G_{\text{ref}}(t)_{i,j} - G_{\text{pert}}(t)_{i,j} \right|$$

```python
def get_log_diff(ch1, ch2, cell_size):
    g1 = raster(ch1)
    g2 = raster(ch2)
    diff = np.count_nonzero(g1 != g2)
    return np.log(diff) if diff > 0 else np.nan

log_norms = np.array([
    get_log_diff(chb_ref.channels[t], chb_pert.channels[t], cell_size=50.0)
    for t in range(0, NIT, 1)
])
```

A linear trend in $\log(d_H)$ indicates exponential growth of the initial perturbation — the signature of deterministic chaos.

<p align="center">
  <img src="figures/hamming.png" width="80%" alt="Hamming distance divergence">
</p>

The finite-time Lyapunov exponent is estimated from the linear growth window:

$$\lambda_{\mathrm{FT}}(t_1,t_2) = \frac{1}{t_2 - t_1} \ln\!\left(\frac{d_H(t_2)}{d_H(t_1)}\right)$$

---

## Repository Structure

```
MeanderChaos/
├── README.md
├── MeanderChaos_Tutorial.ipynb        # Interactive tutorial (start here)
├── scripts/
│   ├── MeanderChaos_Benettin.py       # Benettin algorithm for Lyapunov exponents
│   ├── Gridsize_and_Perturbation_Test.py   # Grid resolution & perturbation sensitivity
│   └── Recurrence_Plot.py             # Recurrence quantification analysis
└── figures/
    ├── codes/                         # Scripts to reproduce paper figures
    │   ├── figure1code.py
    │   ├── figure2code.py
    │   ├── figure3code.py
    │   ├── figure4code.py
    │   └── figure5code.py
    ├── fig1.png
    └── hamming.png
```

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/MeanderChaos_Benettin.py` | Benettin algorithm for Lyapunov exponents on the Eulerian grid |
| `scripts/Gridsize_and_Perturbation_Test.py` | Sensitivity analysis: grid cell size and perturbation magnitude |
| `scripts/Recurrence_Plot.py` | Recurrence quantification analysis of planform evolution |
| `figures/codes/figure[1-5]code.py` | Reproduction scripts for each paper figure |

## Interactive Demo

**[braydennoh.github.io/chaotic-rivers](https://braydennoh.github.io/chaotic-rivers.html)**

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

## License

MIT

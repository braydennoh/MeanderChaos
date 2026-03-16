# MeanderChaos

**Cutoffs as a sufficient condition for chaos in kinematic river channel evolution**

Noh, B. & Wani, O. (2025). *Communications Earth & Environment*.

<p align="center">
  <img src="figures/ensemble.png" width="90%" alt="Ensemble of diverging river planforms">
</p>

> *Ensemble of 100 simulations initialized from near-identical conditions. With cutoffs enabled, planforms diverge exponentially.*

## The Problem

Lowland rivers evolve through gradual curvature-driven migration punctuated by abrupt cutoff events that reorganize planform topology. We ask: **are cutoffs alone sufficient to make river evolution chaotic?**

Two completely different planforms can have the same sinuosity -- so tracking global parameters like sinuosity cannot detect chaos. Instead, we map each centerline onto a **fixed Eulerian grid** of binary cells (channel or floodplain), giving a fixed-dimensional state for comparison.

## Key Result

**When cutoffs are disabled, two nearly identical rivers stay identical forever. When cutoffs are enabled, they diverge exponentially.**

<p align="center">
  <img src="figures/lagrangian.png" width="90%" alt="Lagrangian centerline evolution">
</p>

The finite-time Lyapunov exponent converges with grid refinement, remains independent of initial perturbation size, and appears consistent across different river shapes.

## Method

### 1. Paired Simulations

Two simulations are initialized from identical sine-generated curves, differing by a single-node perturbation of magnitude $\delta = 10^{-5}$ m. Both are evolved using the deterministic Howard--Knutson (1984) curvature-driven migration model via [`meanderpy`](https://github.com/zsylvester/meanderpy).

### 2. Eulerian Rasterization

Since Lagrangian node positions drift independently, direct coordinate comparison is meaningless. Each centerline is rasterized onto a fixed binary occupancy grid:

$$G_{i,j}(t) = \begin{cases} 1 & \text{if channel occupies cell } (i,j) \text{ at time } t \\\\ 0 & \text{otherwise} \end{cases}$$

<p align="center">
  <img src="figures/eulerian.png" width="80%" alt="Eulerian grid comparison">
</p>

> *Purple: agreement between reference and perturbed runs. Red/blue: spatial divergence.*

### 3. Hamming Distance

Divergence is measured via the Hamming distance between occupancy fields:

$$d_H(t) = \sum_{i,j} \left| G_{\text{ref}}(t)_{i,j} - G_{\text{pert}}(t)_{i,j} \right|$$

A linear trend in $\log(d_H)$ indicates exponential growth of the initial perturbation.

<p align="center">
  <img src="figures/hamming.png" width="80%" alt="Hamming distance divergence">
</p>

### 4. Lyapunov Exponent

The maximum Lyapunov exponent is estimated via the Benettin algorithm on the Eulerian grid, confirming positive exponents **only** when cutoffs are enabled.

<p align="center">
  <img src="figures/PNAS_Supp.gif" width="80%" alt="Supplementary animation">
</p>

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
    ├── ensemble.png
    ├── lagrangian.png
    ├── eulerian.png
    ├── hamming.png
    ├── meander_evolution.gif
    ├── PNAS_Supp.gif
    └── PNAS_Supp.mp4
```

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/MeanderChaos_Benettin.py` | Benettin algorithm for Lyapunov exponents on the Eulerian grid |
| `scripts/Gridsize_and_Perturbation_Test.py` | Sensitivity analysis: grid cell size and perturbation magnitude |
| `scripts/Recurrence_Plot.py` | Recurrence quantification analysis of meander planform evolution |
| `figures/codes/figure[1-5]code.py` | Reproduction scripts for each paper figure |

## Getting Started

### Requirements

```
numpy
matplotlib
scipy
meanderpy
cmocean
```

```bash
pip install numpy matplotlib scipy cmocean meanderpy
```

### Quick Start

Open the interactive tutorial:

```bash
jupyter notebook MeanderChaos_Tutorial.ipynb
```

Or run the Lyapunov exponent computation directly:

```bash
python scripts/MeanderChaos_Benettin.py
```

## Interactive Demo

An interactive Eulerian occupancy heatmap with adjustable grid resolution is available at:

**[braydennoh.github.io/chaotic-rivers](https://braydennoh.github.io/chaotic-rivers.html)**

## Citation

```bibtex
@article{noh2025cutoffs,
  title={Cutoffs as a sufficient condition for chaos in kinematic river channel evolution},
  author={Noh, Brayden and Wani, Omar},
  journal={Communications Earth \& Environment},
  year={2025}
}
```

## License

MIT

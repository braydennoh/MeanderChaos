# MeanderChaos

Deterministic chaos in curvature-driven river meander models. Code and data for Noh et al. (in review).

<p align="center">
  <img src="figures/meander_evolution.gif" width="80%" alt="Meander evolution">
</p>

## Overview

River meanders evolve through smooth curvature-driven migration punctuated by abrupt cutoff events that reorganize planform topology. Within a fully deterministic kinematic model, we show that cutoff events alone are sufficient to generate sensitive dependence on initial conditions -- sustained exponential divergence from infinitesimal perturbations -- when all other sources of stochasticity are excluded.

We quantify this divergence by rasterizing paired simulations onto a fixed Eulerian grid and computing the Hamming distance between binary occupancy fields over time.

## Repository Structure

```
MeanderChaos/
├── README.md
├── MeanderChaos_Tutorial.ipynb      # Interactive tutorial (start here)
├── scripts/
│   ├── MeanderChaos_Benettin.py     # Benettin algorithm for Lyapunov exponents
│   ├── Gridsize_and_Perturbation_Test.py  # Grid resolution & perturbation sensitivity
│   └── Recurrence_Plot.py           # Recurrence quantification analysis
├── figures/
│   ├── codes/                       # Scripts to reproduce paper figures
│   │   ├── figure1code.py
│   │   ├── figure2code.py
│   │   ├── figure3code.py
│   │   ├── figure4code.py
│   │   └── figure5code.py
│   ├── lagrangian.png
│   ├── eulerian.png
│   ├── hamming.png
│   ├── meander_evolution.gif
│   ├── PNAS_Supp.gif
│   └── PNAS_Supp.mp4
```

## Getting Started

### Requirements

```
numpy
matplotlib
scipy
meanderpy
cmocean
```

Install dependencies:

```bash
pip install numpy matplotlib scipy cmocean
pip install meanderpy
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

## Method

### 1. Paired Simulations

Two simulations are initialized from identical sine-generated curves, differing by a single-node perturbation of magnitude $\delta = 10^{-5}$ m:

```python
# Reference run
chb_ref = run_sim(pert=0.0, crdist=2*W)

# Perturbed run (one node shifted by 1e-5 m)
chb_pert = run_sim(pert=1e-5, crdist=2*W)
```

### 2. Eulerian Rasterization

Since Lagrangian node positions drift independently, direct coordinate comparison is meaningless. Instead, each centerline is rasterized onto a fixed binary occupancy grid:

$$G_{i,j}(t) = \begin{cases} 1 & \text{if channel occupies cell } (i,j) \text{ at time } t \\ 0 & \text{otherwise} \end{cases}$$

<p align="center">
  <img src="figures/eulerian.png" width="80%" alt="Eulerian grid comparison">
</p>

### 3. Hamming Distance

Divergence is measured via the Hamming distance between the two occupancy fields:

$$d_H(t) = \sum_{i,j} \left| G_{\text{ref}}(t)_{i,j} - G_{\text{pert}}(t)_{i,j} \right|$$

A linear trend in $\log(d_H)$ indicates exponential growth of the initial perturbation.

<p align="center">
  <img src="figures/hamming.png" width="80%" alt="Hamming distance divergence">
</p>

### 4. Lyapunov Exponent

The maximum Lyapunov exponent is estimated via the Benettin algorithm (`scripts/MeanderChaos_Benettin.py`), confirming positive exponents only when cutoffs are enabled.

## Key Result

By comparing paired simulations with cutoffs enabled vs. disabled, we show that **only cutoff-enabled runs develop sustained exponential divergence**. Cutoffs act as the topological mechanism that amplifies infinitesimal perturbations into macroscopic planform differences.

<p align="center">
  <img src="figures/PNAS_Supp.gif" width="80%" alt="Supplementary animation">
</p>

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/MeanderChaos_Benettin.py` | Benettin algorithm for computing Lyapunov exponents on the Eulerian grid |
| `scripts/Gridsize_and_Perturbation_Test.py` | Sensitivity analysis: grid cell size and perturbation magnitude |
| `scripts/Recurrence_Plot.py` | Recurrence quantification analysis of meander planform evolution |
| `figures/codes/figure[1-5]code.py` | Reproduction scripts for each paper figure |

## Citation

If you use this code, please cite:

```
Noh, B., et al. (in review). [Title]. [Journal].
```

## License

MIT

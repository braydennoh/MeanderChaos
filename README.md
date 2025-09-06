# MeanderChaos

![PNAS Supplementary Animation](https://github.com/braydennoh/MeanderChaos/blob/main/Supplement/PNAS_Supp.gif)

## Overview
River meanders evolve smoothly most of the time, but **neck cutoffs** abruptly change planform topology. In this repo we test a simple question inside a deterministic geometric model: **Are cutoffs alone sufficient to produce sensitive dependence on initial conditions (chaos)?**  
By comparing paired simulations with identical settings—once with cutoffs **enabled** and once with cutoffs **disabled**—we show that only the cutoff-enabled runs develop sustained exponential divergence on a fixed-dimension Eulerian grid.

## Key Result
With all else held constant, **enabling cutoffs yields positive finite-time Lyapunov exponents (FTLEs)**; disabling cutoffs eliminates measurable growth. The growth rate:
- is **independent of perturbation size** once the perturbation is small relative to the grid,  
- **recurs across distinct initial planforms**, and  
- **converges under grid refinement**, indicating it reflects model dynamics rather than discretization.

## Method in Brief
- **Lagrangian model**: curvature-driven lateral migration with upstream influence (via an exponential kernel), implemented with `meanderpy`.  
- **Cutoff rule**: excise a loop when the minimum distance between non-adjacent segments falls below a threshold \(d_c\).  
- **Eulerian embedding**: rasterize each centerline to a fixed \(N\times N\) binary image; measure separation  
  \( d(t)=\|\mathbf{S}^\ast(t)-\mathbf{S}(t)\|_2 \).  
- **FTLE estimate**: fit the slope of \(\ln d(t)\) over the linear window before saturation on the fixed grid.

## Reproducing the Main Experiments
```bash
# clone
git clone https://github.com/braydennoh/MeanderChaos.git
cd MeanderChaos

# create environment (example)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

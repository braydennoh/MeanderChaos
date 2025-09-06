## MeanderChaos

River meanders evolve smoothly most of the time, but neck cutoffs abruptly change planform topology. We test a simple question inside a deterministic geometric model: are cutoffs alone sufficient to produce sensitive dependence on initial conditions (chaos)?

![PNAS Supplementary Animation](https://github.com/braydennoh/MeanderChaos/blob/main/Supplement/PNAS_Supp.gif)

By comparing paired simulations with identical settings—once with cutoffs enabled and once with cutoffs disabled—we show that only the cutoff-enabled runs develop sustained exponential divergence on a fixed-dimension Eulerian grid. Just add the equation for the Lyapunov on the Eulerian grid.


### Lyapunov on a fixed eulerian grid

Let each planform at time t be rasterized to a binary occupancy field  
$\mathbf{S}(t) \in \{0,1\}^{N_x \times N_y}$, and let $\mathbf{S}^\ast(t)$ be the paired run.
The separation is

$$
d(t) = \left\| \mathbf{S}^\ast(t) - \mathbf{S}(t) \right\|_2
$$

The finite-time Lyapunov exponent over a linear-growth window \([t_1,t_2]\) is

$$
\lambda_{\mathrm{FT}}(t_1,t_2) =
\frac{1}{t_2 - t_1}
\ln\left(\frac{d(t_2)}{d(t_1)}\right),
\quad t_0 \leq t_1 < t_2
$$

## Tutorial Code Overview

The full workflow is implemented in [TutorialCode.ipynb](TutorialCode.ipynb).  
Below is a high-level description of the main steps:

1. **Imports and setup**  
   Load numerical, plotting, and meanderpy libraries. Configure figure styles and helper functions.

2. **Kinoshita initialization**  
   - Define the heading angle θ(s) and integrate to construct initial sinusoidal/Kinoshita centerlines.  
   - Apply a small transverse perturbation at a chosen node to generate paired runs (reference vs. perturbed).

3. **Resampling and rasterization tools**  
   - Densify polylines to ensure resolution consistency.  
   - Rasterize centerlines onto a fixed Eulerian grid with spacing Δ.  
   - Overlay binary occupancy fields to highlight overlap and divergence.

4. **Simulation with meanderpy**  
   - Run both channels forward using curvature-driven migration with lateral migration coefficient \(k_\ell\).  
   - Toggle the cutoff threshold \(d_c\) to enable or disable neck cutoffs.  
   - Store centerline geometry at each timestep.

5. **Visualization**  
   - Plot Lagrangian centerlines for reference and perturbed runs at selected timesteps.  
   - Plot Eulerian occupancy overlays (red = reference only, blue = perturbed only, purple = overlap).  

6. **Separation metrics**  
   - Compute the binary-grid separation \(d(t)\) between channels through time.  
   - Plot the log separation \(\ln d(t)\) to estimate finite-time Lyapunov exponents.  

The notebook contains all functions and plotting routines for reproducing the figures and analyses in the manuscript.

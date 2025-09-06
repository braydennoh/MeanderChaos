## MeanderChaos

River meanders evolve smoothly most of the time, but neck cutoffs abruptly change planform topology. We test a simple question inside a deterministic geometric model: are cutoffs alone sufficient to produce sensitive dependence on initial conditions (chaos)?

![PNAS Supplementary Animation](https://github.com/braydennoh/MeanderChaos/blob/main/Supplement/PNAS_Supp.gif)

By comparing paired simulations with identical settings—once with cutoffs enabled and once with cutoffs disabled—we show that only the cutoff-enabled runs develop sustained exponential divergence on a fixed-dimension Eulerian grid. Just add the equation for the Lyapunov on the Eulerian grid.


### lyapunov on a fixed eulerian grid

Let each planform at time t be rasterized to a binary occupancy field  
\(\mathbf{S}(t) \in \{0,1\}^{N_x \times N_y}\), and let \(\mathbf{S}^\ast(t)\) be the paired run.  
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

#!/usr/bin/env python
# coding: utf-8



import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import meanderpy as mp
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 500




SECONDS_PER_YEAR = 365.25 * 24 * 3600
NIT = 1001
W = 100.0
MAG = 1e-10
EXTENT = (0, 1000, -2500, 2500)

KL_M_PER_YR = 100.0       
CRDIST = 5*W


KV_M_PER_YR = 3.16e-5    
DT_YEARS = 0.1         

kl = KL_M_PER_YR / SECONDS_PER_YEAR
kv = KV_M_PER_YR / SECONDS_PER_YEAR
dt = DT_YEARS * SECONDS_PER_YEAR

def run_sim(pert, crdist):
    theta0 = 1.5
    n_nodes = 1000
    total_length = 10000.0
    lamb = 500.0
    base_angle = 0.0

    s = np.linspace(0.0, total_length, n_nodes)
    theta = theta0 * np.sin(2.0 * np.pi * s / lamb) + base_angle

    x = np.zeros_like(s)
    y = np.zeros_like(s)
    for i in range(1, n_nodes):
        ds = s[i] - s[i-1]
        x[i] = x[i-1] + ds * np.cos(theta[i])
        y[i] = y[i-1] + ds * np.sin(theta[i])

    if pert != 0.0:
        mid = len(x) // 2
        y[mid] += pert

    z = np.zeros_like(x)
    ch = mp.Channel(x, y, z, W, 1.0)
    chb = mp.ChannelBelt([ch], [], [0.0], [])

    saved_ts = 1
    deltas = 50.0
    pad = 0
    depths = np.ones(NIT)
    Cfs = 0.0065 * np.ones(NIT)
    dens = 1000
    t1 = t2 = t3 = 0
    aggr_factor = 0.0

    chb.migrate(
        NIT, saved_ts, deltas, pad, crdist,
        depths, Cfs,
        kl, kv, dt,
        dens, t1, t2, t3, aggr_factor
    )
    return chb

chb_ref  = run_sim(0.0, CRDIST)
chb_pert = run_sim(MAG, CRDIST)




import numpy as np
import matplotlib.pyplot as plt

t_idx = 400
cell_size = 10.0

ch = chb_ref.channels[t_idx]

# --------- rasterize this single channel ----------
x = ch.x
y = ch.y

xmin = np.floor(x.min() / cell_size) * cell_size
xmax = np.ceil(x.max() / cell_size) * cell_size
ymin = np.floor(y.min() / cell_size) * cell_size
ymax = np.ceil(y.max() / cell_size) * cell_size

cols = int((xmax - xmin) / cell_size)
rows = int((ymax - ymin) / cell_size)

g = np.zeros((rows, cols), dtype=int)

xs = np.linspace(x[:-1], x[1:], 10)
ys = np.linspace(y[:-1], y[1:], 10)

col_idx = ((xs - xmin) / cell_size).astype(int)
row_idx = ((ys - ymin) / cell_size).astype(int)

mask = (0 <= col_idx) & (col_idx < cols) & (0 <= row_idx) & (row_idx < rows)
g[row_idx[mask], col_idx[mask]] = 1

# --------- PLOT CENTERLINE + RASTER GRID ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4))

# left: centerline
ax[0].plot(ch.x, ch.y, lw=2, color='red')
ax[0].set_title(f"Centerline at iteration {t_idx}")
ax[0].set_aspect('equal')

# right: raster grid
ax[1].imshow(g, cmap="viridis", origin="lower")
ax[1].set_title("Rasterized grid (cell size = 10 m)")

plt.tight_layout()
plt.show()
0




def get_log_diff(ch1, ch2, cell_size):
    if ch1 is None or ch2 is None:
        return np.nan

    x_all = np.concatenate([ch1.x, ch2.x])
    y_all = np.concatenate([ch1.y, ch2.y])

    xmin = x_all.min()
    xmax = x_all.max()
    ymin = y_all.min()
    ymax = y_all.max()

    xmin_snap = np.floor(xmin / cell_size) * cell_size
    xmax_snap = np.ceil(xmax / cell_size) * cell_size
    ymin_snap = np.floor(ymin / cell_size) * cell_size
    ymax_snap = np.ceil(ymax / cell_size) * cell_size

    cols = int((xmax_snap - xmin_snap) / cell_size)
    rows = int((ymax_snap - ymin_snap) / cell_size)

    def raster(c):
        g = np.zeros((rows, cols), dtype=bool)

        xs = np.linspace(c.x[:-1], c.x[1:], 10)
        ys = np.linspace(c.y[:-1], c.y[1:], 10)

        col_idx = ((xs - xmin_snap) / cell_size).astype(int)
        row_idx = ((ys - ymin_snap) / cell_size).astype(int)

        mask = (0 <= col_idx) & (col_idx < cols) & (0 <= row_idx) & (row_idx < rows)
        g[row_idx[mask], col_idx[mask]] = True
        return g

    g1 = raster(ch1).astype(int)
    g2 = raster(ch2).astype(int)

    return_val = np.linalg.norm(g2 - g1)
    return np.log(return_val) if return_val > 0 else np.nan

cell_size = 10.0

# sample every 100 timesteps
step = 1
time_phys = np.arange(NIT) * DT_YEARS
indices = np.arange(0, NIT, step)

log_norms = np.array([
    get_log_diff(chb_ref.channels[t], chb_pert.channels[t], cell_size=cell_size)
    for t in indices
])

mask = np.isfinite(log_norms)
t_data = time_phys[indices][mask]   # still in YEARS
y_data = log_norms[mask]
plt.plot(t_data,y_data)




cell_size = 10.0  # you can increase this to 20, 50, ... to save memory

# --------------------------------------------------
# 1) Choose which run to analyse (ONE channel belt)
# --------------------------------------------------
chb = chb_ref  

step = 10                    # every step; try 5 or 10 if memory is tight
indices = np.arange(0, NIT, step)

# --------------------------------------------------
# 2) Compute a global bounding box for this run
#    so every raster has the same grid
# --------------------------------------------------
xs_all = []
ys_all = []
for t in indices:
    ch = chb.channels[t]
    if ch is None:
        continue
    xs_all.append(ch.x)
    ys_all.append(ch.y)

x_all = np.concatenate(xs_all)
y_all = np.concatenate(ys_all)

xmin = x_all.min()
xmax = x_all.max()
ymin = y_all.min()
ymax = y_all.max()

xmin_snap = np.floor(xmin / cell_size) * cell_size
xmax_snap = np.ceil(xmax / cell_size) * cell_size
ymin_snap = np.floor(ymin / cell_size) * cell_size
ymax_snap = np.ceil(ymax / cell_size) * cell_size

cols = int((xmax_snap - xmin_snap) / cell_size)
rows = int((ymax_snap - ymin_snap) / cell_size)

print(f"Raster grid: {rows} rows × {cols} cols")

# --------------------------------------------------
# 3) Rasterizer for a single channel
# --------------------------------------------------
def rasterize_channel(c):
    g = np.zeros((rows, cols), dtype=bool)

    # interpolate between nodes along the centerline
    xs = np.linspace(c.x[:-1], c.x[1:], 10)
    ys = np.linspace(c.y[:-1], c.y[1:], 10)

    col_idx = ((xs - xmin_snap) / cell_size).astype(int)
    row_idx = ((ys - ymin_snap) / cell_size).astype(int)

    mask = (0 <= col_idx) & (col_idx < cols) & (0 <= row_idx) & (row_idx < rows)
    g[row_idx[mask], col_idx[mask]] = True
    return g

# --------------------------------------------------
# 4) Build matrix of raster states over time
# --------------------------------------------------
states = []
t_list = []

for t in indices:
    ch = chb.channels[t]
    if ch is None:
        continue
    g = rasterize_channel(ch)         # (rows, cols) binary grid
    states.append(g.ravel().astype(float))   # flatten to 1D vector
    t_list.append(t * DT_YEARS)       # physical time in years

X = np.vstack(states)   # shape: (n_times, rows*cols)
t_rp = np.array(t_list)

print("X shape:", X.shape)

def recurrence_plot(X, eps=None, percentile=10):
    X = np.asarray(X, float)
    N = X.shape[0]
    diff = X[:, None, :] - X[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)

    if eps is None:
        dvec = dist[np.triu_indices(N, k=1)]
        eps = np.percentile(dvec, percentile)

    R = (dist <= eps).astype(int)
    return R, eps

R, eps = recurrence_plot(X, percentile=10)

# --------------------------------------------------
# 6) Plot RP with time axes instead of indices
# --------------------------------------------------
plt.figure(figsize=(5, 5))
plt.imshow(R, origin='lower', cmap='binary', interpolation='nearest',
           extent=[t_rp[0], t_rp[-1], t_rp[0], t_rp[-1]])
plt.xlabel("time (years)")
plt.ylabel("time (years)")
plt.title(f"Spatial RP of chb_ref (ε={eps:.3g})")
plt.tight_layout()
plt.show()




import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import meanderpy as mp
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext, FixedLocator
from matplotlib.cm import ScalarMappable

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.dpi"] = 1000


states = []
t_list = []

for t in indices:
    ch = chb.channels[t]
    if ch is None:
        continue
    g = rasterize_channel(ch)                # (rows, cols) binary grid
    states.append(g.ravel().astype(float))   # flatten to 1D vector
    t_list.append(t * DT_YEARS)              # physical time in years

X = np.vstack(states)   # shape: (n_times, rows*cols)
t_rp = np.array(t_list)

print("X shape:", X.shape)

# --------------------------------------------------
# 5) Compute Percentile Matrix
# --------------------------------------------------
def get_distance_matrix(X):
    # Euclidean distance between all pairs of rows
    # Broadcasting: (N, 1, D) - (1, N, D) -> (N, N, D) -> norm axis 2 -> (N, N)
    diff = X[:, None, :] - X[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    return dist

def to_percentiles(D):
    """
    Convert a distance matrix D into a percentile matrix P.
    Smallest distance = 0th percentile
    Largest distance = 100th percentile
    """
    flat_D = D.flatten()
    # Compute rank (1 to N^2)
    ranks = rankdata(flat_D, method='average')
    # Normalize to 0-100
    percentiles = (ranks - 1) / (len(ranks) - 1) * 100
    return percentiles.reshape(D.shape)

D = get_distance_matrix(X)
P = to_percentiles(D)

import cmocean
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(figsize=(2.6,2.4))

im = ax.imshow(
    P,
    origin='lower',
    cmap=cmocean.cm.deep,
    interpolation='nearest',
    extent=[t_rp[0], t_rp[-1], t_rp[0], t_rp[-1]],
    vmin=0,
    vmax=50
)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.15)

cb = plt.colorbar(im, cax=cax)
cb.set_label("Relative Distance (%)")

ax.set_xlabel("Time (years)")
ax.set_ylabel("Time (years)")

ax.set_xticks([0,20,40,60,80,100])
ax.set_yticks([0,20,40,60,80,100])


plt.tight_layout()
fig.savefig("/Users/braydennoh/Downloads/fig2_chutecutoff.svg", format="svg", bbox_inches="tight")
plt.show()




import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# -------------------
# Parameters
# -------------------
theta_values = [0.5, 1, 1.5, 2]
n_nodes = 500
total_length = 10000.0
lamb = 500.0
base_angle = 0.0

# Create 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(6,5))
axes = axes.flatten()

for i, theta0 in enumerate(theta_values):
    ax = axes[i]

    # -------------------
    # Build orientation and centerline
    # -------------------
    s = np.linspace(0.0, total_length, n_nodes)
    ds = s[1] - s[0]

    theta = theta0 * np.sin(2.0 * np.pi * s / lamb) + base_angle

    # Integrate orientation to get x(s), y(s)
    x = np.cumsum(np.cos(theta)) * ds
    y = np.cumsum(np.sin(theta)) * ds

    # Recenter
    x = x - x.mean()
    y = y - y.mean()

    # -------------------
    # Curvature calculation
    # -------------------
    dx_ds = np.gradient(x, ds)
    dy_ds = np.gradient(y, ds)
    d2x_ds2 = np.gradient(dx_ds, ds)
    d2y_ds2 = np.gradient(dy_ds, ds)

    numerator = dx_ds * d2y_ds2 - dy_ds * d2x_ds2
    denominator = (dx_ds**2 + dy_ds**2)**1.5
    kappa = numerator / denominator

    # -------------------
    # Plotting
    # -------------------
    # Use fixed vmin/vmax as requested
    sc = ax.scatter(x, y, c=kappa, cmap="viridis", s=2, vmin=-0.02, vmax=0.02)

    ax.set_aspect("equal", "box")
    ax.set_title(f"$\\theta_0 = {theta0}$")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_xlim([0,600])
    ax.set_ylim([-200,200])

    # Divider logic to match colorbar height for each subplot
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label("Curvature $\kappa$ [1/m]")

plt.tight_layout()
plt.show()




import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# -------------------
# Parameters
# -------------------
theta_values = [0.5, 1, 1.5]
n_nodes = 500
total_length = 10000.0
lamb = 500.0
base_angle = 0.0

# Create 2x2 subplots
fig, axes = plt.subplots(1,3, figsize=(9, 3))
axes = axes.flatten()

for i, theta0 in enumerate(theta_values):
    ax = axes[i]

    # -------------------
    # Build orientation and centerline
    # -------------------
    s = np.linspace(0.0, total_length, n_nodes)
    ds = s[1] - s[0]

    theta = theta0 * np.sin(2.0 * np.pi * s / lamb) + base_angle

    # Integrate orientation
    x = np.cumsum(np.cos(theta)) * ds
    y = np.cumsum(np.sin(theta)) * ds

    # Recenter
    x = x - x.mean()
    y = y - y.mean()

    # -------------------
    # Curvature calculation
    # -------------------
    dx_ds = np.gradient(x, ds)
    dy_ds = np.gradient(y, ds)
    d2x_ds2 = np.gradient(dx_ds, ds)
    d2y_ds2 = np.gradient(dy_ds, ds)

    numerator = dx_ds * d2y_ds2 - dy_ds * d2x_ds2
    denominator = (dx_ds**2 + dy_ds**2)**1.5
    kappa = numerator / denominator

    # -------------------
    # Tangent and normals (exclude endpoints)
    # -------------------
    # Central-difference tangent at interior nodes (1 .. n_nodes-2)
    tx = x[2:] - x[:-2]
    ty = y[2:] - y[:-2]

    tnorm = np.sqrt(tx**2 + ty**2)
    tx = tx / tnorm
    ty = ty / tnorm

    # *** Flipped unit normals ***
    # Previously: nx = -ty, ny = tx
    # Now:       nx =  ty,  ny = -tx
    nx = ty
    ny = -tx

    # Interior coordinates and curvature aligned with these normals
    x_mid = x[1:-1]
    y_mid = y[1:-1]
    kappa_mid = kappa[1:-1]

    # -------------------
    # Curvature-weighted normal vectors
    # -------------------
    vx = kappa_mid * nx
    vy = kappa_mid * ny

    # Optional global scaling factor for visibility
    scale_factor = 100.0  # tweak if arrows are too small/large
    vx_plot = vx * scale_factor
    vy_plot = vy * scale_factor

    # -------------------
    # Plotting
    # -------------------
    sc = ax.scatter(
        x, y,
        c=kappa,
        cmap="bwr",
        s=4,
        vmin=-0.02,
        vmax=0.02, zorder=1000
    )

    # Plot curvature-weighted normals
    ax.quiver(
        x_mid, y_mid,
        vx_plot, vy_plot,
        angles="xy",
        scale_units="xy",
        scale=0.03,
        width=0.008,
        color="black",
        alpha=1, zorder=10
    )

    ax.set_aspect("equal", "box")
    ax.set_title(f"$\\theta_0 = {theta0}$")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_xlim([0, 500])
    ax.set_ylim([-200, 200])
    ax.set_xticks([0,250,500])

    # Colorbar setup
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label("Curvature $\\kappa$ [1/m]")

plt.tight_layout()
fig.savefig("/Users/braydennoh/Downloads/sinousitydata.svg", format="svg", bbox_inches="tight")
plt.show()







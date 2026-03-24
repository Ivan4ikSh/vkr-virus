import numpy as np
import matplotlib.pyplot as plt
import os

# Model parameters
dt = 0.5
Tmax = 2000.0
L = 100
R0 = 1.8
Ub = 1e-3
a = 7.0
N = 1e8
data_dir = "out"

# Read Observables from file
obs_file = os.path.join(data_dir, "observables.txt")
observables = {"speed": 0.0, "sdI": 0.0, "sdR": 0.0, "finf": 0.0}

if os.path.exists(obs_file):
    with open(obs_file, "r") as f:
        for line in f:
            if "=" in line:
                key, val = line.strip().split("=")
                observables[key] = float(val)

# Format text for headers
params_str = f"Params: R0={R0}, Ub={Ub}, a={a}, L={L}, N={N:.0e}, dt={dt}"
obs_str = f"Observables: speed={observables['speed']:.4f}, sdI={observables['sdI']:.4f}, sdR={observables['sdR']:.4f}, finf={observables['finf']:.2e}"
title_text = f"{params_str}\n{obs_str}"

# Load main data
norm = np.loadtxt(os.path.join(data_dir, "norm.txt"))
finf = np.loadtxt(os.path.join(data_dir, "finf.txt"))
slice_times = np.atleast_1d(np.loadtxt(os.path.join(data_dir, "slice_times.txt")))
time = np.arange(dt, Tmax + 1e-9, dt)
n_slices = len(slice_times)

# Figure 1: I and R Distributions
fig1, ax1 = plt.subplots(figsize=(10, 6))
x_vals = np.arange(L)

for idx, t in enumerate(slice_times):
    I_slice = np.loadtxt(os.path.join(data_dir, f"I_slice_{idx}.txt"))
    R_slice = np.loadtxt(os.path.join(data_dir, f"R_slice_{idx}.txt"))
    
    alpha = 0.3 + 0.7 * idx / (n_slices - 1) if n_slices > 1 else 1.0
    
    ax1.plot(x_vals, 10 * I_slice, '-', color='red', alpha=alpha, label='10*I' if idx == 0 else "")
    ax1.plot(x_vals, R_slice, '-', color='blue', alpha=alpha, label='R' if idx == 0 else "")

ax1.set_xlabel('Antigenic coordinate (Mutation number)')
ax1.set_ylabel('Density')

# Adjust title position and add top padding
ax1.set_title(f'Distributions at different times\n{title_text}', fontsize=12, pad=15)
ax1.grid(True, linestyle=':', alpha=0.6)

handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys())

# Reserve space for the header
fig1.tight_layout(rect=[0, 0, 1, 0.95])
file_dist = os.path.join(data_dir, "1D_distributions.png")
fig1.savefig(file_dist, dpi=150)
plt.close(fig1)

# Figure 2: finf behavior (fraction infected)
fig2, ax_finf = plt.subplots(figsize=(10, 6))

ax_finf.plot(time, finf, 'r-')
ax_finf.set_xlabel('Time')
ax_finf.set_ylabel('Fraction infected (finf)')

full_title = f'Fraction infected over time\n{title_text}'
ax_finf.set_title(full_title, fontsize=12, pad=15)
ax_finf.grid(True, linestyle=':', alpha=0.6)

fig2.tight_layout()
file_diag = os.path.join(data_dir, "1D_diagnostics.png")
fig2.savefig(file_diag, dpi=150)
plt.close(fig2)

print(f"Plots saved successfully:\n1) {file_dist}\n2) {file_diag}")
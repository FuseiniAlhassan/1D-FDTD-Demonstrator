#1D FDTD demonstrator (single-file)

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import time

# Create output directory

FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

# Physical / numerical parameters (normalized units)

c = 1.0
eps0 = 1.0
mu0 = 1.0

# Grid
Nx = 800
dx = 0.5
x = np.arange(Nx) * dx

# Time step (Courant stability: c*dt/dx <= 1)
S = 0.99
dt = S * dx / c

# Simulation time
n_steps = 1800

# Field arrays (staggered grid)
Ez = np.zeros(Nx, dtype=float)
Hy = np.zeros(Nx-1, dtype=float)

# Material: relative permittivity array for Ez grid points
eps_r = np.ones(Nx, dtype=float)

# Define a dielectric slab in the middle
slab_center = int(Nx * 0.6)
slab_width = int(Nx * 0.08)
eps_r[slab_center - slab_width//2 : slab_center + slab_width//2] = 4.0

# Mur boundary storage
Ez_left_old = 0.0
Ez_right_old = 0.0

# Source: Gaussian pulse
src_pos = int(Nx * 0.2)
t0 = 60.0
spread = 18.0
def source_time(n):
    return math.exp(-0.5 * ((n - t0) / spread)**2)

# Precompute coefficients
ce = dt / (eps0 * eps_r * dx)
ch = dt / (mu0 * dx)

# Snapshots and animation
snapshot_times = [40, 100, 220, 400, 800, 1400]
snapshots = {}
frame_interval = 4
frames_data = []

# Time loop
start_time = time.time()
print("Starting FDTD run: Nx =", Nx, "dx =", dx, "dt =", dt, "n_steps =", n_steps)

for n in range(1, n_steps + 1):
    Hy += ch * (Ez[1:] - Ez[:-1])
    Ez[1:-1] += ce[1:-1] * (Hy[1:] - Hy[:-1])
    Ez[src_pos] += source_time(n)

    # Mur boundaries
    Ez[0]  = Ez_left_old + ((S - 1)/(S + 1)) * (Ez[1] - Ez[0])
    Ez_left_old = Ez[1]
    Ez[-1] = Ez_right_old + ((S - 1)/(S + 1)) * (Ez[-2] - Ez[-1])
    Ez_right_old = Ez[-2]

    # Save snapshots
    if n in snapshot_times:
        snapshots[n] = (Ez.copy(), Hy.copy())

    # Save frames for animation
    if (n % frame_interval) == 0:
        frames_data.append((n, Ez.copy(), Hy.copy()))

elapsed = time.time() - start_time
print(f"FDTD run completed in {elapsed:.2f} s. Collected {len(snapshots)} snapshots and {len(frames_data)} frames.")

# Save snapshot figures

for tstep, (Ez_snap, Hy_snap) in snapshots.items():
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax[0].plot(x, Ez_snap, label=f"Ez (t={tstep})")
    ax[0].plot(x, eps_r - 1.0, 'k--', lw=1, label='εᵣ(x)-1 (scaled)')
    ax[0].set_ylabel("Ez")
    ax[0].legend(loc='upper right', fontsize='small')
    ax[0].grid(True)

    x_h = x[:-1] + 0.5*dx
    ax[1].plot(x_h, Hy_snap, color='tab:orange', label=f"Hy (t={tstep})")
    ax[1].set_ylabel("Hy")
    ax[1].set_xlabel("x")
    ax[1].legend(loc='upper right', fontsize='small')
    ax[1].grid(True)

    fig.suptitle(f"FDTD fields at time-step {tstep}")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    fname = os.path.join(FIG_DIR, f"snapshot_t{tstep:04d}.png")
    plt.savefig(fname, dpi=200)
    plt.close(fig)
    print("Saved snapshot:", fname)

# Reflection / transmission estimate

Ez_final = Ez.copy()
left_window = slice(int(Nx*0.05), int(Nx*0.15))
right_window = slice(int(Nx*0.75), int(Nx*0.95))
E_left = np.trapz(Ez_final[left_window]**2, x[left_window])
E_right = np.trapz(Ez_final[right_window]**2, x[right_window])
E_total = np.trapz(Ez_final**2, x)
print(f"Estimated energies (integral |E|^2): left={E_left:.6f}, right={E_right:.6f}, total={E_total:.6f}")

# Final figure with measurement windows
fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
ax[0].plot(x, Ez_final, label="Ez final")
ax[0].axvspan(left_window.start*dx, left_window.stop*dx, color='green', alpha=0.15, label='left window')
ax[0].axvspan(right_window.start*dx, right_window.stop*dx, color='red', alpha=0.15, label='right window')
ax[0].plot(x, eps_r - 1.0, 'k--', lw=1, label='εᵣ(x)-1 (scaled)')
ax[0].legend(loc='upper right', fontsize='small')
ax[0].grid(True)

x_h = x[:-1] + 0.5*dx
ax[1].plot(x_h, Hy, color='tab:orange', label="Hy final")
ax[1].set_xlabel("x")
ax[1].legend(loc='upper right', fontsize='small')
ax[1].grid(True)

fig.suptitle("Final fields and measurement windows")
fig.tight_layout(rect=[0, 0, 1, 0.96])
fname = os.path.join(FIG_DIR, "final_fields_and_windows.png")
plt.savefig(fname, dpi=200)
plt.close(fig)
print("Saved final fields figure:", fname)

# Create animation

print("Creating animation (this may take a moment)...")
fig, ax = plt.subplots(figsize=(10,4))
line_e, = ax.plot([], [], lw=1.5, label='Ez')
line_eps, = ax.plot(x, eps_r - 1.0, 'k--', lw=1.0, label='εᵣ(x)-1 (scaled)')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-1.2*np.max(np.abs(Ez)), 1.2*np.max(np.abs(Ez)))
ax.set_xlabel("x")
ax.set_ylabel("Ez")
ax.set_title("1D FDTD: Ez field propagation")
ax.legend(loc='upper right', fontsize='small')
ax.grid(True)

def init_anim():
    line_e.set_data([], [])
    return line_e, line_eps

def animate_frame(frame_tuple):
    nstep, Ez_frame, Hy_frame = frame_tuple
    line_e.set_data(x, Ez_frame)
    ax.set_title(f"1D FDTD: Ez field (time-step {nstep})")
    return line_e, line_eps

anim = animation.FuncAnimation(fig, animate_frame, frames=frames_data,
                               init_func=init_anim, blit=True, interval=40)

anim_fname = os.path.join(FIG_DIR, "wave_propagation.mp4")
anim.save(anim_fname, writer='ffmpeg', fps=25, dpi=200)
plt.close(fig)
print("Saved animation:", anim_fname)
anim.save(os.path.join(FIG_DIR, "wave_propagation.gif"), writer='pillow', fps=25, dpi=200)

# Save diagnostic energy bar chart

fig, ax = plt.subplots(figsize=(5,3))
bars = ax.bar(['Left (reflected)', 'Right (transmitted)'], [E_left, E_right], color=['green','red'])
ax.set_ylabel("Integrated |E|^2")
ax.set_title("Estimated reflected and transmitted energy")
for bar in bars:
    height = bar.get_height()
    ax.annotate(f"{height:.4f}", xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize='small')
fname = os.path.join(FIG_DIR, "reflected_transmitted_energy.png")
plt.tight_layout()
plt.savefig(fname, dpi=200)
plt.close(fig)
print("Saved energy bar chart:", fname)

print("\nAll figures and animation saved into the folder:", FIG_DIR)
